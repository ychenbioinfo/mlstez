from multiprocessing import Process, Queue, Manager
from time import time, sleep
from sys import stderr

class AlignRecord(object):
    """class to store single align result"""
    def __init__(self):
        self.id = ''
        self.des = ''
        self.ls = 0
        self.le = 0
        self.lscore = 0
        self.rs = 0
        self.re = 0
        self.rscore = 0
        self.dir = ''

class AlignedSeq(object):
    """class to store multiple align results for one seq"""
    def __init__(self,seq):
        self.seq = seq
        self.seqid = seq.id
        self.barcode = ""
        self.strain = ""
        self.gene = ""
        self.alnBarcode = ""
        self.alnPrimer = ""
    
    def TrimBarcode(self):
        if(self.alnBarcode == ""):
            raise ValueError ("Read %s has not been aligned with barcode" %(self.seqid))
        else:
            seq_s = self.alnBarcode.le + 1
            seq_e = self.alnBarcode.rs
            trimseq = self.seq[seq_s:seq_e]
            if(self.alnBarcode.dir == "-"):
                trimseq = trimseq.reverse_complement(id=True,name=True,description=True)
            return trimseq
    
    def TrimPrimer(self):
        if(self.alnPrimer == ""):
            raise ValueError ("Read %s has not been aligned with primer" %(self.seqid))
        else:
            seq_s = self.alnPrimer.le + 1
            seq_e = self.alnPrimer.rs + 1
            trimseq = self.seq[seq_s:seq_e]
            if(self.alnPrimer.dir == "-"):
                trimseq = trimseq.reverse_complement(id=True,name=True,description=True)
            return trimseq
    
    def BarcodeFreeRegion(self):
        seq_s = self.alnBarcode.le + 1
        seq_e = self.alnBarcode.rs
        return (seq_s,seq_e)
    
    def LocusLength(self):
        seq_s = self.alnPrimer.le + 1
        seq_e = self.alnPrimer.rs + 1
        length = seq_e - seq_s
        return (length)


class AlignRes(object):
    """class to store several align results for single sequence"""
    def __init__(self,seq):
        self.seq = seq
        self.aligns = []
        
    def __len__(self):
        return len(self.aligns)
    
    def __getitem__(self, index):
        return self.aligns[index]
    
    def append(self,alignRes):
        self.aligns.append(alignRes)
        
    def __best(self):
        if(len(self.aligns) < 1):
            raise ValueError ("Error: AlignRes does not contain any record")
        else:
            best = 0
            bestscore = self.aligns[0].lscore + self.aligns[0].rscore
            for i in xrange(0,len(self.aligns)):
                score = self.aligns[i].lscore + self.aligns[i].rscore
            if(score > bestscore):
                best = i
                bestscore = score
        return i
    
    def BestAlign(self):
        index = self.__best()
        alnRec = self.aligns[index]
        return alnRec

class SeqAlignments(object):
    """class to store all aligned sequences"""
    def __init__(self, projenv):
        self.alignedseqs = []
        self.paras = projenv.parameters
        self.seqs = projenv.Seqs
        self.barcodes = projenv.Barcodes
        self.primers = projenv.Primers
        self.msgHandle = projenv
        self.threadNum = self.paras.Threads
        self.symbarcode = projenv.SymBarcode
    
    
    def Run(self):
        self.msgHandle.showMsg("Searching for barcodes in reads...")
        groupnum = int(len(self.seqs)/self.threadNum)
        self.seqgroups = []
        for i in range(self.threadNum):
            startnum = groupnum * i
            endnum = groupnum * (i+1)
            if(i == self.threadNum - 1):
                endnum = len(self.seqs)
            self.seqgroups.append(self.seqs[startnum:endnum])
        
        self.workers = []
        manager = Manager()
        alignedseqs = manager.list()
        primeredseqs = manager.list()
        stats = manager.list()
        unbarcode = manager.list()
        unprimer = manager.list()
    
        for i in range(self.threadNum):
            child = Process(target=BarcodeSearch,
                            args=(self.paras, self.seqgroups[i],
                                  self.barcodes, self.symbarcode, stats, alignedseqs, unbarcode))
            child.start()
            self.workers.append(child)
        
        totalcount = 0
        while any(i.is_alive() for i in self.workers):
            sleep(0.1)
            while len(stats) > 0:
                totalcount += 100
                if(totalcount % 1000 == 0):
                    self.msgHandle.showMsg("%s reads have been processed..."
                                           %(totalcount))
                stats.pop()
        self.msgHandle.showMsg('Done!')
    
        for self.worker in self.workers:
            self.worker.join()
        
        self.workers = []
        self.msgHandle.showMsg("Searching for primers in reads...")
    
        for i in range(self.threadNum):
            child = Process(target=PrimerSearch,
                            args=(self.paras, alignedseqs[i],
                                  self.primers, stats, unprimer,primeredseqs))
            child.start()
            self.workers.append(child)
        
        totalcount = 0
        
        while any(i.is_alive() for i in self.workers):
            sleep(0.1)
            while len(stats) > 0:
                totalcount += 100
                if(totalcount % 1000 == 0):
                    self.msgHandle.showMsg("%s reads have been processed..."
                                           %(totalcount))
                stats.pop()
        
        for self.worker in self.workers:
            self.worker.join()
        
        self.alignedseqs = [aligned for aligns in primeredseqs for aligned in aligns]
        self.num_unbarcode = sum(unbarcode)
        self.num_unprimer = sum(unprimer)
        self.msgHandle.showMsg('Done!')
    
    def Stop(self):
        for worker in self.workers:
            if(worker.is_alive()):
                worker.terminate()
            worker.join()
        return 
    
#def BarcodeSearch(self):
def BarcodeSearch(paras, seqs, barcodes, symbarcode, stats, result, unbarcode):
    
    seqcount = 0
    unmapcount = 0
    alignedseqs = []
    
    padlen = paras.PadLength - paras.FlankingLength
    if(padlen < 0):
        padlen = 0
    reglen = 2*paras.FlankingLength + paras.BarcodeLen
    
    for seq in seqs:
        seqcount += 1
        if(seqcount % 100 == 0):
            stats.append(1)
        
        isMatch,alnrec = _SeqSearch(paras, seq, barcodes, padlen, reglen)

        if(isMatch):
            alignSeq = AlignedSeq(seq)
            alignSeq.barcode = alnrec.id
            alignSeq.strain = alnrec.des
            alignSeq.alnBarcode = alnrec
            alignedseqs.append(alignSeq)
        else:
            if(symbarcode):
                alignSeq = AlignedSeq(seq)
                alignSeq.barcode = ''
                alignSeq.strain = ''
                alignSeq.alnBarcode = ''
                alignedseqs.append(alignSeq)
                unmapcount += 1
            else:
                seq_rev = seq.reverse_complement()
                isMatch,alnrec = _SeqSearch(paras, seq_rev, barcodes, padlen, reglen)
                if(isMatch):
                    alignSeq.barcode = alnrec.id
                    alignSeq.strain = alnrec.des
                    alignSeq.alnBarcode = alnrec
                    alignedseqs.append(alignSeq)
                else:
                    alignSeq = AlignedSeq(seq)
                    alignSeq.barcode = ''
                    alignSeq.strain = ''
                    alignSeq.alnBarcode = ''
                    alignedseqs.append(alignSeq)
                    unmapcount += 1

    result.append(alignedseqs)
    unbarcode.append(unmapcount)
    
    
def PrimerSearch(paras, alignedseqs, primers, stats, unprimer, primeredseqs):
    #stderr.write ('\nSearching for self.primers in reads...\n')
    seqcount = 0
    unmappcount = 0
    maxprimerlen = _MaxPrimerLen(primers)
    padlen = paras.UniLength - paras.FlankingLength
    if(padlen < 0):
        padlen = 0
    reglen = 2*paras.FlankingLength + maxprimerlen
    primered = []
    
    for barcodedseq in alignedseqs:
        seqcount += 1
        if(seqcount % 100 == 0):
            stats.append(1)
            #stderr.write (unicode(seqcount) + " reads have been processed...\n")
            
        if(barcodedseq.barcode == ''): continue
        
        seqs,seqe = barcodedseq.BarcodeFreeRegion()
        seqr = len(barcodedseq.seq) - seqe
        trimed_seq = barcodedseq.TrimBarcode()

        isMatch,alnrec = _SeqSearch(paras,trimed_seq,primers,padlen,reglen)
        if(isMatch):
            alnrec = _AlignAddPad(alnrec,seqs)
            barcodedseq.alnPrimer = alnrec
            barcodedseq.gene = alnrec.id
        else:
            trimed_seq_rv = trimed_seq.reverse_complement()
            #print (trimed_seq_rv.format("fasta"))
            isMatch,alnrec = _SeqSearch(paras,trimed_seq_rv,primers,padlen,reglen)
            if(isMatch):
                alnrec = _AlignAddPad(alnrec,seqr)
                seq_len = len(barcodedseq.seq) - 1
                ls = seq_len - alnrec.re
                le = seq_len - alnrec.rs
                rs = seq_len - alnrec.le
                re = seq_len - alnrec.ls
                alnrec.lscore, alnrec.rscore = alnrec.rscore,alnrec.lscore
                alnrec.ls = ls
                alnrec.le = le
                alnrec.rs = rs
                alnrec.re = re
                alnrec.dir = "-"
                barcodedseq.gene = alnrec.id
                barcodedseq.alnPrimer = alnrec
        primered.append(barcodedseq)
        if(not isMatch): unmappcount += 1
    primeredseqs.append(primered)
    unprimer.append(unmappcount)
    
def _SeqSearch(paras,seq,refseqs,padlen,reglen):
    import SWAlign
    seq_len = len(seq)
    totallen = padlen + reglen
    seq_L = seq[padlen:totallen]
    seq_R = seq[(seq_len - totallen):(seq_len - padlen)]
    seq_Rr = seq_R.reverse_complement()
    trim_len = seq_len - 2*padlen
    
    alignRes = AlignRes(seq)
    sw = SWAlign.LocalAlignment(SWAlign.NucleotideScoringMatrix
                                (paras.MatchScore, paras.MismatchScore), paras.GapScore)
    isMatch = False       
    #isMatch,qrec = _QuickSearch(seq_L,seq_Rr,refseqs,sw,paras.MatchScore,trim_len)
    
    #if(isMatch):
        #qrec = _AlignAddPad(qrec,padlen)
        #alignRes.append(qrec)
    #else:
    #remove quicksearch function
    
    for refseq in refseqs:
        alnrec = AlignRecord()
        l_align = _SeqAlign(refseq.f,seq_L,sw)
        r_align = _SeqAlign(refseq.r,seq_Rr,sw)
        if(l_align['mismatch'] <= paras.MaxMisMatch or r_align['mismatch'] <= paras.MaxMisMatch):
            #if((l_align['score'] + r_align['score'])/2 >= self.paras.MinBarcodeScore):
            if(l_align['score'] >= refseq.fs or r_align['score'] >= refseq.rs):
                #print ("###",barcode)
                isMatch = True
                r_s = seq_len - r_align['e'] - 1
                r_e = seq_len - r_align['s'] - 1
                alnrec.id = refseq.id
                alnrec.des = refseq.des
                alnrec.ls = l_align['s']
                alnrec.le = l_align['e']
                alnrec.lscore = l_align['score']
                alnrec.rs = r_s
                alnrec.re = r_e
                alnrec.rscore = r_align['score']
                alnrec.dir = '+'
                alnrec = _AlignAddPad(alnrec,padlen)
                alignRes.append(alnrec)
                    
    if(isMatch):
        alnrec = alignRes.BestAlign()
        return (isMatch,alnrec)
    else:
        return (isMatch,"")
    
def _QuickSearch(seqL,seqR,refseqs,sw,matchscore,seqlen):
    alignRes = AlignRecord()
    Matched = False
    lseq = seqL.seq.upper()
    rseq = seqR.seq.upper()
    for refseq in refseqs:
        posL = lseq.find(refseq.f)
        posR = rseq.find(refseq.r)
        if(posL >= 0):
            Matched = True
            alignRes.id = refseq.id
            alignRes.des = refseq.des
            alignRes.ls = posL
            alignRes.le = posL + refseq.fl -1
            alignRes.lscore = refseq.fl * matchscore
            alignRes.dir = '+'
            if(posR >= 0):
                rs = posR
                re = posR + refseq.rl -1
                alignRes.rs = seqlen - re -1
                alignRes.re = seqlen - rs -1
                alignRes.rscore = refseq.rl * matchscore
            else:
                align = _SeqAlign(refseq.r,seqR,sw)
                alignRes.rs = seqlen - align['e'] -1
                alignRes.re = seqlen - align['s'] -1
                alignRes.rscore = align['score']
            break
        else:
            if(posR >= 0):
                alignRes.id = refseq.id
                alignRes.des = refseq.des
                Matched = True
                alignRes.dir = '+'
                rs = posR
                re = posR + refseq.rl -1
                alignRes.rs = seqlen - re -1
                alignRes.re = seqlen - rs -1
                alignRes.rscore = refseq.rl * matchscore
                align = _SeqAlign(refseq.f,seqL,sw)
                alignRes.ls = align['s']
                alignRes.le = align['e']
                alignRes.lscore = align['score']
                break

    return (Matched,alignRes)    
    
def _AlignAddPad(alnrec,padlen):
    alnrec.ls = alnrec.ls + padlen
    alnrec.le = alnrec.le + padlen
    alnrec.rs = alnrec.rs + padlen
    alnrec.re = alnrec.re + padlen
    return alnrec


def _SeqAlign(ref,query,sw):
        query = query.seq
        align = sw.align(ref,query)
        start = align.q_pos
        end = align.q_end - 1
        #align.dump()
        return {'s':start,'e':end,'score':align.score,'mismatch':align.mismatches}

def _MaxPrimerLen(primers):
    maxlen = 0
    for primer in primers:
        if (primer.fl > maxlen): maxlen = primer.fl
        if (primer.rl > maxlen): maxlen = primer.rl
    return maxlen
