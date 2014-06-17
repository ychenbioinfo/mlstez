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
    
    def __len__(self):
        return len(self.seq)
    
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
        self.msghandle = projenv.msgHandle
        self.alignedseqs = []
    
    def Run(self):
        self.BarcodeSearch()
        self.PrimerSearch()
    
    def BarcodeSearch(self):
    #def BarcodeSearch(paras,seqs,barcodes):
        #stderr.write ('\nSearching for barcodes in reads...\n')
        #self.msghandle.showMsg("Searching for barcodes in reads...")
        seqcount = 0
        unmapcount = 0
        padlen = self.paras.PadLength - self.paras.FlankingLength
        if(padlen < 0):
            padlen = 0
        reglen = 2*self.paras.FlankingLength + self.paras.BarcodeLen
        
        for seq in self.seqs:
            seqcount += 1
            #if(seqcount % 1000 == 0):
            #    stderr.write (unicode(seqcount) + " reads have been processed...\n")
            
            isMatch,alnrec = self._SeqSearch(seq,self.barcodes,padlen,reglen)
    
            if(isMatch):
                alignSeq = AlignedSeq(seq)
                alignSeq.barcode = alnrec.id
                alignSeq.strain = alnrec.des
                alignSeq.alnBarcode = alnrec
                self.alignedseqs.append(alignSeq)
            else:
                alignSeq = AlignedSeq(seq)
                alignSeq.barcode = ''
                alignSeq.strain = ''
                alignSeq.alnBarcode = ''
                self.alignedseqs.append(alignSeq)
                unmapcount += 1

    def PrimerSearch(self):
        #stderr.write ('\nSearching for self.primers in reads...\n')
        seqcount = 0
        unmappcount = 0
        maxprimerlen = self._MaxPrimerLen()
        padlen = self.paras.UniLength - self.paras.FlankingLength
        if(padlen < 0):
            padlen = 0
        reglen = 2*self.paras.FlankingLength + maxprimerlen
        
        for barcodedseq in self.alignedseqs:
            seqcount += 1
            #if(seqcount % 1000 == 0):
                #stderr.write (unicode(seqcount) + " reads have been processed...\n")
                
            if(barcodedseq.barcode == ''): continue
            
            seqs,seqe = barcodedseq.BarcodeFreeRegion()
            seqr = len(barcodedseq.seq) - seqe
            trimed_seq = barcodedseq.TrimBarcode()
    
            isMatch,alnrec = self._SeqSearch(trimed_seq,self.primers,padlen,reglen)
            if(isMatch):
                alnrec = self._AlignAddPad(alnrec,seqs)
                barcodedseq.alnPrimer = alnrec
                barcodedseq.gene = alnrec.id
            else:
                trimed_seq_rv = trimed_seq.reverse_complement()
                #print (trimed_seq_rv.format("fasta"))
                isMatch,alnrec = self._SeqSearch(trimed_seq_rv,self.primers,padlen,reglen)
                if(isMatch):
                    alnrec = self._AlignAddPad(alnrec,seqr)
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
            if(not isMatch): unmappcount += 1
    
        #stderr.write ('Done!\n')
        #stderr.write ("%s reads cannot be aligned to primers\n" %(unmappcount))
        #return self.alignedseqs
        
    def _SeqSearch(self,seq,refseqs,padlen,reglen):
        import SWAlign
        seq_len = len(seq)
        totallen = padlen + reglen
        seq_L = seq[padlen:totallen]
        seq_R = seq[(seq_len - totallen):(seq_len - padlen)]
        seq_Rr = seq_R.reverse_complement()
        trim_len = seq_len - 2*padlen
        
        alignRes = AlignRes(seq)
        sw = SWAlign.LocalAlignment(SWAlign.NucleotideScoringMatrix(self.paras.MatchScore, self.paras.MismatchScore),self.paras.GapScore)
        isMatch = False       
        isMatch,qrec = self._QuickSearch(seq_L,seq_Rr,refseqs,sw,self.paras.MatchScore,trim_len)
        
        if(isMatch):
            qrec = self._AlignAddPad(qrec,padlen)
            alignRes.append(qrec)
        else:
            for refseq in refseqs:
                alnrec = AlignRecord()
                l_align = self._SeqAlign(refseq.f,seq_L,sw)
                r_align = self._SeqAlign(refseq.r,seq_Rr,sw)
                if(l_align['mismatch'] <= self.paras.MaxMisMatch or r_align['mismatch'] <= self.paras.MaxMisMatch):
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
                        alnrec = self._AlignAddPad(alnrec,padlen)
                        alignRes.append(alnrec)
                        
        if(isMatch):
            alnrec = alignRes.BestAlign()
            return (isMatch,alnrec)
        else:
            return (isMatch,"")
        
    def _QuickSearch(self,seqL,seqR,refseqs,sw,matchscore,seqlen):
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
                    align = self._SeqAlign(refseq.r,seqR,sw)
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
                    align = self._SeqAlign(refseq.f,seqL,sw)
                    alignRes.ls = align['s']
                    alignRes.le = align['e']
                    alignRes.lscore = align['score']
                    break
    
        return (Matched,alignRes)    
        
        #stderr.write ('Done!\n')
        #stderr.write ("%s reads cannot be aligned to self.barcodes.\n" %(unmapcount))
        #self.barcodeseq
        
    def _AlignAddPad(self,alnrec,padlen):
        alnrec.ls = alnrec.ls + padlen
        alnrec.le = alnrec.le + padlen
        alnrec.rs = alnrec.rs + padlen
        alnrec.re = alnrec.re + padlen
        return alnrec
    
    
    def _SeqAlign(self,ref,query,sw):
            query = query.seq
            align = sw.align(ref,query)
            start = align.q_pos
            end = align.q_end - 1
            #align.dump()
            return {'s':start,'e':end,'score':align.score,'mismatch':align.mismatches}
    
    def _MaxPrimerLen(self):
        maxlen = 0
        for primer in self.primers:
            if (primer.fl > maxlen): maxlen = primer.fl
            if (primer.rl > maxlen): maxlen = primer.rl
        return maxlen

class ConsensusSeqs(object):
    def __init__(self, parameters, sortedSeqs, lengthRanges):
        self.parameters = parameters
        self.SortedSeqs = sortedSeqs
        self.locusLengthRange = lengthRanges
        self.ConsSeqs = {}
    
    def makeConsensus(self):
        import tempfile
        from Bio.Align.Applications import MuscleCommandline
        from Bio import AlignIO
        from io import StringIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        #stderr.write ('\nGenerating consensus sequences...')
        ConsSeqs = {}
        MUSCLE = self.parameters.MuscleCMD
        try:
            for strain in self.SortedSeqs:
                strainSeqs = self.SortedSeqs[strain]
                for gene in strainSeqs:
                    if(gene == "unmapped"):
                        continue
                    geneSeqs = strainSeqs[gene]
                    lenRange = self.locusLengthRange[gene]
                    sortedseqs = self._SortSeqs(geneSeqs,lenRange,self.parameters.MinReadNum,self.parameters.MaxReadNum)
                    if(sortedseqs != ""):
                        tmpfile = tempfile.NamedTemporaryFile('w')
                        tmpname = tmpfile.name
                        for seq in sortedseqs:
                            tmpfile.write(seq.format("fasta"))
                        tmpfile.flush()
                        cmdline = MuscleCommandline(MUSCLE,input=tmpname)
                        STDOUT, STDERR  = cmdline()
                        align = AlignIO.read(StringIO(STDOUT.decode('utf-8')), "fasta")
                        #print (align)
                        consensus = self._AlignConsensus(align)
                        #print (consensus)
                        tmpfile.close()
                        seqrec = SeqRecord(Seq(consensus,generic_dna),id=strain,description=gene)
                        if(gene not in ConsSeqs):
                            genecons = []
                            genecons.append(seqrec)
                            ConsSeqs[gene] = genecons
                        else:
                            genecons = ConsSeqs[gene]
                            genecons.append(seqrec)
            self.ConsSeqs = ConsSeqs
            return (True, None)
        except Exception as e:
            return (False, e)

    def _SortSeqs(self,alnseqs,lenRange,minReadNum,maxReadNum):

        if(len(alnseqs) < minReadNum):return ""
        
        seqrec = {}
        scorerec = {}
        for alnseq in alnseqs:
            seq = alnseq.TrimPrimer()
            if(len(seq) < lenRange['s1'] or len(seq) > lenRange['s2']) : continue
            scores = seq.letter_annotations["phred_quality"]
            avescore = sum(scores)/len(scores)
            seqrec[seq.id] = seq
            scorerec[seq.id] = avescore
            
        if(len(seqrec) < minReadNum):return ""
        
        sorted_ids = sorted(scorerec, key=scorerec.get,reverse=True)
        sortedseqs = []
        for seqid in sorted_ids:
            sortedseqs.append(seqrec[seqid])
        
        return sortedseqs[:maxReadNum-1]
    
    def _AlignConsensus(self, alignment):
    
        consensus = ''
        con_len = alignment.get_alignment_length()
        gapchar = '-'
        #consuscut = self.parameters.ConsensusCut
        
        for n in xrange(con_len):
            
            base_dict = {}
            num_bases = 0
    
            for record in alignment._records:
    
                if n < len(record.seq):
                    if record.seq[n] not in base_dict:
                        base_dict[record.seq[n]] = 1
                    else:
                        base_dict[record.seq[n]] += 1
                    num_bases = num_bases + 1

            max_bases = []
            max_size = 0

            for base in base_dict:
                if(base_dict[base] > max_size):
                    max_size = base_dict[base]
                    max_bases = [base]
                elif(base_dict[base] == max_size):
                    max_bases.append(base)
            
            if gapchar in max_bases: max_bases.remove(gapchar)
            
            if(len(max_bases) == 1):
                    consensus += max_bases[0]
            elif(len(max_bases) > 1):
                base = self._IUPACambiguity(sorted(max_bases))
                consensus += base
                
        return consensus
    
    def _IUPACambiguity(self, bases):
        base = ''
        if(bases == ['A','G']): base = 'R'
        elif(bases == ['A','C']): base = 'M'
        elif(bases == ['A','T']): base = 'W'
        elif(bases == ['C','T']): base = 'Y'
        elif(bases == ['C','G']): base = 'S'
        elif(bases == ['G','T']): base = 'K'
        elif(bases == ['A','C','G']): base = 'V'
        elif(bases == ['A','C','T']): base = 'H'
        elif(bases == ['A','G','T']): base = 'D'
        elif(bases == ['C','G','T']): base = 'B'
        elif(bases == ['A','T','C','G']): base = 'N'
        else:
            raise ValueError ("%s not defined!" %('/'.join(bases)))
        
        return base