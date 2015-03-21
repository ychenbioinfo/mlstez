import csv
import SeqAlignParallel
import ConsensusSeqs
import HetSearchParallel

from Bio import SeqIO

class Primer(object):
    
    def __init__(self,seqF,seqR,lscore,rscore,id,des):
        self.f = seqF
        self.r = seqR
        self.id = id
        self.des = des
        self.fs = lscore
        self.rs = rscore
        self.fl = len(seqF)
        self.rl = len(seqR)

class ProjectEnviroment(object):
    
    def __init__(self, parameters, msgHandle):
        self.parameters = parameters
        self.msgHandle = msgHandle
        self.Seqs = []
        self.Barcodes = []
        self.Primers = []
        self.AlignedSeqs = []
        self.SortedSeqs = {}
        self.locusLengthRange = {}
        self.consSeqs = {}
        self.HetSeqs = {}
        self.HetInfo = {}
        self.status = 0
        self.ismultirun = 0
        self.locusLengthsInfo = None
        self.num_unbarcode = 0
        self.num_unprimer = 0
        self.StrainStats = None
        self.HetStats = None
        self.SymBarcode = True

    def loadFiles(self):
        isokay, errMsg = self.__readBarcodes()
        if(not isokay):
            return (False, errMsg)
        isokay, errMsg = self.__readPrimers()
        if(not isokay):
            return (False, errMsg)
        isokay, errMsg = self.__readSeqs()
        if(not isokay):
            return (False, errMsg)
        self.status = 1
        return (True, None)

    def numprimers(self):
        return len(self.Primers)
    
    def numbarcodes(self):
        return len(self.Barcodes)
    
    def numseqs(self):
        return len(self.Seqs)
    
    def seqlengths(self):
        seqlens = []
        for seq in self.Seqs:
            seqlens.append(len(seq))
        return seqlens

    def __readSeqs(self):
        files = self.parameters.Seq_Files
        filetype = self.parameters.Filetype
        scoretype = self.parameters.ScoringSys
        
        self.showMsg('Loading sequences...', end="")
        #if(self.msgHandle is not None):
        #    self.msgHandle.showMsg('Loading sequences...', end="")

        Seqs = []
        for file in files:
            try: 
                handle = open(file,'rU')
                if(filetype == "FASTA"):
                    for seq in SeqIO.parse(handle,"fasta"):
                        seq = seq.upper()
                        seqlen = len(seq)
                        quality = [50] * seqlen
                        seq.letter_annotations["phred_quality"] = quality
                        Seqs.append(seq)
                else:
                    if(scoretype == "phred33"):
                        for seq in SeqIO.parse(handle,"fastq-sanger"):
                            seq = seq.upper()
                            Seqs.append(seq)
                    else:
                        for seq in SeqIO.parse(handle,"fastq-solexa"):
                            seq = seq.upper()
                            Seqs.append(seq)
            except Exception as e:
                return (False, e)
        self.showMsg('done!')
        #if(self.msgHandle is not None):                
        #    self.msgHandle.showMsg('done!')
        
        self.Seqs = Seqs
        return (True, None)

    def __readPrimers(self):
        primers = []
        primerfile = self.parameters.Primer_File
        self.showMsg('Loading primers...', end="")
        #if(self.msgHandle is not None):
        #    self.msgHandle.showMsg('Loading primers...', end="")
        try: 
            with open(primerfile,'rU') as csv_primer:
                reader = csv.reader(csv_primer)
                currlocus = {}
                for line in reader:
                    locus = line[0].strip()
                    locus = locus.replace(' ','.')
                    primer_F = line[1].strip()
                    primer_R = line[2].strip()
                    primer_F_score = (len(primer_F) - self.parameters.MaxMisMatch) * self.parameters.MatchScore + \
                    max(self.parameters.GapScore,self.parameters.MismatchScore) * self.parameters.MaxMisMatch
                    primer_R_score = (len(primer_R) - self.parameters.MaxMisMatch) * self.parameters.MatchScore + \
                    max(self.parameters.GapScore,self.parameters.MismatchScore) * self.parameters.MaxMisMatch           
                    primer = Primer(primer_F,primer_R,primer_F_score,primer_R_score,locus,locus)
                    primers.append(primer)
        except Exception as e:
            return (False, e)
        self.showMsg('done!')
        #if(self.msgHandle is not None):
        #    self.msgHandle.showMsg('done!')
        self.Primers = primers 
        return (True, None)
        
    def __readBarcodes(self):
        Barcodes = []
        barcodefile = self.parameters.Barcode_File
        self.showMsg('Loading barcodes...', end="")
        #if(self.msgHandle is not None):
        #    self.msgHandle.showMsg('Loading barcodes...', end="")
        #
        try:
            with open(barcodefile,'rU') as csv_barcode:
                reader = csv.reader(csv_barcode)
                firstline = True
                
                for line in reader:
                    if(firstline): #detect barcode mode
                        numcol = len(line)
                        if(numcol == 2):
                            self.SymBarcode = True
                        elif(numcol == 3):
                            self.SymBarcode = False
                        else:
                            raise ValueError ("Incorrect format of barcode file!")
                    firstline = False
                    
                    species = line[0].strip()
                    species = species.replace(' ','.')
                    barcode = line[1].strip()
                    barcode = barcode.upper()
                    minscore = (len(barcode) - self.parameters.MaxMisMatch) * self.parameters.MatchScore + \
                    max(self.parameters.GapScore,self.parameters.MismatchScore) * self.parameters.MaxMisMatch
                    if (self.SymBarcode):
                        Barcode = Primer(barcode,barcode,minscore,minscore,barcode,species)
                    else:
                        barcodeR = line[2].strip()
                        barcodeR = barcode.upper()
                        Barcode = Primer(barcode,barcodeR,minscore,minscore,barcode,species)
                    Barcodes.append(Barcode)

        except Exception as e:
            return (False, e)
        self.showMsg('done!')
        #if(self.msgHandle is not None):
        #    self.msgHandle.showMsg('done!')
        self.Barcodes = Barcodes
        return (True, None)
    
    def alignSeqs(self):
        self.Aligns = SeqAlignParallel.SeqAlignments(self)
        self.ismultirun = 1
        self.Aligns.Run()
        #self.Aligns.AlignBarcodes()
        #self.Aligns.AlignPrimers()
        self.ismultirun = 0
        self.AlignedSeqs = self.Aligns.alignedseqs
        self.num_unbarcode = self.Aligns.num_unbarcode
        self.num_unprimer = self.Aligns.num_unprimer
        self.status = self.status + (1<<2)
        self.Aligns = None
        return (True, None)
    
    def alignStop(self):
        self.Aligns.Stop()
        self.Aligns = None
    
    def sortAlignedSeqs(self):
        self.showMsg('Sorting the aligned reads by barcodes and primers...', end="")
        #if(self.msgHandle is not None):
        #    self.msgHandle.showMsg('Sorting the aligned reads by barcodes and primers...', end="")
        for seq in self.AlignedSeqs:
            if(seq.strain == ""): continue
            strainname = seq.strain
            genename = seq.gene
            if(genename == ""): genename = "unmapped"
            if(strainname in self.SortedSeqs):
                strainseq = self.SortedSeqs[strainname]
                if(genename in strainseq):
                    tmpseq = strainseq[genename]
                    tmpseq.append(seq)
                else:
                    tmpseq = []
                    tmpseq.append(seq)
                    strainseq[genename] = tmpseq
            else:
                tmpseq = []
                tmpseq.append(seq)
                strainseq = {}
                strainseq[genename] = tmpseq
                self.SortedSeqs[strainname] = strainseq
        self.showMsg('done!')
        #if(self.msgHandle is not None):
        #    self.msgHandle.showMsg('done!')
        self.status = self.status + (1<<3)
        return (True, None)
        
    def locusLengths(self):
        self.showMsg('Calculate length distribution of each locus...', end="")
        locuslens = {}
        for seq in self.AlignedSeqs:
            if(seq.gene == ""):continue
            if(seq.gene in locuslens):
                tmplengths = locuslens[seq.gene]
                tmplengths.append(seq.LocusLength())
            else:
                tmplengths = []
                tmplengths.append(seq.LocusLength())
                locuslens[seq.gene] = tmplengths
        
        locuslenList = []
        locusname = sorted(list(locuslens.keys()))
        for loci in locusname:
            locuslenList.append(locuslens[loci])
        
        self.locusLengthsInfo = {}
        self.locusLengthsInfo['name'] = locusname
        self.locusLengthsInfo['data'] = locuslenList
                
        locusRange = {}
        #lendiver = 0.05
        for gene in locuslens:
            lengths = locuslens[gene]
            totallength = 0
            #for i in lengths:
            #    totallength += i
            #meanlength = totallength/len(lengths)
            #s1 = meanlength * (1 - lendiver)
            #s2 = meanlength * (1 + lendiver)
            sortedlength = sorted(lengths)
            listlen = len(lengths)
            q1 = sortedlength[int(listlen*0.25)]
            q3 = sortedlength[int(listlen*0.75)]
            iqd = q3-q1
            #whisker = 1.5*iqd
            whisker = 2*iqd
            s1 = q1 - whisker
            s2 = q3 +  whisker
            lenRange = {'s1':s1,'s2':s2}
            locusRange[gene] = lenRange
            
        self.locusLengthRange = locusRange
        self.showMsg('done!')
        return (True, None)
    
    def hetStats(self):
        strain_ids = self._UniqueIDs(self.Barcodes,1)
        locus_ids = self._UniqueIDs(self.Primers,2)
        stats = []
        outline = ['Strain']
        outline += locus_ids
        stats.append(list(outline))
        #self.HetSeqs = {}
        #self.HetInfo = {}
        for strain in strain_ids:
            outline = [strain]
            if(strain in self.HetInfo):
                strainHet = self.HetInfo[strain]
                for locus in locus_ids:
                    if(locus in strainHet):
                        if(strainHet[locus] == 1):
                            outline.append("Yes")
                        elif(strainHet[locus] == 0):
                            outline.append("No")
                        else:
                            outline.append("NA")
                    else:
                        outline.append("NA")
            else:
                for locus in locus_ids:
                    outline.append("NA")
            stats.append(list(outline))
        self.HetStats = stats
        return (True, None)
            
    def strainStats(self):
        #sortedSeqs,barcodes,primers,minReadNum,lengthRange
        self.showMsg('Generate statistical information of each sample...', end="")
        strain_ids = self._UniqueIDs(self.Barcodes,1)
        locus_ids = self._UniqueIDs(self.Primers,2)
        locus_ids.append("unmapped")
        locustotal = {}
        stats = []
        outline = ['Strain']
        outline += locus_ids
        outline.append('Sum')
        stats.append(list(outline))
        for strain in strain_ids:
            outline = [strain]
            straintotal = 0
            if(strain in self.SortedSeqs):
                strainAlns = self.SortedSeqs[strain]
                for locus in locus_ids:
                    if(locus in strainAlns):
                        seqs = strainAlns[locus]
                        seqpassed = 0
                        if(locus != 'unmapped'):
                            seqpassed = self._SeqPassed(seqs,self.locusLengthRange[locus])
                            seqcount = str(seqpassed) + '(' + str(len(seqs)) + ')'
                            outline.append(seqcount)
                        else:
                            outline.append(str(len(seqs)))
                        straintotal += len(seqs)
                        
                        if(locus in locustotal):
                            locustotal[locus] += len(seqs)
                        else:
                            locustotal[locus] = len(seqs)
                    else:
                        if(locus != 'unmapped'):
                            outline.append('0(0)')
                        else:
                            outline.append('0')
            else:
                for locus in locus_ids[0:-1]:
                    outline.append('0(0)')
                outline.append('0')
            outline.append(str(straintotal))
            stats.append(list(outline))
            #HTMLTable += _DList2HTMLRow(outline,minReadNum)
            #csvwriter.writerow(outline)
            
        outline = ["Total"]
        totalreads = 0
        totalpass = 0
        for locus in locus_ids:
            if(locus in locustotal):
                outline.append(str(locustotal[locus]))
                totalreads += locustotal[locus]
            else:
                outline.append('0')
        outline.append(totalreads)
        stats.append(list(outline))
        self.StrainStats = stats
        self.showMsg('done!')
        return (True, None)
    
    def GenerateConsensus(self):
        #self.showMsg('Generating consensus sequences...', end="")
        #if(self.msgHandle is not None):
        #    self.msgHandle.showMsg('Generating consensus sequences...', end="")
        consseqs = ConsensusSeqs.ConsensusSeqs(self)
        (isokay, errMsg) = consseqs.makeConsensus()
        if(isokay):    
            self.consSeqs = consseqs.ConsSeqs
            try:
                for gene in self.consSeqs:
                    outfile = self.parameters.Out_Folder + "/" + gene + ".cons.fasta"
                    fh_out = open(outfile,'w')
                    seqs = self.consSeqs[gene]
                    for seq in seqs:
                        fh_out.write(seq.format("fasta").decode('utf-8'))
                    fh_out.close()
                self.status = self.status + (1<<4)
                #self.showMsg('done!')
                #if(self.msgHandle is not None):
                #    self.msgHandle.showMsg('done!')
                return (True, None)
            except Exception as errMsg:
                return (False, errMsg)
        else:
            return (False, errMsg)
    
    def DumpUnmappedReads(self):
        self.showMsg('Dumping unaligned reads to file...', end="")
        #if(self.msgHandle is not None):
            #self.msgHandle.showMsg('Dumping unaligned reads to file...', end="")
            #print ('Dumping unaligned reads to file...', end = "")

        unmapfile = self.parameters.Out_Folder + '/UnmappedReads.seq'
        unmapSeqs = {}
        for seq in self.AlignedSeqs:
            if(seq.barcode == ""):
                if('0barcode' in unmapSeqs):
                    tmpseqs = unmapSeqs['0barcode']
                    tmpseqs.append(seq)
                else:
                    tmpseqs = []
                    tmpseqs.append(seq)
                    unmapSeqs['0barcode'] = tmpseqs
            else:
                if(seq.gene == ""):
                    if(seq.barcode in unmapSeqs):
                        tmpseqs = unmapSeqs[seq.barcode]
                        tmpseqs.append(seq)
                    else:
                        tmpseqs = []
                        tmpseqs.append(seq)
                        unmapSeqs[seq.barcode] = tmpseqs

        barcodes = sorted(list(unmapSeqs.keys()))
        
        try:
            fh_out = open(unmapfile,'w')
            
            for barcode in barcodes:
                if(barcode == '0barcode'):
                    fh_out.write('<NoBarcode>\n')
                else:
                    outline = '<' + barcode + '>\n'
                    fh_out.write(outline)
                seqs = unmapSeqs[barcode]
                for seq in seqs:
                    fh_out.write(seq.seq.format('fasta'))
            fh_out.close()
            self.showMsg('done!')
            #if(self.msgHandle is not None):
            #    self.msgHandle.showMsg('done!')
            return (True, None)
        except Exception as e:
            return (False, e)
    
    def HetSearch(self):
        hetsearch = HetSearchParallel.HetSearch(self)
        (isokay, errMsg) = hetsearch.Run()
        if(isokay):
            try:
                self.HetSeqs = hetsearch.HetSeqs
                self.HetInfo = hetsearch.HetInfo
                outfile = self.parameters.Out_Folder + "/" + "Het.cons.fasta"
                fh_out = open(outfile,'w')
                for strain in self.HetSeqs:
                    strainseqs = self.HetSeqs[strain]
                    for gene in strainseqs:
                        hetseqs = strainseqs[gene]
                        fh_out.write(hetseqs['seq1'].format("fasta").decode('utf-8'))
                        fh_out.write(hetseqs['seq2'].format("fasta").decode('utf-8'))
                fh_out.close()
                return (True, None)
            except Exception as e:
                return (False, errMsg)
        else:
            return (False, errMsg)
    
    def _UniqueIDs(self, refseqs, mode):
    ##mode 1 for barcode, mode 2 for primer
        ids = []
        for refinfo in refseqs:
            if(mode == 1): ids.append(refinfo.des)
            elif(mode == 2): ids.append(refinfo.id)
        uniqueids = sorted(set(ids))
        return list(uniqueids)
    
    def _SeqPassed(self, seqs, lenRange):
        passcount = 0
        for seq in seqs:
            seqlen = seq.LocusLength()
            if(seqlen >= lenRange['s1'] and seqlen <= lenRange['s2']):
                passcount += 1
        return passcount

    def showMsg(self, msg, end="\n"):
        if(self.msgHandle is not None):
            self.msgHandle.emit(msg, end)
        else:
            stderr.write(unicode(msg + end))