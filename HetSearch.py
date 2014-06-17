#!/usr/bin/env python
from __future__ import division
import tempfile
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from io import StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from scipy.stats import ttest_1samp
import sys

class NucleotideScoringMatrix(object):
    def __init__(self, match=2, mismatch=-2, gapscore=-1):
        self.match = match
        self.mismatch = mismatch
        self.gapscore = gapscore
    def score(self, one, two):
        if one == two:
            return self.match
        if(one == '-' or two == '-'):
            return self.gapscore
        return self.mismatch

class HetSearch(object):
    def __init__(self, projenv):
        self.parameters = projenv.parameters
        self.SortedSeqs = projenv.SortedSeqs
        self.locusLengthRange = projenv.locusLengthRange
        self.msgHandle = projenv

        self.MinReadNum = 5
        self.MinReadRatio = 0.2
        self.MinVariantRatio = 0.3
        self.HeteroPvalue = 0.001
        
        self.HetSeqs = {}
        self.HetInfo = {}

    def Run(self):
        self.msgHandle.showMsg ('Searching for heterozygote...', "")
        #stderr.write ('\nGenerating consensus sequences...')
        
        MUSCLE = self.parameters.MuscleCMD
        #try:
        for strain in self.SortedSeqs:
            strainSeqs = self.SortedSeqs[strain]
            for gene in strainSeqs:
                if(gene == "unmapped"):
                    continue
                geneSeqs = strainSeqs[gene]
                lenRange = self.locusLengthRange[gene]
                filteredseqs = self.__SeqFilter(geneSeqs,lenRange,self.MinReadNum)
                if(filteredseqs != ""):
                    alignSeqs = self.__MuscleAlignment(filteredseqs, MUSCLE)
                    variantBases = self.__getVariants(alignSeqs, self.MinVariantRatio)
                    (isHet, index) = self.__isHetero(variantBases, self.HeteroPvalue, self.MinReadRatio)
                    if(isHet):
                        (seqN1, seqN2) = self.__getSeqGroups(alignSeqs,index)
                        seqs1 = self.__getSeqs(filteredseqs, seqN1)
                        seqs2 = self.__getSeqs(filteredseqs, seqN2)
                        aligns1 = self.__MuscleAlignment(seqs1, MUSCLE)
                        aligns2 = self.__MuscleAlignment(seqs2, MUSCLE)
                        cons1 = self.__AlignConsensus(aligns1)
                        cons2 = self.__AlignConsensus(aligns2)
                        gname1 = gene + "_allele1"
                        gname2 = gene + "_allele2"
                        seqrec1 = SeqRecord(Seq(cons1,generic_dna),id=strain,description=gname1)
                        seqrec2 = SeqRecord(Seq(cons2,generic_dna),id=strain,description=gname2)
                        heteSeqs = {}
                        heteSeqs['seq1'] = seqrec1
                        heteSeqs['seq2'] = seqrec2
                        straininfo = {}
                        if(strain in self.HetSeqs):
                            straininfo = self.HetSeqs[strain]
                        straininfo[gene] = heteSeqs
                        self.HetSeqs[strain] = straininfo
                        hetinfo = {}
                        if(strain in self.HetInfo):
                            hetinfo = self.HetInfo[strain]
                        hetinfo[gene] = 'TRUE'
                        self.HetInfo = hetinfo
                    else:
                        hetinfo = {}
                        if(strain in self.HetInfo):
                            hetinfo = self.HetInfo[strain]
                        hetinfo[gene] = 'FALSE'
                        self.HetInfo = hetinfo
                else:
                    hetinfo = {}
                    if(strain in self.HetInfo):
                        hetinfo = self.HetInfo[strain]
                    hetinfo[gene] = "NA"
                    self.HetInfo = hetinfo
        self.msgHandle.showMsg ('done!')
        return (True, None)
        #except Exception as e:
        #    return (False, e)
    
    def __isHetero(self, variantBases, maxpvalue, minreadratio):
        recnum = len(variantBases)
        if(recnum > 0):
            #print ("variantsite is %s" %variantBases[0])
            #print ("variantlength is %s" %(len(variantBases[0])))
            scores = []
            for i in range(recnum-1):
                seq1 = variantBases[i]
                seq2 = variantBases[i+1]
                score = self.__alignscore(seq1, seq2)
                scores.append(score)
                #print ("seq1 %s - seq2 %s: %s" %(i, i+1, score))
            isunique, scoreinfo = self.__minScore(scores, minreadratio)
            if(not isunique):
                return (False, 0)
            else:
                t_statistic, p_value = ttest_1samp(scoreinfo['scores'], scoreinfo['minscore'])
                #print("minscore: %s\tminindex: %s\tpvalue: %s" %(scoreinfo['minscore'],scoreinfo['minindex'],p_value))
                if(p_value <= maxpvalue):
                    #print ("Is hetero!")
                    return (True, scoreinfo['minindex'])
                return (False, 0)
        
    
    def __MuscleAlignment(self, sortedseqs, MUSCLE):
        tmpfile = tempfile.NamedTemporaryFile('w')
        tmpname = tmpfile.name
        for seq in sortedseqs:
            tmpfile.write(seq.format("fasta"))
        tmpfile.flush()
        cmdline = MuscleCommandline(MUSCLE,input=tmpname)
        STDOUT, STDERR  = cmdline()
        align = AlignIO.read(StringIO(STDOUT.decode('utf-8')), "fasta")
        tmpfile.close()
        return align
    
    def __SeqFilter(self,alnseqs,lenRange,minReadNum):
    
        if(len(alnseqs) < minReadNum):return ""
        
        seqs = []
        for alnseq in alnseqs:
            seq = alnseq.TrimPrimer()
            if(len(seq) < lenRange['s1'] or len(seq) > lenRange['s2']) : continue
            seqs.append(seq)
            
        if(len(seqs) < minReadNum):return ""
        
        return seqs

    def __getVariants(self, Aligns, ratio = 0.3):
        """Search variant nuclotide in the alignment"""
        alignlen = Aligns.get_alignment_length()
        lastvariant = 0
        allvariant = []
        variantBase = {}
        recnum = len(Aligns)
        for n in range(alignlen):
            base_dict = {}
            for i in range(recnum):
                curbase = Aligns[i].seq[n]
                upbase = curbase.upper()
                if(curbase not in base_dict):
                    base_dict[curbase] = 1
                else:
                    base_dict[curbase] += 1

            filterbases = self.__baseFilter(base_dict, recnum, ratio)
            if(len(filterbases) > 1):
                for j in range(recnum):
                    curbase = Aligns[j].seq[n]
                    upbase = curbase.upper()
                    if(j in variantBase):
                        variantBase[j] += upbase
                    else:
                        variantBase[j] = upbase
    
        return variantBase
    
    def __baseFilter(self, basedict, seqnum, minratio):
        filterdict = {}
        for base in basedict:
            ratio = basedict[base]/seqnum
            if(ratio < minratio):
                continue
            filterdict[base] = basedict[base]
        return filterdict

    def __minScore(self, scores, ratio = 0.3):
        numscores = len(scores)
        startindex = int(numscores*ratio)
        endindex = numscores - startindex
        minscore = scores[startindex]
        minindex = startindex
        scoredict = {}
        for i in range(startindex,endindex + 1):
            if(scores[i] not in scoredict):
                scoredict[scores[i]] = 1
            else:
                scoredict[scores[i]] += 1
            if(minscore > scores[i]):
                minindex = i
                minscore = scores[i]
        scores.pop(minindex)
        scoreinfo = {}
        isUnique = True
        if(scoredict[minscore] > 1):
            isUnique = False
        scoreinfo['minscore'] = minscore
        scoreinfo['minindex'] = minindex
        scoreinfo['scores'] = scores
        return(isUnique,scoreinfo)

    def __getSeqGroups(self, alignSeqs, index):
        seqcount = 0
        groupA = []
        groupB = []
        for i in range(0,index+1):
            groupA.append(alignSeqs[i].id)
        for j in range(index+1,len(alignSeqs)):
            groupB.append(alignSeqs[j].id)
        return (groupA, groupB)
    
    def __getSeqs(self, seqs, seqnames):
        selSeqs = []
        for i in range(len(seqs)):
            if(seqs[i].id in seqnames):
                selSeqs.append(seqs[i])
        return selSeqs
    
    def __alignscore(self, seq1, seq2, match = 2, mismatch=-2, gap=-1):
        scorematrix = NucleotideScoringMatrix(match,mismatch,gap)
        #(seq1t, seq2t) = trimgap(seq1, seq2)
        score = 0
        seqlen = len(seq1)
        for i in range(seqlen):
            matchscore = scorematrix.score(seq1[i],seq2[i])
            score += matchscore
        return score
    
    def __AlignConsensus(self, alignment):
        consensus = ''
        con_len = alignment.get_alignment_length()
        gapchar = '-'
        #consuscut = self.parameters.ConsensusCut
        
        for n in range(con_len):
            
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
                base = self.__IUPACambiguity(sorted(max_bases))
                consensus += base
                
        return consensus
    
    def __IUPACambiguity(self, bases):
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