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
from multiprocessing import Pool
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

def __isHetero(variantBases, maxpvalue, minreadratio):
    recnum = len(variantBases)
    if(recnum > 0):
        #print ("variantsite is %s" %variantBases[0])
        #print ("variantlength is %s" %(len(variantBases[0])))
        scores = []
        for i in range(recnum-1):
            seq1 = variantBases[i]
            seq2 = variantBases[i+1]
            score = __alignscore(seq1, seq2)
            scores.append(score)
            print ("seq1 %s - seq2 %s: %s" %(i, i+1, score))
        isunique, scoreinfo = __minScore(scores, minreadratio)
        if(isunique):
            t_statistic, p_value = ttest_1samp(scoreinfo['scores'], scoreinfo['minscore'])
            print("minscore: %s\tminindex: %s\tpvalue: %s" %(scoreinfo['minscore'],scoreinfo['minindex'],p_value))
            if(p_value <= maxpvalue):
                print ("Is hetero!")
                return (True, scoreinfo['minindex'])
    return (False, 0)

def __MuscleAlignment(sortedseqs, MUSCLE):
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

def __SeqFilter(alnseqs,lenRange,minReadNum):

    if(len(alnseqs) < minReadNum):return ""
    
    seqs = []
    for alnseq in alnseqs:
        seq = alnseq.TrimPrimer()
        if(len(seq) < lenRange['s1'] or len(seq) > lenRange['s2']) : continue
        seqs.append(seq)
        
    if(len(seqs) < minReadNum):return ""
    
    return seqs

def __isConsIden(seq1,seq2,MUSCLE,MinHetVariants):
    NuCoding = ['A','T','C','G']
    conSeqs = []
    conSeqs.append(seq1)
    conSeqs.append(seq2)
    align = __MuscleAlignment(conSeqs, MUSCLE)
    trimSeq1, trimSeq2 = __trimGap(align[0].seq, align[1].seq)
    if(str(trimSeq1) == str(trimSeq2)):
        return True
    vCount = 0
    for i in range(len(trimSeq1)):
        if(trimSeq1[i] != trimSeq2[i]):
            if(trimSeq1[i] != '-' and trimSeq2[i] != '-'):
                if((trimSeq1[i] in NuCoding) and (trimSeq2[i] in NuCoding)):
                    vCount += 1
                elif((trimSeq1[i] in NuCoding) and (trimSeq2[i] not in NuCoding)):
                    bases = __IUPARevCambiguity(trimSeq2[i])
                    if(trimSeq1[i] not in bases):
                        vCount += 1
                elif((trimSeq1[i] not in NuCoding) and (trimSeq2[i] in NuCoding)):
                    bases = __IUPARevCambiguity(trimSeq1[i])
                    if(trimSeq2[i] not in bases):
                        vCount += 1
    
    vCount += __gapcount(trimSeq1)
    vCount += __gapcount(trimSeq2)
    
    if(vCount <= MinHetVariants):
        return True
    else:
        return False

def __gapcount(seq):
    gaps = re.findall(r'-(-+)',str(seq))
    gapcount = 0
    for i in range(len(gaps)):
        gaplen = len(gaps[i]) + 1
        gapcount += gaplen
    return gapcount

def __getVariants(Aligns, ratio = 0.3):
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

        filterbases = __baseFilter(base_dict, recnum, ratio)
        if(len(filterbases) > 1):
            for j in range(recnum):
                curbase = Aligns[j].seq[n]
                upbase = curbase.upper()
                if(j in variantBase):
                    variantBase[j] += upbase
                else:
                    variantBase[j] = upbase

    return variantBase

def __baseFilter(basedict, seqnum, minratio):
    filterdict = {}
    for base in basedict:
        ratio = basedict[base]/seqnum
        if(ratio < minratio):
            continue
        filterdict[base] = basedict[base]
    return filterdict

def __minScore(scores, ratio = 0.3):
    numscores = len(scores)
    startindex = int(numscores*ratio)
    endindex = numscores - startindex
    if(endindex == numscores):
        endindex = numscores - 1
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

def __getSeqGroups(alignSeqs, index):
    seqcount = 0
    groupA = []
    groupB = []
    for i in range(0,index+1):
        groupA.append(alignSeqs[i].id)
    for j in range(index+1,len(alignSeqs)):
        groupB.append(alignSeqs[j].id)
    return (groupA, groupB)

def __getSeqs(seqs, seqnames):
    selSeqs = []
    for i in range(len(seqs)):
        if(seqs[i].id in seqnames):
            selSeqs.append(seqs[i])
    return selSeqs

def __alignscore(seq1, seq2, match = 2, mismatch=-2, gap=-1):
    scorematrix = NucleotideScoringMatrix(match,mismatch,gap)
    #(seq1t, seq2t) = trimgap(seq1, seq2)
    score = 0
    seqlen = len(seq1)
    for i in range(seqlen):
        matchscore = scorematrix.score(seq1[i],seq2[i])
        score += matchscore
    return score

def __AlignConsensus(alignment):
    consensus = ''
    con_len = alignment.get_alignment_length()
    gapchar = '-'
    #consuscut = parameters.ConsensusCut
    
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
            base = __IUPACambiguity(sorted(max_bases))
            consensus += base
            
    return consensus

def __trimGap(seq1, seq2):
    (seq1_s, seq1_e) = __gapRange(seq1)
    (seq2_s, seq2_e) = __gapRange(seq2)
    startpos = 0
    endpos = 0
    if(seq2_s > seq1_s):
        startpos = seq2_s
    else:
        startpos = seq1_s
    
    if(seq2_e < seq1_e):
        endpos = seq2_e
    else:
        endpos = seq1_e

    trimmed1 = seq1[startpos:endpos]
    trimmed2 = seq2[startpos:endpos]
    return (trimmed1, trimmed2)

def __gapRange(seq):
    startpos = 0
    endpos = 0
    seqlen = len(seq)
    i = 0
    while(seq[i] == '-'):
        startpos += 1
        i += 1
    revseq = seq.reverse_complement()
    i = 0
    pos = 0
    while(revseq[i] == '-'):
        pos += 1
        i += 1
    endpos = seqlen - pos
    return (startpos, endpos)

def __IUPACambiguity(bases):
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

def __IUPARevCambiguity(base):
    bases = []
    if(base == 'R'): bases = ['A','G']
    elif(base == 'M'): bases = ['A','C']
    elif(base == 'W'): bases = ['A','T']
    elif(base == 'Y'): bases = ['C','T']
    elif(base == 'S'): bases = ['C','G'] 
    elif(base == 'K'): bases = ['G','T']
    elif(base == 'V'): bases = ['A','C','G']
    elif(base == 'H'): bases = ['A','C','T'] 
    elif(base == 'D'): bases = ['A','G','T']
    elif(base == 'B'): bases = ['C','G','T']
    elif(base == 'N'): bases = ['A','T','C','G']
    else:
        raise ValueError ("%s not defined!" %('/'.join(bases)))
    return bases

if(__name__ == '__main__'):
    infile = sys.argv[1]
    MinReadNum = 5
    MinReadRatio = 0.2
    MinVariantRatio = 0.2
    HeteroPvalue = 0.001
    alignSeqs = AlignIO.read(open(infile, 'rU'), "fasta")
    variantBases = __getVariants(alignSeqs, MinVariantRatio)
    (isHet, index) = __isHetero(variantBases, HeteroPvalue, MinReadRatio)
    if(isHet):
        print ("Is Hetero!")
    #    (seqN1, seqN2) = __getSeqGroups(alignSeqs,index)
    #    seqs1 = __getSeqs(filteredseqs, seqN1)
    #    seqs2 = __getSeqs(filteredseqs, seqN2)
    #    aligns1 = __MuscleAlignment(seqs1, MUSCLE)
    #    aligns2 = __MuscleAlignment(seqs2, MUSCLE)
    #    cons1 = __AlignConsensus(aligns1)
    #    cons2 = __AlignConsensus(aligns2)
    #    gname1 = gene + "_allele1"
    #    gname2 = gene + "_allele2"
    #    seqrec1 = SeqRecord(Seq(cons1,generic_dna),id=strain,description=gname1)
    #    seqrec2 = SeqRecord(Seq(cons2,generic_dna),id=strain,description=gname2)
    #    isIdent = __isConsIden(seqrec1,seqrec2,MUSCLE)
    #    if(isIdent):
    #        print ("strain:%s gene:%s not Hetero!" %(strain, gene))
    #        #hetinfo.het = 0
    #    else:
    #        print ("strain:%s gene:%s is Hetero!" %(strain, gene))
    #        #hetinfo.het = 1
    #        #hetinfo.hetSeqs['seq1'] = seqrec1
    #        #hetinfo.hetSeqs['seq2'] = seqrec2
    #    #return (2,hetSeqs)
    else:
        print ("Not Hetero!")
    #    print ("strain:%s gene:%s not Hetero!" %(strain, gene))
    #    #hetinfo.het = 0
    #    #return (1,None