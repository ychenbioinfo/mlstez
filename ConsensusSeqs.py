#!/usr/bin/env python


class ConsensusSeqs(object):
    def __init__(self, projenv):
        self.parameters = projenv.parameters
        self.SortedSeqs = projenv.SortedSeqs
        self.locusLengthRange = projenv.locusLengthRange
        self.msgHandle = projenv
        self.ConsSeqs = {}
    
    def makeConsensus(self):
        import tempfile
        from Bio.Align.Applications import MuscleCommandline
        from Bio import AlignIO
        from io import StringIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        import os
        self.msgHandle.showMsg('Generating consensus sequences...', "")
        #stderr.write ('\nGenerating consensus sequences...')
        ConsSeqs = {}
        MUSCLE = self.parameters.MuscleCMD
        try:
            for strain in self.SortedSeqs:
                strainSeqs = self.SortedSeqs[strain]
                for gene in strainSeqs:
                    if gene == "unmapped":
                        continue
                    geneSeqs = strainSeqs[gene]
                    lenRange = self.locusLengthRange[gene]
                    sortedseqs = self._SortSeqs(geneSeqs, lenRange, self.parameters.MinReadNum,
                                                self.parameters.MaxReadNum)
                    if sortedseqs != "":
                        tmpfile = tempfile.NamedTemporaryFile('w', delete=False)
                        tmpname = tmpfile.name
                        for seq in sortedseqs:
                            tmpfile.write(seq.format("fasta"))
                        tmpfile.flush()
                        tmpfile.close()
                        cmdline = MuscleCommandline(MUSCLE, input=tmpname)
                        # print(cmdline)
                        stdout, stderr = cmdline()
                        os.remove(tmpname)
                        align = AlignIO.read(StringIO(stdout), "fasta")
                        #print (align)
                        consensus = self._AlignConsensus(align)
                        #print (consensus)
                        seqrec = SeqRecord(Seq(consensus, generic_dna), id=strain, description=gene)
                        if gene not in ConsSeqs:
                            genecons = []
                            genecons.append(seqrec)
                            ConsSeqs[gene] = genecons
                        else:
                            genecons = ConsSeqs[gene]
                            genecons.append(seqrec)
            self.ConsSeqs = ConsSeqs
            self.msgHandle.showMsg('done!')
            # print("Consensus Done")
            return True, None
        except Exception as e:
            return False, e
    
    def _SortSeqs(self,alnseqs,lenRange,minReadNum,maxReadNum):
    
        if len(alnseqs) < minReadNum : return ""
        
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
            raise ValueError("%s not defined!" %('/'.join(bases)))
        
        return base
