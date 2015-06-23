<div id="table-of-contents">
<h2>Table of Contents</h2>
<div id="text-table-of-contents">
<ul>
<li><a href="#sec-1">1. Introduction</a></li>
<li><a href="#sec-2">2. System requirements</a></li>
<li><a href="#sec-3">3. Create a new project</a>
<ul>
<li><a href="#sec-3-1">3.1. Sequencing Files (FASTA format; FASTQ format)</a>
<ul>
<li><a href="#sec-3-1-1">3.1.1. FASTA Format</a></li>
<li><a href="#sec-3-1-2">3.1.2. FASTQ Format</a></li>
</ul>
</li>
<li><a href="#sec-3-2">3.2. Barcode File</a>
<ul>
<li><a href="#sec-3-2-1">3.2.1. Example of barcode file for MLSTEZ 2.0</a></li>
<li><a href="#sec-3-2-2">3.2.2. Example of barcode file MLSTEZ 1.0</a></li>
</ul>
</li>
<li><a href="#sec-3-3">3.3. Primer File</a>
<ul>
<li><a href="#sec-3-3-1">3.3.1. Example of primer file</a></li>
</ul>
</li>
<li><a href="#sec-3-4">3.4. Output Folder</a></li>
<li><a href="#sec-3-5">3.5. Advanced parameters</a></li>
</ul>
</li>
<li><a href="#sec-4">4. Run project</a>
<ul>
<li><a href="#sec-4-1">4.1. Barcode and primer identification</a></li>
<li><a href="#sec-4-2">4.2. Generate consensus sequences</a></li>
<li><a href="#sec-4-3">4.3. Dump unmapped reads</a></li>
<li><a href="#sec-4-4">4.4. Search for heterozygous locus</a></li>
</ul>
</li>
<li><a href="#sec-5">5. Open Project</a></li>
<li><a href="#sec-6">6. Merge Projects</a></li>
</ul>
</div>
</div>


# Introduction<a id="sec-1" name="sec-1"></a>

Efficient methods for estimating genetic diversity among the microorganisms are essential for understanding evolutionary history, geographic distribution, pathogenicity and virulence. In the past decades, numerous methods have been developed for typing of bacteria and fungi. Multilocus sequence typing (MLST) based DNA sequencing results, which can be easily archived and shared among different laboratories. MLST is one of the most reliable and informative method for molecular genotyping, and it has been adopted in many bacterial and fungal studies. 
MLSTEasy was designed for next generation sequencing technology (PacBio CCS or Roche 454 platform) based MSLT methods. MLST-Easy, can automatically identify the barcodes and primers used in the PCR reaction, corrects sequencing errors, generates the MLST profile for each isolate, predicts the potential heterozygous locus, and outputs different alleles. 

# System requirements<a id="sec-2" name="sec-2"></a>

MLSTEasy was written in Python, version 2.7.6. The graphic user interface (GUI) was created by PyQt4 (<http://www.riverbankcomputing.com/software/pyqt/download>) and Qt Designer (<http://qt-project.org/doc/qt-4.8/designer-manual.html>). Mac version were tested under Mac OS X 10.9, and Windows version were tested under Windows Vista and Windows 7. The software runs on IBM-compatible PC under 32/64-bit Windows, and Mac OS X 10.6+. The minimum hardware requirements for the program are:
  a processor based on the Intel Pentium 4/AMD Athlon
  200 MB of RAM memory
  hard drive with more than 200 MB available space
  Windows XP or later version, Mac OS X 10.6+ with MUSCLE installed (<http://www.drive5.com/muscle/downloads.htm>)

The recommended hardware requirements for the program are:
  a multiple core processor based on the Intel Core 2 Due/AMD Athlon II (or higher)
  4 GB of RAM memory
  hard drive with more than 1 GB available space
  Windows XP or later version, Mac OS X 10.6+ with MUSCLE installed (<http://www.drive5.com/muscle/downloads.htm>)

# Create a new project<a id="sec-3" name="sec-3"></a>

A new project can be created by clicking the "New Project" button on toolbar or by selecting "Project->New Project" on menu. Creating a new project requires "Sequencing Files", "Barcode File", "Primer File" and properly set of the "Parameters".

## Sequencing Files (FASTA format; FASTQ format)<a id="sec-3-1" name="sec-3-1"></a>

MLSTEZ supports FASTA and FASTQ data file formats. Corresponding data format needs to be selected in the "Advanced" parameters. When you use FASTQ file as input, the corresponding scoring system (phred33 or phred64) needs to be selected in the "Advanced" parameters. Users can obtain the scoring system using FastQC (<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>) or asking the information from sequencing facility. Mutilple files with the same file format and scoring system can be used at one time.

### FASTA Format<a id="sec-3-1-1" name="sec-3-1-1"></a>

FASTA file format must begin with the symbol '>' in the first line of the file; the sequence name is the first word after that symbol. Additional characters in this line are considered to be comments. The sequence data starts in the second line. Nucleotide data can be written in one or more lines. For more detail informaiton of FASTA Format please check <http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml> for more details. 

1.  Example of FASTA format

        >m140505
        AGACTGGACCCACAGCGGGCGAGAGAAGTACAAGCCCCCGCGACTGGAGTAATTCTTAGT
        GATCATTAATCTTTTCTAGACTTTGCTTGACTGAGCTTGACTCAACTTAAAACGTTTGCT
        TGACCAGCCTATTAGAGCCACCGTCAGGTCGGGTCAACAACTATTCAAAGTTTGATTTGC
        CCATCCCCTCTTTGACTATGCTATAAGCACACCCACTGCATACACTTGGCAGCCCCCCCC
        >m140507
        AAGACTGGACCCACAGCGGGCGAGAGAAGTTACAAGCCCCCCGCGACTGGAGTAATTCTT
        AGTGATCATTAATCTTTTCTAGACTTTGCTTGACTGAGCTTGACTCAACTTAAAACGTTT
        GCTTGACCAGCCTATTAGAGCCACCGTCAGGTCGGGTCAACAACTATTCAAAGTTTGATT
        GCCCATCCCCTCTTGACTATGCTATAAGCACACCCACTGCATACACTTGGCAGCCCCCCT
        >m140510
        AGACTGGACCCACAGCGGGCGAGAGAAGTTACAGGCCCCCCGCGACTGGAGTAATCCTTA
        TGATCATTAATCTTTTCTAGACTTTGCTTGACTGAGCTTGACTCAACTTAAAACGTTTG
        CTTGACCAGCCTATTAGAGCCACCGTCAGGTCGGGTCAACAACTATTCAAAGTTTGATTG
        CCCATCCCCTCTTGACTATGCTATAAGCACACCCACTGCATACACTTGGCAGCCCCCCCT
        CTCACCATCCATACCGCATTTACCCATTTTTCATTCCGGCTCACTACCACTATCAAAGTC
        CCCCACGACTGGAAAAGTAACAAATACTAGTACTATAATACTAATACTAATTGACTTGCT

### FASTQ Format<a id="sec-3-1-2" name="sec-3-1-2"></a>

A FASTQ file normally uses four lines per sequence. The first line begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line); the second line is the rw sequence letters; the third line begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again; and the last line encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence. 

1.  Example of FASTQ format

        @NCYC361-11a03.q1k bases 1 to 1576
        GCGTGCCCGAAAAAATGCTTTTGGAGCCGCGCGTGAAAT...
        +NCYC361-11a03.q1k bases 1 to 1576
        !)))))****(((***%%((((*(((+,**(((+**+,-...

## Barcode File<a id="sec-3-2" name="sec-3-2"></a>

MLSTEZ uses CSV file for importing barcode information. Users can use Microsoft Excel to generate CSV file using "Save As&#x2026;" and selecting "Comma Seperated Values .csv" as output format. MLSTEZ 2.0 supports asymmetric barcode design, so you can save more money and put more samples into one batch. The new barcode file contains two (symmetric design) or three columns (asymmetric design). MLSTEZ 2.0 will switch between different modes based on the columns provided in the barcode file. The second column contains the names of strains/samples related to the barcodes in the first column, and the second column contains barcode sequences, which do not contains padding sequence or universal primer.

### Example of barcode file for MLSTEZ 2.0<a id="sec-3-2-1" name="sec-3-2-1"></a>

    Isolates_A,gcgctctgtgtgcagc,gcgctctgtgtgcagc
    Isolates_B,agagtactacatatga,agagtactacatatga
    Isolates_C,cgtgtgcatagatcgc,cgtgtgcatagatcgc
    Isolates_D,atgtatctcgactgca,atgtatctcgactgca

### Example of barcode file MLSTEZ 1.0<a id="sec-3-2-2" name="sec-3-2-2"></a>

    gcgctctgtgtgcagc,StrainA
    agagtactacatatga,StrainB
    tcatgagtcgacacta,StrainC

## Primer File<a id="sec-3-3" name="sec-3-3"></a>

CSV file is used for importing primer information in MLSTEZ. The primer file contains three columns, which are locus name, upper primer of locus and lower primer of locus. Universal primer sequence should be removed from upper/lower primer sequences. 

### Example of primer file<a id="sec-3-3-1" name="sec-3-3-1"></a>

    LOCUS1,TCTAATCGAAATGGTCAAGG,CGCAGCTGTTCGTCTGGATA
    LOCUS2,AATCGTCAAGGAGACCAACG,CGTCACCAGACTTGACGAAC
    LOCUS3,GATGGTTATGAACGAGAGGT,CTTACAGTCAGTATCGGACT

## Output Folder<a id="sec-3-4" name="sec-3-4"></a>

The output folder is used to store the project file, project.nma, and other results including consensus sequence for each loci, unmapped reads and allele sequences for heterozygous loci. Same output folder cannot be used for different project, otherwise, the project file and other output files will be overwritten by the newly built project. 

## Advanced parameters<a id="sec-3-5" name="sec-3-5"></a>

All the advanced parameters are automatically saved in system after first time use. User can click "Advanced" button to expand the parameter panel in order to change the settings. The parameters are set according to following instructions:
1.  File Format (default "FASTQ"): File format of input sequence file.
2.  Score Type (default "Phred 32"): Phred quality score of FASTQ file. If FASTA file is used as input sequence file, this option will not be valid.
3.  MUSCLE: Full path name (including file name) of MUSCLE. For example: /Users/YOURUSERNAME/bin/muscle-3.6/muscle on Mac OS or c:\muscle-3.6\muscle<sub>i86win32</sub>.exe on Windows. MUSCLE can be downloaded from <http://www.drive5.com/muscle/downloads.htm> for free.
4.  Padding Sequence (default "GGTAG"): The padding sequence of the barcode primer for the second PCR round. Please check the reference for more details.
5.  Universal Primer (default "CTGGAGCACGAGGACACTGA"): The universal primer sequence of the first and second PCR rounds. Please check the reference for more details.
6.  Barcode Length (default 16): Length barcode sequence in the barcode primer for the second PCR rounds. This length should be barcode sequence only, which does not count padding sequence and universal primer sequence in.
7.  Min Read Depth (default 3): The minimal number of reads used to generate the consensus sequence for each locus. With smaller number, the user can obtain more consensus sequence on low coverage locus, but some of them may have higher sequencing error rate.
8.  Max Read Depth (default 10): The maximal number of reads used to generate the consensus sequence for each locus. Larger number usage can lower the software effiency and may or may not increase the accuracy of the consensus sequence. Please check the reference for more details.
9.  Flanking Length (default 5): The flanking region length that the software used to search for barcode and primer sequences in sequencing reads. For example, when flanking length 5 is used, if the padding sequence length is 5 bp, and the barcode sequence length is 16 bp, the software will use the region between 0 (5-5) - 21 (16+5) bp on the 5' of sequencing read and the corresponding region on 3' of the read to look for the barcode sequence.
10. Match Score (default 2): The match score of the Smith-Waterman algorithm for barcode and primer identification.
11. Mismatch Score (default -1): The mismatch score of the Smith-Waterman algorithm for barcode and primer identification.
12. Gap Score (default -1): The gap score of the Smith-Waterman algorithm for barcode and primer identification.
13. Max Mismatch (default 3): The maximal number of the mismatches occurs in one alignment.
14. Threads (default 1): The thread number that can be used for parallel search barcode/primer sequence and heterozygous locus. Each thread runs on different processor/core and all threads run in parallel. Please select proper number based on processors/cores based on the computer's hardware.

# Run project<a id="sec-4" name="sec-4"></a>

MLSTEZ has four major functions in data analysis:
1.  Barcode and primer identification
2.  Generate consensus sequences
3.  Dump unmapped reads
4.  Heterozygous loci identification

After user has created a project, user can click "project setting" button in toolbar or select "Job settings" in "Project" menu to choose programs in the analysis. User can select "Run the whole process", which will run all four programs one by one. Alternatly, user also can selected one to several programs they interested in. After programs are selected, the "Run" button/menu will be enabled. 

## Barcode and primer identification<a id="sec-4-1" name="sec-4-1"></a>

This function is used to identify barcodes and primers in the sequencing reads, and it is necessary for the other three functions. Smith-Waterman algorithm is used for identification the barcode and primer sequences in the reads based on user settings (see "Advanced parameters"). After the analysis, MLSTEZ will show "Read Length", "Alignment Ratio", "Length Range" and "Sample Stats" in the software interface.
-   Read length dist: Barplot shows length distribution of all sequencing reads in the project
-   Alignment ratio: Pie chart shows the ratio of barcode and primer identified reads in all reads
-   Length ranges: Boxplot shows the length distribution of each locus. Blue "+" stands for outlier.
-   Read stats: Table shows the read number of each locus of each sample and total read number of each sample. The number in the brackets is the total read number identified for certain locus of the sample, and the number outside the brackets is the valid read number, which stands for the number of reads after filtered by the length range. If the valid read number is less than the minimal read number for generate the consensus (see "Advanced parameters"), this locus will shows with grey background, which means no consensus sequence will be generated for this locus. Double click on the locus grid will open a window, which shows the "Trimmed" reads (reads without primer sequences) and "Untrimmed" reads (raw reads).

## Generate consensus sequences<a id="sec-4-2" name="sec-4-2"></a>

The consensus sequences will be generated for the locus with more than minimal number of valid reads. All the consensus sequences will be output into the "Output Folder" in FASTA format automatically. Each locus will be output as single file named with "cons.LOCUSNAME.fasta", and each consensus is named with the locus name and sample name. The barcode and primer sequences have removed from the consensus, and the sequence directions have been adjusted based on the given upper and lower primers. A table of all sample loci will be shown in the main interface after analysis is completed. Double click on the the sample name shows all the generated consensus sequences for the sample, and double click on single locus shows the corresponding consensus sequence. 

## Dump unmapped reads<a id="sec-4-3" name="sec-4-3"></a>

The reads failed to identify by barcode or primer sequences can be output by this function for further analysis. The output file is stored in "Output Folder" named with "UnmappedReads.seq". Sequences after <NoBarcode> are the reads that are failed to identify barcode on one or both ends. Sequences after certain barcode sequences are the reads that are failed to identify primer sequences on one or both ends after barcode indentification.

## Search for heterozygous locus<a id="sec-4-4" name="sec-4-4"></a>

MLSTEZ can identify possible heterozygous locus based on the sequence differences among the reads. Five valid reads is the minimal requirements for the analysis. If two different sequence clusters are identified, the software will generate consensus sequences for both clusters. A possible heterozygous locus requires more than three nucleotide differences among the two concensus sequences. A table view of the result are shown in the main interface. The locus labelled and "Yes" and orange background indicate this locus might be a heterozygous locus. Double click the grid will show the consensus sequences for both alleles. The locus labelled with "NA" indicates the read number is less than the minimal requirements for this analysis. All the consensus sequences are saved in "Het.cons.fasta" under "Output Folder". The consensus sequences are named with SAMPLENAME<sub>LOCUSNAME</sub><sub>allele1</sub>/2.

# Open Project<a id="sec-5" name="sec-5"></a>

All of the project information is saved in "project.nma" under "Output Folder" automatically. User can open a existing project by clicking "Open Project" button in toolbar or by selecting "Open Project" in "Project" menu and then select the folder that was used as "Output Folder" before. All the results will be loaded in the software, and the user can even run the analysis functions that have not been processed. 

# Merge Projects<a id="sec-6" name="sec-6"></a>

The project merge function is designed for the samples that have been sequenced more than once in different batches in order to get higher read coverage. Different project can be merged together based on the sample names. The sequence reads identified by the sample name in different project will be merged together. After the merge step, user can "generate consensus sequences", "dump the unmapped reads" and identify the heteozygous locus using the merged data.
