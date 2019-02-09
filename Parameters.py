#!/usr/bin/env

import configparser
class Parameters(object):
  
    def __init__(self):
        self.ProjectName = ""
        self.Seq_Files = []
        self.Out_Folder = ""
        self.Primer_File = ""
        self.Barcode_File = ""
        self.Filetype = "FASTQ"
        self.MuscleCMD = ""
        self.ScoringSys = "phred33"
        self.PadSeq = "GGTAG"
        self.PadLength = 5
        self.UniPrimer = "CTGGAGCACGAGGACACTGA"
        self.UniLength = 20
        self.BarcodeLen = 16
        self.MinReadNum = 3
        self.MaxReadNum = 10
        self.FlankingLength = 5
        self.MatchScore = 2
        self.MismatchScore = -1
        self.GapScore = -1
        self.MaxMisMatch = 3
        self.Threads = 1
        #self.ConsensusCut = 0.5
        self.EndLength = 0

    
    def update(self):
        self.PadLength = len(self.PadSeq)
        self.UniLength = len(self.UniPrimer)
        self.EndLength = self.PadLength + self.BarcodeLen + self.FlankingLength
    
    def isEssential(self):
        if(self.MuscleCMD == ""):
            return False
        if(self.PadSeq == ""):
            return False
        if(self.UniPrimer == ""):
            return False
        return True

    def savefile(self,fn):
        config = configparser.RawConfigParser()
        config.add_section('PROJECT')
        config.set('PROJECT','Project_Name', self.ProjectName)
        config.add_section('FILES')
        config.set('FILES', 'Sequencing_Files', ','.join(self.Seq_Files))
        config.set('FILES', 'Output_Folder', self.Out_Folder)
        config.set('FILES', 'Primer_File', self.Primer_File)
        config.set('FILES', 'Barcode_File', self.Barcode_File)
        config.add_section('SETTINGS')
        config.set('SETTINGS', 'File_Type', self.Filetype)
        config.set('SETTINGS', 'Muscle_Command', self.MuscleCMD)
        config.set('SETTINGS', 'Score_Type', self.ScoringSys)
        config.set('SETTINGS', 'Padding_Seq', self.PadSeq)
        config.set('SETTINGS', 'Universal_Primer', self.UniPrimer)
        config.set('SETTINGS', 'Barcode_Length', self.BarcodeLen)
        config.set('SETTINGS', 'Min_ReadNum', self.MinReadNum)
        config.set('SETTINGS', 'Max_ReadNum', self.MaxReadNum)
        config.set('SETTINGS', 'Flanking_Length', self.FlankingLength)
        config.set('SETTINGS', 'Match_Score', self.MatchScore)
        config.set('SETTINGS', 'Mismatch_Score', self.MismatchScore)
        config.set('SETTINGS', 'Gap_Score', self.GapScore)
        config.set('SETTINGS', 'Max_Mismatch', self.MaxMisMatch)
        config.set('SETTINGS', 'Threads', self.Threads)
        #config.set('SETTINGS', 'Consensus_Cut', self.ConsensusCut)
        
        try:
            with open(fn, 'w') as configfile:
                config.write(configfile)
                return (True, None)
        except IOError as e:
            return (False, e)
            
    def openfile(self,fn):
        config = configparser.ConfigParser()
        try:
            config.read(fn)
            #---python3 
            #ProjectConfig = config['PROJECT']
            #FileConfig = config['FILES']
            #SettingConfig = config['SETTINGS']
            
            #self.ProjectName = ProjectConfig.get('Project_Name','MLST_Project')
            #Files = FileConfig['Sequencing_Files']
            #self.Seq_Files = Files.split(',')
            #self.Out_Folder = config.get('FILES','Output_Folder','output')
            #self.Primer_File = config.get('FILES','Primer_File','primers')
            #self.Barcode_File = config.get('FILES','Barcode_File','barcodes')
            #self.Filetype = config.get('SETTINGS','File_Type','FASTQ')
            #self.Filetype = self.Filetype.upper()
            #self.MuscleCMD = config.get('SETTINGS','Muscle_Command','muscle')
            #self.ScoringSys = config.get('SETTINGS','Score_Type','phred33')
            #self.PadSeq = config.get('SETTINGS','Padding_Seq','GGTAG')
            #self.PadSeq = self.PadSeq.upper()
            #self.PadLength = len(self.PadSeq)
            #self.UniPrimer = config.get('SETTINGS','Universal_Primer','CTGGAGCACGAGGACACTGA')
            #self.UniPrimer = self.UniPrimer.upper()
            #self.UniLength = len(self.UniPrimer)
            #self.BarcodeLen = int(config.get('SETTINGS','Barcode_Length','16'))
            #self.MinReadNum = int(config.get('SETTINGS','Min_ReadNum','5'))
            #self.MaxReadNum = int(config.get('SETTINGS','Max_ReadNum','10'))
            #self.FlankingLength = int(config.get('SETTINGS','Flanking_Length','5'))
            #self.MatchScore = int(config.get('SETTINGS','Match_Score','2'))
            #self.MismatchScore = int(config.get('SETTINGS','Mismatch_Score','-1'))
            #self.GapScore = int(config.get('SETTINGS','Gap_Score','-1'))
            #self.MaxMisMatch = int(config.get('SETTINGS','Max_Mismatch','3'))
            #self.ConsensusCut = float(config.get('SETTINGS','Consensus_Cut','0.5'))
            #self.EndLength = self.PadLength + self.BarcodeLen + self.FlankingLength
            
            #--python 2.7
            self.ProjectName = config.get('PROJECT','Project_Name','MLST_Project')
            Files = config.get('FILES','Sequencing_Files')
            self.Seq_Files = Files.split(',')
            self.Out_Folder = config.get('FILES','Output_Folder','output')
            self.Primer_File = config.get('FILES','Primer_File','primers')
            self.Barcode_File = config.get('FILES','Barcode_File','barcodes')
            self.Filetype = config.get('SETTINGS','File_Type','FASTQ')
            self.Filetype = self.Filetype.upper()
            self.MuscleCMD = config.get('SETTINGS','Muscle_Command','muscle')
            self.ScoringSys = config.get('SETTINGS','Score_Type','phred33')
            self.PadSeq = config.get('SETTINGS','Padding_Seq','GGTAG')
            self.PadSeq = self.PadSeq.upper()
            self.PadLength = len(self.PadSeq)
            self.UniPrimer = config.get('SETTINGS','Universal_Primer','CTGGAGCACGAGGACACTGA')
            self.UniPrimer = self.UniPrimer.upper()
            self.UniLength = len(self.UniPrimer)
            self.BarcodeLen = int(config.get('SETTINGS','Barcode_Length','16'))
            self.MinReadNum = int(config.get('SETTINGS','Min_ReadNum','5'))
            self.MaxReadNum = int(config.get('SETTINGS','Max_ReadNum','10'))
            self.FlankingLength = int(config.get('SETTINGS','Flanking_Length','5'))
            self.MatchScore = int(config.get('SETTINGS','Match_Score','2'))
            self.MismatchScore = int(config.get('SETTINGS','Mismatch_Score','-1'))
            self.GapScore = int(config.get('SETTINGS','Gap_Score','-1'))
            self.MaxMisMatch = int(config.get('SETTINGS','Max_Mismatch','3'))
            self.Threads = int(config.get('SETTINGS','Threads','1'))
            #self.ConsensusCut = float(config.get('SETTINGS','Consensus_Cut','0.5'))
            self.EndLength = self.PadLength + self.BarcodeLen + self.FlankingLength
            
            return True
        except IOError as e:
            return (False,e)