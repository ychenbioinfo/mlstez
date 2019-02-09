from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QFileDialog, QDialog, QListWidgetItem, QMessageBox, QApplication
import ui_parameterdlg
import Parameters
import sys
from os.path import isfile
from os.path import isdir
from os.path import dirname


class ParameterDlg(QDialog, ui_parameterdlg.Ui_ParameterDlg):
    def __init__(self, parameters, parent=None):
        super(ParameterDlg, self).__init__(parent)
        self.parameters = parameters
        self.setupUi(self)
        self.workingpath = dirname(__file__)
        self.initUi()

    def initUi(self):
        """load parameters into ui"""
        # define window size based on parameters we have
        if not isinstance(self.parameters.MuscleCMD, str):
            self.parameters.MuscleCMD = ""
        if not isinstance(self.parameters.PadSeq, str):
            self.parameters.PadSeq = ""
        if not isinstance(self.parameters.UniPrimer, str):
            self.parameters.UniPrimer = ""
            
        if self.parameters.isEssential():
            self.moreButton.setChecked(False)
            self.setFixedSize(560,240)
        else:
            self.moreButton.setChecked(True)
            self.setFixedSize(560,455)

        self.removeseqfileButton.setEnabled(False)
        
        # init all values based on parameter information
        self.projnameLineEdit.setText(self.parameters.ProjectName)
        if len(self.parameters.Seq_Files) > 0:
            for filen in self.parameters.Seq_Files:
                self.__addSeqFiles(filen)
        self.barcodefileLineEdit.setText(self.parameters.Barcode_File)
        self.primerfileLineEdit.setText(self.parameters.Primer_File)
        self.outfolderLineEdit.setText(self.parameters.Out_Folder)
        if self.parameters.Filetype == "FASTQ":
            self.fastqRadioBtn.setChecked(True)
        else:
            self.fastaRadioBtn.setChecked(True)
        if self.parameters.ScoringSys == "phred33":
            self.phred33RadioBtn.setChecked(True)
        else:
            self.phred64RadioBtn.setChecked(True)
        
        self.muscleLineEdit.setText(self.parameters.MuscleCMD)
        self.paddingEditLine.setText(self.parameters.PadSeq)
        self.uniprimerLineEdit.setText(self.parameters.UniPrimer)
        self.barcodelenSpinBox.setValue(self.parameters.BarcodeLen)
        self.mindepthSpinBox.setValue(self.parameters.MinReadNum)
        self.maxdepthSpinBox.setValue(self.parameters.MaxReadNum)
        self.flankinglengthSpinBox.setValue(self.parameters.FlankingLength)
        self.matchscoreSpinBox.setValue(self.parameters.MatchScore)
        self.mismatchscireSpinBox.setValue(self.parameters.MismatchScore)
        self.gapscoreSpinBox.setValue(self.parameters.GapScore)
        self.maxmismatchSpinBox.setValue(self.parameters.MaxMisMatch)
        self.threadSpinBox.setValue(self.parameters.Threads)
        #self.consensuscutDoubleSpinBox.setValue(self.parameters.ConsensusCut)
        if isdir(self.parameters.Out_Folder):
            self.workingpath = self.parameters.Out_Folder
    
    @pyqtSlot()
    def on_moreButton_clicked(self):
        if self.moreButton.isChecked():
            self.setFixedSize(560,455)
        else:
            self.setFixedSize(560,240)
    
    @pyqtSlot()
    def on_addseqfileButton_clicked(self):
        filename, _ = QFileDialog.getOpenFileName(parent=self, caption='Add file', directory=self.workingpath)
        if filename:
            self.__addSeqFiles(filename)
            self.workingpath = dirname(filename)
    
    @pyqtSlot()
    def on_barcodefileButton_clicked(self):
        filename, _ = QFileDialog.getOpenFileName(self, 'Barcode file',
                self.workingpath)
        if filename:
            self.barcodefileLineEdit.setText(filename)
            self.workingpath = dirname(filename)
    
    @pyqtSlot()
    def on_primerfileButton_clicked(self):
        filename, _ = QFileDialog.getOpenFileName(self, 'Primer file',
                self.workingpath)
        if filename:
            self.primerfileLineEdit.setText(filename)
            self.workingpath = dirname(filename)
    
    @pyqtSlot()
    def on_outfolderButton_clicked(self):
        foldername = QFileDialog.getExistingDirectory(self, 'Output Folder',
                self.workingpath)
        if foldername:
            self.outfolderLineEdit.setText(foldername)
            
    @pyqtSlot()
    def on_muscleButton_clicked(self):
        filename, _ = QFileDialog.getOpenFileName(self, 'MUSCLE',
                self.workingpath)
        if filename:
            self.muscleLineEdit.setText(filename)
            self.workingpath = dirname(filename)
    
    @pyqtSlot()
    def on_removeseqfileButton_clicked(self):
        self.seqfileListWidget.takeItem(self.seqfileListWidget.currentRow())
        if self.seqfileListWidget.currentItem() is None:
            self.removeseqfileButton.setEnabled(False)
    
    @pyqtSlot()
    def on_okayButton_clicked(self):
        if not self.__checkEditLine(self.projnameLineEdit, "project name"):
            return
        if self.seqfileListWidget.count() < 1:
            self.__sendMessage("Sequencing file")
            self.addseqfileButton.setFocus()
            return
        if not self.__checkEditLine(self.barcodefileLineEdit, "barcode file"):
            return
        if not self.__checkEditLine(self.primerfileLineEdit, "primer file"):
            return
        if not self.__checkEditLine(self.outfolderLineEdit, "output folder"):
            return
        if not self.__checkEditLine(self.muscleLineEdit, "MUSCLE command"):
            return
        if not self.__checkEditLine(self.paddingEditLine, "padding sequence"):
            return
        if not self.__checkEditLine(self.uniprimerLineEdit, "unviersal primer"):
            return
        
        self.parameters.Seq_Files, isOkay = self.__getSeqFiles()
        if not isOkay:
            return
        
        if isdir((self.outfolderLineEdit.text())):
            self.parameters.Out_Folder = self.outfolderLineEdit.text()
        else:
            return
        
        if self.__checkFileExist(self.primerfileLineEdit.text()):
            self.parameters.Primer_File = self.primerfileLineEdit.text()
        else:
            return
        
        if self.__checkFileExist(self.barcodefileLineEdit.text()):
            self.parameters.Barcode_File = self.barcodefileLineEdit.text()
        else:
            return
        
        if self.__checkFileExist(self.muscleLineEdit.text()):
            self.parameters.MuscleCMD = self.muscleLineEdit.text()
        else:
            return
        
        self.parameters.ProjectName = self.projnameLineEdit.text()
        if self.fastqRadioBtn.isChecked():
            self.parameters.Filetype = "FASTQ"
        else:
            self.parameters.Filetype = "FASTA"
        if self.phred33RadioBtn.isChecked():
            self.parameters.ScoringSys = "phred33"
        else:
            self.parameters.ScoringSys = "phred64"

        self.parameters.PadSeq = self.paddingEditLine.text().upper()
        self.parameters.PadLength = len(self.parameters.PadSeq)
        self.parameters.UniPrimer = self.uniprimerLineEdit.text().upper()
        self.parameters.UniLength = len(self.parameters.UniPrimer)
        self.parameters.BarcodeLen = self.barcodelenSpinBox.value()
        self.parameters.MinReadNum = self.mindepthSpinBox.value()
        self.parameters.MaxReadNum = self.maxdepthSpinBox.value()
        self.parameters.FlankingLength = self.flankinglengthSpinBox.value()
        self.parameters.MatchScore = self.matchscoreSpinBox.value()
        self.parameters.MismatchScore = self.mismatchscireSpinBox.value()
        self.parameters.GapScore = self.gapscoreSpinBox.value()
        self.parameters.MaxMisMatch = self.maxmismatchSpinBox.value()
        self.parameters.Threads = self.threadSpinBox.value()
        #self.parameters.ConsensusCut = self.consensuscutDoubleSpinBox.value()
        self.parameters.update()
        
        self.accept()
    
    def on_seqfileListWidget_itemClicked(self):
        self.removeseqfileButton.setEnabled(True)

    def __addSeqFiles(self,filename):
        fileItem = QListWidgetItem()
        fileItem.setText(filename)
        self.seqfileListWidget.addItem(fileItem)
    
    def __sendMessage(self, msg):
        QMessageBox.warning(self, "Parameter missing",
                                    ("%s is empty" %(msg)), QMessageBox.Ok|QMessageBox.Default)
        
    def __checkEditLine(self,edtline,msg):
        if edtline.text() == "":
            self.__sendMessage(msg)
            edtline.setFocus()
            return False
        else:
            return True
    
    def __getSeqFiles(self):
        filecount = self.seqfileListWidget.count()
        filenames = []
        for count in range(filecount):
            item = self.seqfileListWidget.item(count)
            filename = item.text()
            if not self.__checkFileExist(filename):
                return (None, False)
            filenames.append(filename)
        return filenames, True
    
    def __checkFileExist(self, fn):
        if isfile(fn):
            return True
        else:
            QMessageBox.warning(self, "Error",
                                ("Cannot find file: %s" %fn))
            return False


if __name__ == "__main__":
    app = QApplication(sys.argv)
    paras = Parameters.Parameters()
    #paras.Seq_Files = ["/home/aa.fastq","/home/bb.fastq"]
    paras.Filetype = "FASTA"
    paras.MuscleCMD = "/Users/Alvin/bin/MUSCLE/muscle"
    paras.PadSeq = "AATTG"
    paras.UniPrimer = "ATGAGTAAGTATGATTGAAGT"
    mainwindow = ParameterDlg(paras)
    mainwindow.show()
    app.exec_()
