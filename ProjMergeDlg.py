from PyQt4.QtCore import *
from PyQt4.QtGui import *
import ui_projmergedlg
from os.path import isfile
from os.path import isdir
import sys

class ProjMergeDlg(QDialog, ui_projmergedlg.Ui_ProjMergeDlg):
    def __init__(self, curfolder = "", parent=None):
        super(ProjMergeDlg, self).__init__(parent)
        self.setupUi(self)
        self.projName = ""
        self.outFolder = ""
        self.projFiles = []
        self.projFileName = "project.nma"
        self.workingpath = curfolder
        self.removeButton.setEnabled(False)
        self.setFixedSize(425, 285)
        
    @pyqtSlot()
    def on_selectButton_clicked(self):
        foldername = str(QFileDialog.getExistingDirectory(self, 'Output Folder', 
                self.workingpath))
        if(foldername):
            self.outfolderEdit.setText(foldername)
            self.workingpath = foldername

    @pyqtSlot()
    def on_addButton_clicked(self):
        foldername = str(QFileDialog.getExistingDirectory(self, 'Add Project Folder', 
                self.workingpath))
        if(foldername):
            self.__addProjPath(foldername)
            self.workingpath = foldername
    
    @pyqtSlot()
    def on_removeButton_clicked(self):
        self.projsListWidget.takeItem(self.projsListWidget.currentRow())
        if(self.projsListWidget.currentItem() is None):
            self.removeButton.setEnabled(False)
            
    @pyqtSlot()
    def on_okayButton_clicked(self):
        if(self.__checkEditLine(self.projnameEdit,"Project Name")):
            self.projName = self.projnameEdit.text()
        else:
            return
        
        if(self.__checkEditLine(self.outfolderEdit,"Project Folder")):
            if(not isdir(self.outfolderEdit.text())):
                QMessageBox.warning(self, "Error",
                                    ("Project folder does not exist!"),
                                    QMessageBox.Ok|QMessageBox.Default)
                return
            self.outFolder = self.outfolderEdit.text()
        else:
            return

        self.projFiles = self.__getSeqFiles()
        self.projFiles = list(set(self.projFiles))
        if(len(self.projFiles) < 2):
            QMessageBox.warning(self, "Error",
                                    ("At least two different projects required"),
                                    QMessageBox.Ok|QMessageBox.Default)
            return
        self.accept()


    def on_projsListWidget_itemClicked(self):
        self.removeButton.setEnabled(True)
    
    def __addProjPath(self,foldername):
        projfile = foldername + "/" + self.projFileName
        if(not isfile(projfile)):
            QMessageBox.warning(self, "Error:",
                                    ("Cannot project file in folder %s" %(foldername)), QMessageBox.Ok|QMessageBox.Default)
            return
        fileItem = QListWidgetItem()
        fileItem.setText(foldername)
        self.projsListWidget.addItem(fileItem)
    
    def __getSeqFiles(self):
        foldercount = self.projsListWidget.count()
        projfiles = []
        for count in range(foldercount):
            item = self.projsListWidget.item(count)
            foldername = str(item.text())
            projfile = foldername + "/" + self.projFileName
            projfiles.append(projfile)
        return (projfiles)
    
        
    def __checkEditLine(self,edtline,msg):
        if(edtline.text() == ""):
            self.__sendMessage(msg)
            edtline.setFocus()
            return False
        else:
            return True
        
    def __sendMessage(self, msg):
        QMessageBox.warning(self, "Parameter missing",
                                    ("%s is empty" %(msg)), QMessageBox.Ok|QMessageBox.Default)

if(__name__ == "__main__"):
    app = QApplication(sys.argv)
    mainwindow = ProjMergeDlg()
    mainwindow.show()
    app.exec_()
