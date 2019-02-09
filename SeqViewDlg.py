from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import QDialog
import ui_seqviewdlg
import sys

class SeqViewDlg(QDialog, ui_seqviewdlg.Ui_SeqViewDlg):
    def __init__(self,mode=0,trimmedSeqs="",untrimmedSeqs="",parent=None):
        super(SeqViewDlg,self).__init__(parent)
        self.setupUi(self)
        self.trimmedSeqs = trimmedSeqs
        self.untrimmedSeqs = untrimmedSeqs
        self.trimmedRadioButton.setChecked(True)
        self.textEdit.append(self.trimmedSeqs)
        self.textEdit.setReadOnly(True)
        self.currstat = 0
        if(mode == 1):
            self.trimmedRadioButton.setHidden(True)
            self.untrimmedRadioButton.setHidden(True)
        #self.setSize(610, 425)
    
    def on_trimmedRadioButton_clicked(self):
        if(self.currstat == 0):
            return
        else:
            self.textEdit.clear()
            self.textEdit.append(self.trimmedSeqs)
            self.currstat = 0
    
    def on_untrimmedRadioButton_clicked(self):
        if(self.currstat == 1):
            return
        else:
            self.textEdit.clear()
            self.textEdit.append(self.untrimmedSeqs)
            self.currstat = 1

if(__name__ == "__main__"):
    app = QApplication(sys.argv)
    mainwindow = SeqViewDlg()
    mainwindow.show()
    app.exec_()
            
