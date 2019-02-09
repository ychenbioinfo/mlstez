from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import QDialog
import ui_progconfigdlg
import sys

class ProgConfigDlg(QDialog, ui_progconfigdlg.Ui_ProgConfigDlg):
    def __init__(self, configcode=0, parent=None):
        super(ProgConfigDlg, self).__init__(parent)
        self.setupUi(self)
        self.configcode = configcode
        self.uiarray = [self.alignCheckBox, self.consensusCheckBox,
                   self.unmapCheckBox, self.hybridCheckBox]
        if (configcode > 0):
            self.runallRadioButton.setEnabled(False)
            self.runselectedRadioButton.setChecked(True)
            self.__run_selected()
        else:
            self.runallRadioButton.setChecked(True)
            self.__run_all()
        self.setFixedSize(415, 290)
        
    def __run_all(self):
        for i in range(4):
            self.uiarray[i].setEnabled(False)
    
    def __run_selected(self):
        self.alignCheckBox.setEnabled(False)
        self.alignCheckBox.setChecked(True)
        for i in range(1,4):
            mask = 1 << i
            if(self.configcode & mask):
                self.uiarray[i].setEnabled(False)
                self.uiarray[i].setChecked(True)
            else:
                self.uiarray[i].setEnabled(True)
    
    def on_runselectedRadioButton_clicked(self):
        self.__run_selected()
    
    def on_runallRadioButton_clicked(self):
        self.__run_all()
        
    @pyqtSlot()
    def on_okayButton_clicked(self):
        if(self.runallRadioButton.isChecked()):
            self.retcode = int('1111',2) 
        else:
            retcode = 0
            for i in range(len(self.uiarray)):
                if(self.uiarray[i].isChecked()):
                    retcode = retcode | 1 << i
            if(retcode < 1):
                QMessageBox.warning(self, "Error",("You haven't select any program yet!"))
                return
            else:
                self.retcode = retcode
        self.accept()
    
def main():
    app = QApplication(sys.argv)
    mainwindow = ProgConfigDlg(int('0011',2))
    mainwindow.show()
    app.exec_()

if(__name__ == '__main__'):
    main()
    