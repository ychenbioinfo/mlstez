from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import QDialog
import ui_aboutdlg
import sys

class AboutDlg(QDialog,ui_aboutdlg.Ui_AboutDlg):
    def __init__(self,parent=None):
        super(AboutDlg,self).__init__(parent)
        self.setupUi(self)

if(__name__ == '__main__'):
    app = QApplication(sys.argv)
    mainwindow = AboutDlg()
    mainwindow.exec_()
    app.exec_()