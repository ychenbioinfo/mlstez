# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'aboutdlg.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_AboutDlg(object):
    def setupUi(self, AboutDlg):
        AboutDlg.setObjectName("AboutDlg")
        AboutDlg.resize(452, 278)
        self.label_4 = QtWidgets.QLabel(AboutDlg)
        self.label_4.setGeometry(QtCore.QRect(30, 150, 411, 71))
        self.label_4.setWordWrap(True)
        self.label_4.setObjectName("label_4")
        self.pushButton = QtWidgets.QPushButton(AboutDlg)
        self.pushButton.setGeometry(QtCore.QRect(330, 240, 91, 32))
        self.pushButton.setObjectName("pushButton")
        self.label_3 = QtWidgets.QLabel(AboutDlg)
        self.label_3.setGeometry(QtCore.QRect(30, 120, 201, 16))
        self.label_3.setObjectName("label_3")
        self.label_6 = QtWidgets.QLabel(AboutDlg)
        self.label_6.setGeometry(QtCore.QRect(30, 20, 71, 71))
        self.label_6.setText("")
        self.label_6.setPixmap(QtGui.QPixmap(":/logo.png"))
        self.label_6.setScaledContents(True)
        self.label_6.setObjectName("label_6")
        self.layoutWidget = QtWidgets.QWidget(AboutDlg)
        self.layoutWidget.setGeometry(QtCore.QRect(240, 16, 171, 81))
        self.layoutWidget.setObjectName("layoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.layoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.label_7 = QtWidgets.QLabel(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(17)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.gridLayout.addWidget(self.label_7, 0, 0, 1, 1, QtCore.Qt.AlignHCenter)
        self.label_5 = QtWidgets.QLabel(self.layoutWidget)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 1, 0, 1, 1, QtCore.Qt.AlignHCenter)
        self.label_2 = QtWidgets.QLabel(self.layoutWidget)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 0, 1, 1, QtCore.Qt.AlignHCenter)

        self.retranslateUi(AboutDlg)
        self.pushButton.clicked.connect(AboutDlg.reject)
        QtCore.QMetaObject.connectSlotsByName(AboutDlg)

    def retranslateUi(self, AboutDlg):
        _translate = QtCore.QCoreApplication.translate
        AboutDlg.setWindowTitle(_translate("AboutDlg", "About MLSTEZ"))
        self.label_4.setText(_translate("AboutDlg", "Reference: Next Generation Multilocus Sequence Typing (NGMLST) and the Analytical Software Program MLSTEZ Enable Efficient, Cost-Effective, High-Throughput, Multilocus Sequencing Typing. Y Chen, et al., (submitted)"))
        self.pushButton.setText(_translate("AboutDlg", "OK"))
        self.label_3.setText(_translate("AboutDlg", "Email: ychenbioinfo@gmail.com"))
        self.label_7.setText(_translate("AboutDlg", "MLSTEZ"))
        self.label_5.setText(_translate("AboutDlg", "(version 0.1)"))
        self.label_2.setText(_translate("AboutDlg", "Created by: Yuan Chen"))

import resource_rc

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    AboutDlg = QtWidgets.QDialog()
    ui = Ui_AboutDlg()
    ui.setupUi(AboutDlg)
    AboutDlg.show()
    sys.exit(app.exec_())

