# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'aboutdlg.ui'
#
# Created: Wed Jul 30 16:29:26 2014
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_AboutDlg(object):
    def setupUi(self, AboutDlg):
        AboutDlg.setObjectName(_fromUtf8("AboutDlg"))
        AboutDlg.resize(452, 278)
        self.label_4 = QtGui.QLabel(AboutDlg)
        self.label_4.setGeometry(QtCore.QRect(30, 150, 411, 71))
        self.label_4.setWordWrap(True)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.pushButton = QtGui.QPushButton(AboutDlg)
        self.pushButton.setGeometry(QtCore.QRect(330, 240, 91, 32))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.label_3 = QtGui.QLabel(AboutDlg)
        self.label_3.setGeometry(QtCore.QRect(30, 120, 201, 16))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.label_6 = QtGui.QLabel(AboutDlg)
        self.label_6.setGeometry(QtCore.QRect(30, 20, 71, 71))
        self.label_6.setText(_fromUtf8(""))
        self.label_6.setPixmap(QtGui.QPixmap(_fromUtf8(":/logo.png")))
        self.label_6.setScaledContents(True)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.layoutWidget = QtGui.QWidget(AboutDlg)
        self.layoutWidget.setGeometry(QtCore.QRect(240, 16, 171, 81))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.gridLayout = QtGui.QGridLayout(self.layoutWidget)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_7 = QtGui.QLabel(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(17)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_7.setFont(font)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.gridLayout.addWidget(self.label_7, 0, 0, 1, 1, QtCore.Qt.AlignHCenter)
        self.label_5 = QtGui.QLabel(self.layoutWidget)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout.addWidget(self.label_5, 1, 0, 1, 1, QtCore.Qt.AlignHCenter)
        self.label_2 = QtGui.QLabel(self.layoutWidget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 2, 0, 1, 1, QtCore.Qt.AlignHCenter)

        self.retranslateUi(AboutDlg)
        QtCore.QObject.connect(self.pushButton, QtCore.SIGNAL(_fromUtf8("clicked()")), AboutDlg.reject)
        QtCore.QMetaObject.connectSlotsByName(AboutDlg)

    def retranslateUi(self, AboutDlg):
        AboutDlg.setWindowTitle(_translate("AboutDlg", "About MLSTEZ", None))
        self.label_4.setText(_translate("AboutDlg", "Reference: Next Generation Multilocus Sequence Typing (NGMLST) and the Analytical Software Program MLSTEZ Enable Efficient, Cost-Effective, High-Throughput, Multilocus Sequencing Typing. Y Chen, et al., (submitted)", None))
        self.pushButton.setText(_translate("AboutDlg", "OK", None))
        self.label_3.setText(_translate("AboutDlg", "Email: ychenbioinfo@gmail.com", None))
        self.label_7.setText(_translate("AboutDlg", "MLSTEZ", None))
        self.label_5.setText(_translate("AboutDlg", "(version 0.1)", None))
        self.label_2.setText(_translate("AboutDlg", "Created by: Yuan Chen", None))

import resource_rc
