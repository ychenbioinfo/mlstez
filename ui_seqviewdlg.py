# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'seqviewdlg.ui'
#
# Created: Sat Apr 19 16:45:19 2014
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

class Ui_SeqViewDlg(object):
    def setupUi(self, SeqViewDlg):
        SeqViewDlg.setObjectName(_fromUtf8("SeqViewDlg"))
        SeqViewDlg.resize(625, 426)
        SeqViewDlg.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.verticalLayout_2 = QtGui.QVBoxLayout(SeqViewDlg)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.textEdit = QtGui.QTextEdit(SeqViewDlg)
        self.textEdit.setObjectName(_fromUtf8("textEdit"))
        self.verticalLayout.addWidget(self.textEdit)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.trimmedRadioButton = QtGui.QRadioButton(SeqViewDlg)
        self.trimmedRadioButton.setObjectName(_fromUtf8("trimmedRadioButton"))
        self.buttonGroup = QtGui.QButtonGroup(SeqViewDlg)
        self.buttonGroup.setObjectName(_fromUtf8("buttonGroup"))
        self.buttonGroup.addButton(self.trimmedRadioButton)
        self.horizontalLayout.addWidget(self.trimmedRadioButton)
        self.untrimmedRadioButton = QtGui.QRadioButton(SeqViewDlg)
        self.untrimmedRadioButton.setObjectName(_fromUtf8("untrimmedRadioButton"))
        self.buttonGroup.addButton(self.untrimmedRadioButton)
        self.horizontalLayout.addWidget(self.untrimmedRadioButton)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.closeButton = QtGui.QPushButton(SeqViewDlg)
        self.closeButton.setObjectName(_fromUtf8("closeButton"))
        self.horizontalLayout.addWidget(self.closeButton)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.verticalLayout_2.addLayout(self.verticalLayout)

        self.retranslateUi(SeqViewDlg)
        QtCore.QObject.connect(self.closeButton, QtCore.SIGNAL(_fromUtf8("clicked()")), SeqViewDlg.reject)
        QtCore.QMetaObject.connectSlotsByName(SeqViewDlg)
        SeqViewDlg.setTabOrder(self.trimmedRadioButton, self.untrimmedRadioButton)
        SeqViewDlg.setTabOrder(self.untrimmedRadioButton, self.closeButton)
        SeqViewDlg.setTabOrder(self.closeButton, self.textEdit)

    def retranslateUi(self, SeqViewDlg):
        SeqViewDlg.setWindowTitle(_translate("SeqViewDlg", "Dialog", None))
        self.trimmedRadioButton.setText(_translate("SeqViewDlg", "Trimmed", None))
        self.untrimmedRadioButton.setText(_translate("SeqViewDlg", "Untrimmed", None))
        self.closeButton.setText(_translate("SeqViewDlg", "Close", None))

