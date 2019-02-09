# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'seqviewdlg.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_SeqViewDlg(object):
    def setupUi(self, SeqViewDlg):
        SeqViewDlg.setObjectName("SeqViewDlg")
        SeqViewDlg.resize(625, 426)
        SeqViewDlg.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(SeqViewDlg)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.textEdit = QtWidgets.QTextEdit(SeqViewDlg)
        self.textEdit.setObjectName("textEdit")
        self.verticalLayout.addWidget(self.textEdit)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.trimmedRadioButton = QtWidgets.QRadioButton(SeqViewDlg)
        self.trimmedRadioButton.setObjectName("trimmedRadioButton")
        self.buttonGroup = QtWidgets.QButtonGroup(SeqViewDlg)
        self.buttonGroup.setObjectName("buttonGroup")
        self.buttonGroup.addButton(self.trimmedRadioButton)
        self.horizontalLayout.addWidget(self.trimmedRadioButton)
        self.untrimmedRadioButton = QtWidgets.QRadioButton(SeqViewDlg)
        self.untrimmedRadioButton.setObjectName("untrimmedRadioButton")
        self.buttonGroup.addButton(self.untrimmedRadioButton)
        self.horizontalLayout.addWidget(self.untrimmedRadioButton)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.closeButton = QtWidgets.QPushButton(SeqViewDlg)
        self.closeButton.setObjectName("closeButton")
        self.horizontalLayout.addWidget(self.closeButton)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.verticalLayout_2.addLayout(self.verticalLayout)

        self.retranslateUi(SeqViewDlg)
        self.closeButton.clicked.connect(SeqViewDlg.reject)
        QtCore.QMetaObject.connectSlotsByName(SeqViewDlg)
        SeqViewDlg.setTabOrder(self.trimmedRadioButton, self.untrimmedRadioButton)
        SeqViewDlg.setTabOrder(self.untrimmedRadioButton, self.closeButton)
        SeqViewDlg.setTabOrder(self.closeButton, self.textEdit)

    def retranslateUi(self, SeqViewDlg):
        _translate = QtCore.QCoreApplication.translate
        SeqViewDlg.setWindowTitle(_translate("SeqViewDlg", "Dialog"))
        self.trimmedRadioButton.setText(_translate("SeqViewDlg", "Trimmed"))
        self.untrimmedRadioButton.setText(_translate("SeqViewDlg", "Untrimmed"))
        self.closeButton.setText(_translate("SeqViewDlg", "Close"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    SeqViewDlg = QtWidgets.QDialog()
    ui = Ui_SeqViewDlg()
    ui.setupUi(SeqViewDlg)
    SeqViewDlg.show()
    sys.exit(app.exec_())

