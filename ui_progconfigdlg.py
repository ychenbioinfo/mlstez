# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'progconfigdlg.ui'
#
# Created: Sun Jun 29 19:28:53 2014
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

class Ui_ProgConfigDlg(object):
    def setupUi(self, ProgConfigDlg):
        ProgConfigDlg.setObjectName(_fromUtf8("ProgConfigDlg"))
        ProgConfigDlg.resize(415, 290)
        self.frame = QtGui.QFrame(ProgConfigDlg)
        self.frame.setGeometry(QtCore.QRect(40, 80, 341, 151))
        self.frame.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtGui.QFrame.Raised)
        self.frame.setObjectName(_fromUtf8("frame"))
        self.layoutWidget = QtGui.QWidget(self.frame)
        self.layoutWidget.setGeometry(QtCore.QRect(10, 10, 321, 131))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.layoutWidget)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.alignCheckBox = QtGui.QCheckBox(self.layoutWidget)
        self.alignCheckBox.setObjectName(_fromUtf8("alignCheckBox"))
        self.verticalLayout.addWidget(self.alignCheckBox)
        self.consensusCheckBox = QtGui.QCheckBox(self.layoutWidget)
        self.consensusCheckBox.setObjectName(_fromUtf8("consensusCheckBox"))
        self.verticalLayout.addWidget(self.consensusCheckBox)
        self.unmapCheckBox = QtGui.QCheckBox(self.layoutWidget)
        self.unmapCheckBox.setObjectName(_fromUtf8("unmapCheckBox"))
        self.verticalLayout.addWidget(self.unmapCheckBox)
        self.hybridCheckBox = QtGui.QCheckBox(self.layoutWidget)
        self.hybridCheckBox.setObjectName(_fromUtf8("hybridCheckBox"))
        self.verticalLayout.addWidget(self.hybridCheckBox)
        self.runallRadioButton = QtGui.QRadioButton(ProgConfigDlg)
        self.runallRadioButton.setGeometry(QtCore.QRect(40, 20, 171, 20))
        self.runallRadioButton.setObjectName(_fromUtf8("runallRadioButton"))
        self.buttonGroup = QtGui.QButtonGroup(ProgConfigDlg)
        self.buttonGroup.setObjectName(_fromUtf8("buttonGroup"))
        self.buttonGroup.addButton(self.runallRadioButton)
        self.runselectedRadioButton = QtGui.QRadioButton(ProgConfigDlg)
        self.runselectedRadioButton.setGeometry(QtCore.QRect(40, 50, 171, 20))
        self.runselectedRadioButton.setObjectName(_fromUtf8("runselectedRadioButton"))
        self.buttonGroup.addButton(self.runselectedRadioButton)
        self.layoutWidget1 = QtGui.QWidget(ProgConfigDlg)
        self.layoutWidget1.setGeometry(QtCore.QRect(168, 240, 213, 32))
        self.layoutWidget1.setObjectName(_fromUtf8("layoutWidget1"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.layoutWidget1)
        self.horizontalLayout.setMargin(0)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.cancelButton = QtGui.QPushButton(self.layoutWidget1)
        self.cancelButton.setObjectName(_fromUtf8("cancelButton"))
        self.horizontalLayout.addWidget(self.cancelButton)
        self.okayButton = QtGui.QPushButton(self.layoutWidget1)
        self.okayButton.setObjectName(_fromUtf8("okayButton"))
        self.horizontalLayout.addWidget(self.okayButton)

        self.retranslateUi(ProgConfigDlg)
        QtCore.QObject.connect(self.cancelButton, QtCore.SIGNAL(_fromUtf8("clicked()")), ProgConfigDlg.reject)
        QtCore.QMetaObject.connectSlotsByName(ProgConfigDlg)
        ProgConfigDlg.setTabOrder(self.runallRadioButton, self.runselectedRadioButton)
        ProgConfigDlg.setTabOrder(self.runselectedRadioButton, self.alignCheckBox)
        ProgConfigDlg.setTabOrder(self.alignCheckBox, self.consensusCheckBox)
        ProgConfigDlg.setTabOrder(self.consensusCheckBox, self.unmapCheckBox)
        ProgConfigDlg.setTabOrder(self.unmapCheckBox, self.hybridCheckBox)
        ProgConfigDlg.setTabOrder(self.hybridCheckBox, self.okayButton)
        ProgConfigDlg.setTabOrder(self.okayButton, self.cancelButton)

    def retranslateUi(self, ProgConfigDlg):
        ProgConfigDlg.setWindowTitle(_translate("ProgConfigDlg", "Project Settings", None))
        self.alignCheckBox.setText(_translate("ProgConfigDlg", "Barcode and primer identification", None))
        self.consensusCheckBox.setText(_translate("ProgConfigDlg", "Generate consensus sequences", None))
        self.unmapCheckBox.setText(_translate("ProgConfigDlg", "Dump unmapped reads", None))
        self.hybridCheckBox.setText(_translate("ProgConfigDlg", "Search for heterozygous locus", None))
        self.runallRadioButton.setText(_translate("ProgConfigDlg", "Run the whole process", None))
        self.runselectedRadioButton.setText(_translate("ProgConfigDlg", "Run selected functions", None))
        self.cancelButton.setText(_translate("ProgConfigDlg", "Cancel", None))
        self.okayButton.setText(_translate("ProgConfigDlg", "OK", None))

