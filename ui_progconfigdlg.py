# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'progconfigdlg.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_ProgConfigDlg(object):
    def setupUi(self, ProgConfigDlg):
        ProgConfigDlg.setObjectName("ProgConfigDlg")
        ProgConfigDlg.resize(415, 290)
        self.frame = QtWidgets.QFrame(ProgConfigDlg)
        self.frame.setGeometry(QtCore.QRect(40, 80, 341, 151))
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.layoutWidget = QtWidgets.QWidget(self.frame)
        self.layoutWidget.setGeometry(QtCore.QRect(10, 10, 321, 131))
        self.layoutWidget.setObjectName("layoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.layoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.alignCheckBox = QtWidgets.QCheckBox(self.layoutWidget)
        self.alignCheckBox.setObjectName("alignCheckBox")
        self.verticalLayout.addWidget(self.alignCheckBox)
        self.consensusCheckBox = QtWidgets.QCheckBox(self.layoutWidget)
        self.consensusCheckBox.setObjectName("consensusCheckBox")
        self.verticalLayout.addWidget(self.consensusCheckBox)
        self.unmapCheckBox = QtWidgets.QCheckBox(self.layoutWidget)
        self.unmapCheckBox.setObjectName("unmapCheckBox")
        self.verticalLayout.addWidget(self.unmapCheckBox)
        self.hybridCheckBox = QtWidgets.QCheckBox(self.layoutWidget)
        self.hybridCheckBox.setObjectName("hybridCheckBox")
        self.verticalLayout.addWidget(self.hybridCheckBox)
        self.runallRadioButton = QtWidgets.QRadioButton(ProgConfigDlg)
        self.runallRadioButton.setGeometry(QtCore.QRect(40, 20, 171, 20))
        self.runallRadioButton.setObjectName("runallRadioButton")
        self.buttonGroup = QtWidgets.QButtonGroup(ProgConfigDlg)
        self.buttonGroup.setObjectName("buttonGroup")
        self.buttonGroup.addButton(self.runallRadioButton)
        self.runselectedRadioButton = QtWidgets.QRadioButton(ProgConfigDlg)
        self.runselectedRadioButton.setGeometry(QtCore.QRect(40, 50, 171, 20))
        self.runselectedRadioButton.setObjectName("runselectedRadioButton")
        self.buttonGroup.addButton(self.runselectedRadioButton)
        self.layoutWidget1 = QtWidgets.QWidget(ProgConfigDlg)
        self.layoutWidget1.setGeometry(QtCore.QRect(168, 240, 213, 32))
        self.layoutWidget1.setObjectName("layoutWidget1")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.layoutWidget1)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.cancelButton = QtWidgets.QPushButton(self.layoutWidget1)
        self.cancelButton.setObjectName("cancelButton")
        self.horizontalLayout.addWidget(self.cancelButton)
        self.okayButton = QtWidgets.QPushButton(self.layoutWidget1)
        self.okayButton.setObjectName("okayButton")
        self.horizontalLayout.addWidget(self.okayButton)

        self.retranslateUi(ProgConfigDlg)
        self.cancelButton.clicked.connect(ProgConfigDlg.reject)
        QtCore.QMetaObject.connectSlotsByName(ProgConfigDlg)
        ProgConfigDlg.setTabOrder(self.runallRadioButton, self.runselectedRadioButton)
        ProgConfigDlg.setTabOrder(self.runselectedRadioButton, self.alignCheckBox)
        ProgConfigDlg.setTabOrder(self.alignCheckBox, self.consensusCheckBox)
        ProgConfigDlg.setTabOrder(self.consensusCheckBox, self.unmapCheckBox)
        ProgConfigDlg.setTabOrder(self.unmapCheckBox, self.hybridCheckBox)
        ProgConfigDlg.setTabOrder(self.hybridCheckBox, self.okayButton)
        ProgConfigDlg.setTabOrder(self.okayButton, self.cancelButton)

    def retranslateUi(self, ProgConfigDlg):
        _translate = QtCore.QCoreApplication.translate
        ProgConfigDlg.setWindowTitle(_translate("ProgConfigDlg", "Project Settings"))
        self.alignCheckBox.setText(_translate("ProgConfigDlg", "Barcode and primer identification"))
        self.consensusCheckBox.setText(_translate("ProgConfigDlg", "Generate consensus sequences"))
        self.unmapCheckBox.setText(_translate("ProgConfigDlg", "Dump unmapped reads"))
        self.hybridCheckBox.setText(_translate("ProgConfigDlg", "Search for heterozygous locus"))
        self.runallRadioButton.setText(_translate("ProgConfigDlg", "Run the whole process"))
        self.runselectedRadioButton.setText(_translate("ProgConfigDlg", "Run selected functions"))
        self.cancelButton.setText(_translate("ProgConfigDlg", "Cancel"))
        self.okayButton.setText(_translate("ProgConfigDlg", "OK"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    ProgConfigDlg = QtWidgets.QDialog()
    ui = Ui_ProgConfigDlg()
    ui.setupUi(ProgConfigDlg)
    ProgConfigDlg.show()
    sys.exit(app.exec_())

