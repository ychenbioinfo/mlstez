# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'projmergedlg.ui'
#
# Created: Fri Apr 11 10:58:03 2014
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

class Ui_ProjMergeDlg(object):
    def setupUi(self, ProjMergeDlg):
        ProjMergeDlg.setObjectName(_fromUtf8("ProjMergeDlg"))
        ProjMergeDlg.resize(424, 285)
        self.groupBox = QtGui.QGroupBox(ProjMergeDlg)
        self.groupBox.setGeometry(QtCore.QRect(13, 95, 398, 141))
        self.groupBox.setTitle(_fromUtf8(""))
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.gridLayout = QtGui.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.addButton = QtGui.QPushButton(self.groupBox)
        self.addButton.setObjectName(_fromUtf8("addButton"))
        self.gridLayout.addWidget(self.addButton, 0, 1, 1, 1)
        self.removeButton = QtGui.QPushButton(self.groupBox)
        self.removeButton.setObjectName(_fromUtf8("removeButton"))
        self.gridLayout.addWidget(self.removeButton, 1, 1, 1, 1)
        self.projsListWidget = QtGui.QListWidget(self.groupBox)
        self.projsListWidget.setObjectName(_fromUtf8("projsListWidget"))
        self.gridLayout.addWidget(self.projsListWidget, 0, 0, 3, 1)
        self.layoutWidget = QtGui.QWidget(ProjMergeDlg)
        self.layoutWidget.setGeometry(QtCore.QRect(13, 46, 323, 33))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.layoutWidget)
        self.horizontalLayout.setMargin(0)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.label = QtGui.QLabel(self.layoutWidget)
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout.addWidget(self.label)
        self.outfolderEdit = QtGui.QLineEdit(self.layoutWidget)
        self.outfolderEdit.setObjectName(_fromUtf8("outfolderEdit"))
        self.horizontalLayout.addWidget(self.outfolderEdit)
        self.selectButton = QtGui.QPushButton(self.layoutWidget)
        self.selectButton.setObjectName(_fromUtf8("selectButton"))
        self.horizontalLayout.addWidget(self.selectButton)
        self.widget = QtGui.QWidget(ProjMergeDlg)
        self.widget.setGeometry(QtCore.QRect(13, 13, 240, 23))
        self.widget.setObjectName(_fromUtf8("widget"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.widget)
        self.horizontalLayout_2.setMargin(0)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.label_2 = QtGui.QLabel(self.widget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.horizontalLayout_2.addWidget(self.label_2)
        self.projnameEdit = QtGui.QLineEdit(self.widget)
        self.projnameEdit.setObjectName(_fromUtf8("projnameEdit"))
        self.horizontalLayout_2.addWidget(self.projnameEdit)
        self.widget1 = QtGui.QWidget(ProjMergeDlg)
        self.widget1.setGeometry(QtCore.QRect(200, 242, 213, 32))
        self.widget1.setObjectName(_fromUtf8("widget1"))
        self.horizontalLayout_3 = QtGui.QHBoxLayout(self.widget1)
        self.horizontalLayout_3.setMargin(0)
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem)
        self.cancelButton = QtGui.QPushButton(self.widget1)
        self.cancelButton.setObjectName(_fromUtf8("cancelButton"))
        self.horizontalLayout_3.addWidget(self.cancelButton)
        self.okayButton = QtGui.QPushButton(self.widget1)
        self.okayButton.setObjectName(_fromUtf8("okayButton"))
        self.horizontalLayout_3.addWidget(self.okayButton)
        self.label.setBuddy(self.outfolderEdit)
        self.label_2.setBuddy(self.projnameEdit)

        self.retranslateUi(ProjMergeDlg)
        QtCore.QObject.connect(self.cancelButton, QtCore.SIGNAL(_fromUtf8("clicked()")), ProjMergeDlg.reject)
        QtCore.QMetaObject.connectSlotsByName(ProjMergeDlg)
        ProjMergeDlg.setTabOrder(self.projnameEdit, self.outfolderEdit)
        ProjMergeDlg.setTabOrder(self.outfolderEdit, self.selectButton)
        ProjMergeDlg.setTabOrder(self.selectButton, self.addButton)
        ProjMergeDlg.setTabOrder(self.addButton, self.removeButton)
        ProjMergeDlg.setTabOrder(self.removeButton, self.okayButton)
        ProjMergeDlg.setTabOrder(self.okayButton, self.cancelButton)
        ProjMergeDlg.setTabOrder(self.cancelButton, self.projsListWidget)

    def retranslateUi(self, ProjMergeDlg):
        ProjMergeDlg.setWindowTitle(_translate("ProjMergeDlg", "Dialog", None))
        self.addButton.setText(_translate("ProjMergeDlg", "Add Project", None))
        self.removeButton.setText(_translate("ProjMergeDlg", "Remove", None))
        self.label.setText(_translate("ProjMergeDlg", "Project Folder:", None))
        self.selectButton.setText(_translate("ProjMergeDlg", "Select", None))
        self.label_2.setText(_translate("ProjMergeDlg", "Project Name:", None))
        self.cancelButton.setText(_translate("ProjMergeDlg", "Cancel", None))
        self.okayButton.setText(_translate("ProjMergeDlg", "OK", None))

