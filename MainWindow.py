import sys
from PyQt5.QtCore import pyqtSignal, QThread, QSettings, QVariant, Qt, QSize, QPoint
from PyQt5.QtGui import QKeySequence, QIcon, QTextCursor
from PyQt5.QtWidgets import (
    QMainWindow, QMessageBox, QAction, QTreeWidgetItem, QDialog, QFileDialog, QTreeWidget,
    QTextEdit, QSplitter, QGridLayout, QWidget, QApplication
)

import pickle as pickle
from os.path import isfile
import gzip
import multiprocessing

# -homemade class
import ParameterDlg
import ProgConfigDlg
import ProjMergeDlg
import Parameters
import DataViewer
import ProjectEnviroment
import AboutDlg

__version__ = "0.2.0"


class MainWindow(QMainWindow):
    runinfo = pyqtSignal(object, object)

    def __init__(self, parent=None):
        ##------------
        #  main frame
        ##----------
        super(MainWindow, self).__init__(parent)
        self.paraTree = QTreeWidget()
        self.paraTree.setColumnCount(1)
        self.paraTree.setHeaderLabels([""])
        self.cmdTextBrowser = QTextEdit()
        self.cmdTextBrowser.setReadOnly(True)

        ###---right multifunction viewer---
        self.viewerDict = {}
        self.viewerLayout = QGridLayout()
        viewerWidget = QWidget(self)
        viewerWidget.setLayout(self.viewerLayout)
        # -----

        self.infoSplitter = QSplitter(Qt.Horizontal)
        self.infoSplitter.addWidget(self.paraTree)
        self.infoSplitter.addWidget(viewerWidget)
        self.mainSplitter = QSplitter(Qt.Vertical)
        self.mainSplitter.addWidget(self.infoSplitter)
        self.mainSplitter.addWidget(self.cmdTextBrowser)
        self.infoSplitter.setStretchFactor(0, 1)
        self.infoSplitter.setStretchFactor(1, 3)
        self.mainSplitter.setStretchFactor(0, 4)
        self.mainSplitter.setStretchFactor(1, 2)
        self.setCentralWidget(self.mainSplitter)

        self.status = self.statusBar()
        self.status.setSizeGripEnabled(False)

        # --connections---------
        self.paraTree.itemSelectionChanged.connect(self.treeselectionChanged)
        self.runinfo.connect(self.runinfohandle)

        ##--actions--------

        self.fileNewAction = self.createAction("&New Project", self.projNew,
                                               QKeySequence.New, "filenew", "Create a new project")
        self.fileOpenAction = self.createAction("&Open Project", self.projOpen,
                                                QKeySequence.Open, "fileopen",
                                                "Open an existing project")
        self.fileMergeAction = self.createAction("&Merge Projects", self.projMerge,
                                                 "Ctrl+M", "filemerge", "Merge existing projects")
        self.fileCloseAction = self.createAction("&Close Project", self.projClose,
                                                 QKeySequence.Close, "fileclose",
                                                 "Close current project")
        self.fileQuitAction = self.createAction("&Exit", self.closeApp,
                                                QKeySequence.Quit, "filequit",
                                                "Quit the application")

        self.projRunAction = self.createAction("&Run project", self.projRun,
                                               "Ctrl+U", "projrun", "Run current project")
        self.projStopAction = self.createAction("&Kill running job", self.projStop,
                                                "Ctrl+S", "projstop", "Stop current project")
        self.projSetupAction = self.createAction("&Project settings", self.projSetup,
                                                 "Ctrl+J", "projsetup",
                                                 "Settings for current project")

        # self.plotReadLength = self.createAction("Read &Length", self.plotReadLength,
        #                                        "Ctrl+L", None, "Plot read length distribution")
        #
        # self.plotAlignRatio = self.createAction("Align &Ratio", self.plotMapRatio,
        #                                          "Ctrl+R", None, "Plot ratios of aligned reads")
        # self.plotLocusLength = self.createAction("Locus Length", self.plotLocusLength,
        #                                         None, None, "Plot locus length distribution")
        self.hetSearchAction = self.createAction("Search Het", self.hetSearch,
                                                 None, None, "Search Hetero")
        self.aboutAction = self.createAction("About", self.about, None, None, "About")

        ##--Main Menu setup----
        self.fileMenu = self.menuBar().addMenu("&Project")
        self.fileMenu.addAction(self.fileNewAction)
        self.fileMenu.addAction(self.fileOpenAction)
        self.fileMenu.addAction(self.fileMergeAction)
        self.fileMenu.addAction(self.fileCloseAction)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.projSetupAction)
        self.fileMenu.addAction(self.projRunAction)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.fileQuitAction)

        # self.projMenu = self.menuBar().addMenu("&Job Setting")
        # self.projMenu.addAction(self.projSetupAction)
        # self.projMenu.addSeparator()
        # self.projMenu.addAction(self.projRunAction)

        self.projMenu = self.menuBar().addMenu("&Help")
        self.projMenu.addAction(self.aboutAction)
        # self.projMenu.addAction(self.hetSearchAction)
        # self.projMenu.addAction(self.projStopAction)

        ##--Toolbar setup----
        self.toolbar = self.addToolBar("main toolbar")
        self.toolbar.addAction(self.fileNewAction)
        self.toolbar.addAction(self.fileOpenAction)
        self.toolbar.addAction(self.fileMergeAction)
        self.toolbar.addAction(self.fileCloseAction)
        self.toolbar.addSeparator()
        self.toolbar.addAction(self.projSetupAction)
        self.toolbar.addSeparator()
        self.toolbar.addAction(self.projRunAction)
        # self.toolbar.addAction(self.projStopAction)

        # --------self status-----
        self.projenv = None
        self.dirty = False
        self.__initUi()
        self.status.showMessage("Ready!")
        self.__msgStatus = 1

        # -------------------------
        self.PROJECTFILE = "project.nma"

    def __initUi(self):
        settings = QSettings("CYSoft", "MLSTEZ")
        size = settings.value("MainWindow/Size", QVariant(QSize(850, 600)))
        self.resize(size)
        position = settings.value("MainWindow/Position", QVariant(QPoint(0, 0)))
        self.move(position)
        self.setWindowTitle("MLSTEZ")
        # self.restoreState(settings.value("MainWindow/State").toByteArray())

        self.__btninitStat()

        paras = Parameters.Parameters()
        paras.Filetype = settings.value("Paras/Filetype", "FASTQ")
        paras.MuscleCMD = settings.value("Paras/MuscleCMD", "")
        paras.ScoringSys = settings.value("Paras/ScoringSys", "phred33")
        paras.PadSeq = settings.value("Paras/PadSeq", "GGTAG")
        paras.UniPrimer = settings.value("Paras/UniPrimer", "CTGGAGCACGAGGACACTGA")
        paras.BarcodeLen = settings.value("Paras/BarcodeLen", 16)
        paras.MinReadNum = settings.value("Paras/MinReadNum", 3)
        paras.MaxReadNum = settings.value("Paras/MaxReadNum", 10)
        paras.FlankingLength = settings.value("Paras/FlankingLength", 5)
        paras.MatchScore = settings.value("Paras/MatchScore", 2)
        paras.MismatchScore = settings.value("Paras/MismatchScore", -1)
        paras.GapScore = settings.value("Paras/GapScore", -1)
        paras.MaxMisMatch = settings.value("Paras/MaxMisMatch", 3)
        paras.Threads = settings.value("Paras/Threads", 1)
        # paras.ConsensusCut = settings.value("Paras/ConsensusCut", 0.5).toFloat()[0]
        self.parameters = paras

        # ---renew viewer&tree status
        if self.dirty:
            for i in reversed(list(range(self.viewerLayout.count()))):
                self.viewerLayout.itemAt(i).widget().setParent(None)

            # for viewer in self.viewerDict:
            #    self.viewerLayout.removeWidget(self.viewerDict[viewer])
            self.viewerDict = {}
            self.paraTree.clear()

        self.treeRoot = self.createTreeitem("Project", "project", True)
        self.paraTree.addTopLevelItem(self.treeRoot)
        self.paraTree.expandAll()

        self.viewerCount = 0
        infoview = DataViewer.TextViewer()
        # infoview.setstyle("font: 20pt \"Courier\";")
        # welcomeInfo = "<div align='center'><B><I><font size=20> " \
        #              "<font color='orange'>MLSTEZ</font></font></I></B><p><Next-generation MLST solution></p></div>"
        welcomeInfo = """<p><div align='center'><em><strong><span style="font-size:20px;">MLST<span style="color:#ffa500;">EZ</span></span></strong></em></div></p>
<p><div align='center'>The next-generation MLST analyser</div></p>
<p><div align='center'><I>(ver. 0.2.0)</I></div></p>
<p><div align='center'>Created by: Yuan Chen</div></p>
<p><div align='center'>Duke University Medical Center</div></p>"""
        infoview.addtext(welcomeInfo)
        self.__addViewer(infoview, "MainInfo", True)

        self.cmdTextBrowser.clear()
        self.dirty = False  ##status of job in memory
        self.curRuncode = 0  ##running info code for current project

    ##--signal event----
    def treeselectionChanged(self):
        currItem = self.paraTree.currentItem()
        text = str(currItem.text(0))
        if (text not in self.viewerDict):
            return
        for viewer in self.viewerDict:
            if (viewer != text):
                self.viewerDict[viewer].hide()
            else:
                self.viewerDict[viewer].show()

    def closeEvent(self, event):
        settings = QSettings("CYSoft", "MLSTEZ")
        settings.setValue("MainWindow/Size", QVariant(self.size()))
        settings.setValue("MainWindow/Position", QVariant(self.pos()))
        # settings.setValue("MainWindow/State", QVariant(self.saveState()))
        self.__saveParatoSettings()

    def createAction(self, text, slot=None, shortcut=None, icon=None, tip=None, checkable=False,
                     signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            action.triggered.connect(slot)
        if checkable:
            action.setCheckable(True)
        return action

    def createTreeitem(self, name, icon=None, expanded=False):
        item = QTreeWidgetItem()
        item.setText(0, name)
        if icon is not None:
            item.setIcon(0, QIcon(":/%s.png" % icon))
        if expanded:
            item.setExpanded(True)
        return item

    def __saveParatoSettings(self):
        paras = self.parameters
        settings = QSettings("CYSoft", "MLSTEZ")
        settings.setValue("Paras/Filetype", QVariant(paras.Filetype))
        settings.setValue("Paras/MuscleCMD", QVariant(paras.MuscleCMD))
        settings.setValue("Paras/ScoringSys", QVariant(paras.ScoringSys))
        settings.setValue("Paras/PadSeq", QVariant(paras.PadSeq))
        settings.setValue("Paras/UniPrimer", QVariant(paras.UniPrimer))
        settings.setValue("Paras/BarcodeLen", QVariant(paras.BarcodeLen))
        settings.setValue("Paras/MinReadNum", QVariant(paras.MinReadNum))
        settings.setValue("Paras/MaxReadNum", QVariant(paras.MaxReadNum))
        settings.setValue("Paras/FlankingLength", QVariant(paras.FlankingLength))
        settings.setValue("Paras/MatchScore", QVariant(paras.MatchScore))
        settings.setValue("Paras/MismatchScore", QVariant(paras.MismatchScore))
        settings.setValue("Paras/GapScore", QVariant(paras.GapScore))
        settings.setValue("Paras/MaxMisMatch", QVariant(paras.MaxMisMatch))
        settings.setValue("Paras/Threads", QVariant(paras.Threads))
        # settings.setValue("Paras/ConsensusCut", QVariant(paras.ConsensusCut))

    def __btninitStat(self):
        self.fileNewAction.setEnabled(True)
        self.fileOpenAction.setEnabled(True)
        self.fileCloseAction.setEnabled(False)
        self.fileMergeAction.setEnabled(True)
        self.projSetupAction.setEnabled(False)
        self.projRunAction.setEnabled(False)
        self.projStopAction.setEnabled(False)

    def __btnconfigStat(self):
        self.fileNewAction.setEnabled(True)
        self.fileOpenAction.setEnabled(True)
        self.fileCloseAction.setEnabled(True)
        self.fileMergeAction.setEnabled(True)
        self.projSetupAction.setEnabled(True)
        self.projRunAction.setEnabled(False)
        self.projStopAction.setEnabled(False)

    def __btnrunStat(self):
        self.fileNewAction.setEnabled(False)
        self.fileOpenAction.setEnabled(False)
        self.fileCloseAction.setEnabled(False)
        self.fileMergeAction.setEnabled(False)
        self.projSetupAction.setEnabled(False)
        self.projRunAction.setEnabled(False)
        self.projStopAction.setEnabled(True)

        ##--Action controls-----

    def projNew(self):
        # code
        if self.dirty:
            reply = QMessageBox.question(self, "Project", "Close current project?",
                                         QMessageBox.Yes | QMessageBox.Default,
                                         QMessageBox.No | QMessageBox.Escape)
            if (reply == QMessageBox.Yes):
                self.__initUi()
            else:
                return

        paraDlg = ParameterDlg.ParameterDlg(self.parameters)
        retDlgCode = paraDlg.exec_()
        if (retDlgCode == QDialog.Accepted):
            self.parameters = paraDlg.parameters
            self.parameterfile = self.parameters.Out_Folder + "/config.ini"
            self.projenv = None
            isokay, errmsg = self.parameters.savefile(self.parameterfile)
            if (isokay):
                self.__showParameter()
            else:
                QMessageBox.warning(self, "Configure file save error",
                                    ("Error: %s" % (errmsg)), QMessageBox.Ok | QMessageBox.Default)

    def projOpen(self):
        if (self.dirty):
            if self.dirty:
                reply = QMessageBox.question(self, "Project", "Close current project?",
                                             QMessageBox.Yes | QMessageBox.Default,
                                             QMessageBox.No | QMessageBox.Escape)
                if (reply == QMessageBox.Yes):
                    self.__initUi()
                else:
                    return

        foldername = str(QFileDialog.getExistingDirectory(self, 'Project Folder',
                                                          ""))
        if (foldername):
            projfile = foldername + "/" + self.PROJECTFILE
            if (not isfile(projfile)):
                QMessageBox.warning(self, "Error",
                                    ("Cannot file project information"),
                                    QMessageBox.Ok | QMessageBox.Default)
                return
            self.projenv = ProjectEnviroment.ProjectEnviroment(self.parameters, self.runinfo)
            try:
                fh = gzip.open(projfile, "rb")
                buffer = ""
                BLOCK_SIZE = 2 ** 20
                while True:
                    data = fh.read(BLOCK_SIZE)
                    if data == "":
                        break
                    buffer += data
                projinfo = pickle.loads(buffer)
                self.projenv.parameters = projinfo['parameters']
                self.projenv.Seqs = projinfo['Seqs']
                self.projenv.Primers = projinfo['Primers']
                self.projenv.Barcodes = projinfo['Barcodes']
                self.projenv.AlignedSeqs = projinfo['AlignedSeqs']
                self.projenv.SortedSeqs = projinfo['SortedSeqs']
                self.projenv.locusLengthRange = projinfo['locusLengthRange']
                self.projenv.consSeqs = projinfo['consSeqs']
                self.projenv.locusLengthsInfo = projinfo['locusLengthsInfo']
                self.projenv.num_unbarcode = projinfo['num_unbarcode']
                self.projenv.num_unprimer = projinfo['num_unprimer']
                self.projenv.StrainStats = projinfo['StrainStats']
                self.projenv.HetSeqs = projinfo['HetSeqs']
                self.projenv.HetInfo = projinfo['HetInfo']
                self.projenv.HetStats = projinfo['HetStats']

                self.curRuncode = projinfo['runstats']
                self.dirty = True
                self.parameters = self.projenv.parameters
                self.parameterfile = self.projenv.parameters.Out_Folder + "/config.ini"
                self.__showParameter()
                self.__showSummary()
                self.__showSeqLength()
                self.__showAligned()
                self.__btnconfigStat()
                if (self.curRuncode & 1 << 1):
                    self.__showConsus()
                if (self.curRuncode & 1 << 3):
                    self.__showHet()
            except Exception as e:
                QMessageBox.warning(self.mainframe, "Error", ("Error: %s" % e),
                                    QMessageBox.Ok | QMessageBox.Default)
        return

    def projMerge(self):
        if (self.dirty):
            if self.dirty:
                reply = QMessageBox.question(self, "Project", "Close current project?",
                                             QMessageBox.Yes | QMessageBox.Default,
                                             QMessageBox.No | QMessageBox.Escape)
                if (reply == QMessageBox.Yes):
                    self.__initUi()
                else:
                    return
        projMergeDlg = ProjMergeDlg.ProjMergeDlg(self.parameters.Out_Folder)
        retDlgCode = projMergeDlg.exec_()
        if (retDlgCode == QDialog.Accepted):
            projFiles = projMergeDlg.projFiles
            projName = projMergeDlg.projName
            projFolder = projMergeDlg.outFolder
            self.projenv = ProjectEnviroment.ProjectEnviroment(self.parameters, self.runinfo)
            for projfile in projFiles:
                projinfo = pickle.load(open(projfile, "rb"))
                self.projenv.parameters = projinfo['parameters']
                self.projenv.Seqs += projinfo['Seqs']
                self.projenv.Primers += projinfo['Primers']
                self.projenv.Barcodes += projinfo['Barcodes']
                self.projenv.AlignedSeqs += projinfo['AlignedSeqs']
                self.projenv.num_unbarcode += projinfo['num_unbarcode']
                self.projenv.num_unprimer += projinfo['num_unprimer']
            self.projenv.parameters.ProjectName = projName
            self.projenv.parameters.Out_Folder = projFolder
            self.projenv.locusLengths()
            self.projenv.sortAlignedSeqs()
            self.projenv.strainStats()
            self.treeRoot.setText(0, ("Project - %s" % (self.parameters.ProjectName)))
            self.__showSummary()
            self.__showSeqLength()
            self.__showAligned()
            self.curRuncode = int('0001', 2)
            self.saveStats()
            self.__btnconfigStat()
            self.dirty = True

    def projClose(self):
        reply = QMessageBox.question(self, "Project", "Close current project?",
                                     QMessageBox.Yes | QMessageBox.Default,
                                     QMessageBox.No | QMessageBox.Escape)
        if (reply == QMessageBox.Yes):
            self.__initUi()
        else:
            return

    def projRun(self):
        if (self.projenv is None):
            self.projenv = ProjectEnviroment.ProjectEnviroment(self.parameters, self.runinfo)
        self.threadpool = []
        compRuncode = ~ self.curRuncode
        needRun = self.jobcode & compRuncode

        runproj = ProjectThread(self, needRun)
        self.threadpool.append(runproj)
        runproj.jobstats.connect(self.jobstats)
        # self.connect(runproj, SIGNAL("JobDone"), self.jobstats)
        runproj.start()
        self.__btnrunStat()

    def projStop(self):
        # if(self.projenv.ismultirun):
        #    self.projenv.alignStop()
        for thread in self.threadpool:
            if (thread.isRunning()):
                thread.terminate()
                thread.wait()
        self.threadpool = []
        self.__btnconfigStat()
        self.showMsg("Running job has been terminated!")

    def projSetup(self):
        confDlg = ProgConfigDlg.ProgConfigDlg(self.curRuncode)
        confDlg.setWindowTitle("Job settings...")
        retDlgCode = confDlg.exec_()
        if (retDlgCode == QDialog.Accepted):
            self.jobcode = confDlg.retcode
            self.projRunAction.setEnabled(True)

    def hetSearch(self):
        self.projenv.HetSearch()
        self.projenv.hetStats()
        newTreeItem = self.createTreeitem("Het Stats", "table")
        self.treeRoot.addChild(newTreeItem)
        statsview = DataViewer.HetViewer(self.projenv.parameters.Out_Folder,
                                         self.projenv.HetSeqs)
        statsview.showData(self.projenv.HetStats)
        self.__addViewer(statsview, "Het Stats")
        self.saveStats()

    def about(self):
        aboutDlg = AboutDlg.AboutDlg()
        aboutDlg.setFixedSize(450, 280)
        aboutDlg.setWindowTitle("About MLSTEZ")
        aboutDlg.exec_()

    def closeApp(self):
        # code
        self.status.showMessage("Quit")
        self.close()

    def saveStats(self):
        projinfo = {}
        projinfo['parameters'] = self.projenv.parameters
        projinfo['Seqs'] = self.projenv.Seqs
        projinfo['Primers'] = self.projenv.Primers
        projinfo['Barcodes'] = self.projenv.Barcodes
        projinfo['AlignedSeqs'] = self.projenv.AlignedSeqs
        projinfo['SortedSeqs'] = self.projenv.SortedSeqs
        projinfo['locusLengthRange'] = self.projenv.locusLengthRange
        projinfo['consSeqs'] = self.projenv.consSeqs
        projinfo['num_unbarcode'] = self.projenv.num_unbarcode
        projinfo['num_unprimer'] = self.projenv.num_unprimer
        projinfo['locusLengthsInfo'] = self.projenv.locusLengthsInfo
        projinfo['StrainStats'] = self.projenv.StrainStats
        projinfo['HetSeqs'] = self.projenv.HetSeqs
        projinfo['HetInfo'] = self.projenv.HetInfo
        projinfo['HetStats'] = self.projenv.HetStats
        projinfo['runstats'] = self.curRuncode

        infofile = self.projenv.parameters.Out_Folder + "/" + self.PROJECTFILE
        try:
            fh = gzip.open(infofile, "wb")
            pickle.dump(projinfo, fh)
        except Exception as e:
            QMessageBox.warning(self.mainframe, "Error", ("Error: %s" % e),
                                QMessageBox.Ok | QMessageBox.Default)
        return

    def runinfohandle(self, msg, end):
        if (end is None):
            self.showMsg(msg)
        else:
            self.showMsg(msg, end)

    def jobstats(self, statcode, msg):
        if (statcode == 1):
            if (msg == "Loaded"):
                self.__showSummary()
                self.__showSeqLength()
            if (msg == "Aligned"):
                self.__showAligned()
            if (msg == "Consed"):
                self.__showConsus()
            if (msg == "Heted"):
                self.__showHet()
            if (msg == "Done"):
                self.__btnconfigStat()
                QMessageBox.information(self, "Job finished", "All jobs have finished!")
        # print data
        # self.list_widget.addItem(unicode(data))

    def __addViewer(self, viewer, name, isShow=False):
        if (not isShow):
            viewer.hide()
        self.viewerDict[name] = viewer
        self.viewerLayout.addWidget(viewer, self.viewerCount, 0)
        self.viewerCount += 1

    def __showParameter(self):
        self.treeRoot.setText(0, ("Project - %s" % (self.parameters.ProjectName)))
        self.__btnconfigStat()
        if (isfile(self.parameterfile)):
            newTreeItem = self.createTreeitem("Parameters", "info")
            self.treeRoot.addChild(newTreeItem)
            infoview = DataViewer.TextViewer()
            infoview.readfile(self.parameterfile)
            self.__addViewer(infoview, "Parameters")
            self.__saveParatoSettings()
        self.dirty = True

    def __showSummary(self):
        primercount = self.projenv.numprimers()
        barcodecount = self.projenv.numbarcodes()
        seqcount = self.projenv.numseqs()

        newTreeItem = self.createTreeitem("Summary", "info")
        self.treeRoot.addChild(newTreeItem)
        self.infoview = DataViewer.TextViewer()
        infotext = (
                    "Project: %s\n\nNumber of barcodes: %s\n\nNumber of primer pairs: %s\n\nNumber of reads: %s"
                    % (self.projenv.parameters.ProjectName, barcodecount, primercount, seqcount))

        self.infoview.addtext(infotext)
        self.__addViewer(self.infoview, "Summary")

    def __showSeqLength(self):
        seqlengths = self.projenv.seqlengths()
        newTreeItem = self.createTreeitem("Read length dist", "chart")
        self.treeRoot.addChild(newTreeItem)
        infoview = DataViewer.ImageViewer(self.projenv.parameters.Out_Folder)
        infoview.plotDensity(seqlengths)
        self.__addViewer(infoview, "Read length dist")

    def __showAligned(self):
        # ---viewer for align ratio
        readcount = self.projenv.numseqs()
        unaligned = self.projenv.num_unbarcode
        barcodealigned = readcount - unaligned
        unprimer = self.projenv.num_unprimer
        bothaligned = barcodealigned - unprimer
        alignfracs = [bothaligned, unprimer, unaligned]
        newTreeItem = self.createTreeitem("Alignment ratio", "chart")
        self.treeRoot.addChild(newTreeItem)
        ratioview = DataViewer.ImageViewer(self.projenv.parameters.Out_Folder)
        ratioview.plotAlignRatio(alignfracs)
        self.__addViewer(ratioview, "Alignment ratio")

        # ---viewer for length distribution
        newTreeItem = self.createTreeitem("Length ranges", "chart")
        self.treeRoot.addChild(newTreeItem)
        lengthview = DataViewer.ImageViewer(self.projenv.parameters.Out_Folder)
        lengthview.plotLengthRange(self.projenv.locusLengthsInfo)
        self.__addViewer(lengthview, "Length ranges")

        # ---viewer for data table----
        newTreeItem = self.createTreeitem("Reads stats", "table")
        self.treeRoot.addChild(newTreeItem)
        statsview = DataViewer.TableViewer(self.projenv.parameters.Out_Folder,
                                           self.projenv.AlignedSeqs, self.projenv.locusLengthRange)
        statsview.showData(self.projenv.StrainStats, self.projenv.parameters.MinReadNum)
        self.__addViewer(statsview, "Reads stats")

        # ---update infoview
        infotext = ("\nBarcode Alignment: %s/%s reads have barcodes aligned\n\n"
                    % (barcodealigned, readcount))
        infotext += ("Primer Alignment: %s/%s reads have primers aligned"
                     % (barcodealigned - unprimer, barcodealigned))
        self.infoview.addtext(infotext)

    def __showConsus(self):
        newTreeItem = self.createTreeitem("Consensus", "table")
        self.treeRoot.addChild(newTreeItem)
        consview = DataViewer.ConsensusViewer(self.projenv.consSeqs, self.projenv.Primers,
                                              self.projenv.Barcodes,
                                              self.projenv.parameters.Out_Folder)
        consview.showData()
        self.__addViewer(consview, "Consensus")

    def __showHet(self):
        newTreeItem = self.createTreeitem("Heterozygous loci", "table")
        self.treeRoot.addChild(newTreeItem)
        statsview = DataViewer.HetViewer(self.projenv.parameters.Out_Folder,
                                         self.projenv.HetSeqs)
        statsview.showData(self.projenv.HetStats)
        self.__addViewer(statsview, "Heterozygous loci")

    ##msg port project enviroment
    def showMsg(self, message, end='\n'):
        cursor = self.cmdTextBrowser
        if (self.__msgStatus == 1):
            self.cmdTextBrowser.append(message)
        else:
            self.cmdTextBrowser.moveCursor(QTextCursor.End)
            self.cmdTextBrowser.insertPlainText(message)
            self.cmdTextBrowser.moveCursor(QTextCursor.End)
        if (end == '\n'):
            self.cmdTextBrowser.append("")
            self.__msgStatus = 1
        else:
            self.__msgStatus = 0


class ProjectThread(QThread):
    jobstats = pyqtSignal(int, object)

    def __init__(self, mainframe, jobcode):
        QThread.__init__(self)
        self.jobcode = jobcode
        self.mainframe = mainframe
        self.projenv = mainframe.projenv

    def run(self):
        # self.projenv.msgHandle.showMsg("Start running program...")
        self.projenv.showMsg("Start running program...")
        # t0 = time()
        module0 = [self.__proj_loadFiles, self.__proj_alignSeqs]
        module1 = [self.__proj__genCons]
        module2 = [self.projenv.DumpUnmappedReads]
        module3 = [self.__proj_searchHet]
        runmodule = [module0, module1, module2, module3]
        for i in range(4):
            mask = 1 << i
            if (mask & self.jobcode):
                for job in runmodule[i]:
                    isokay, error = job()
                    if (not isokay):
                        # print ("Error: %s" %error)
                        QMessageBox.warning(self.mainframe, "Error", ("Error: %s" % error),
                                            QMessageBox.Ok | QMessageBox.Default)
                        return
                self.mainframe.curRuncode = self.mainframe.curRuncode + (1 << i)
        self.mainframe.saveStats()

        self.jobstats.emit(1, "Done")

    def __proj_loadFiles(self):
        if (not (self.projenv.status & 1)):
            isokay, error = self.projenv.loadFiles()
            if (isokay):
                self.jobstats.emit(1, "Loaded")
                return (True, None)
            else:
                return (False, error)

    def __proj_alignSeqs(self):
        if (not (self.projenv.status & (1 << 2))):
            isokay, error = self.projenv.alignSeqs()
            if (isokay):
                self.projenv.locusLengths()
                self.projenv.sortAlignedSeqs()
                self.projenv.strainStats()
                self.jobstats.emit(1, "Aligned")
                return (True, None)
            else:
                return (False, error)

    def __proj__genCons(self):
        if (not (self.projenv.status & (1 << 4))):
            isokay, error = self.projenv.GenerateConsensus()
            if (isokay):
                self.jobstats.emit(1, "Consed")
                return (True, None)
            else:
                return (False, error)

    def __proj_searchHet(self):
        isokay, error = self.projenv.HetSearch()
        if (isokay):
            self.projenv.hetStats()
            self.jobstats.emit(1, "Heted")
            return (True, None)
        else:
            return (False, error)


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("MLSTEZ")
    mainwindow = MainWindow()
    mainwindow.show()
    app.exec_()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
