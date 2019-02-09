from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from scipy.stats import gaussian_kde
from numpy import arange
import matplotlib.mlab as mlab
#import numpy as np
#import scipy.stats as stats
import matplotlib.pyplot as plt
import random
import sys
import SeqViewDlg

class ImageViewer(QWidget):
    def __init__(self, outfolder="", parent = None):
        super(ImageViewer, self).__init__(parent)
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.button = QPushButton('Save Image')
        buttonlayout = QHBoxLayout()
        buttonlayout.addStretch()
        buttonlayout.addWidget(self.button)
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addLayout(buttonlayout)
        self.setLayout(layout)
        #self.plot()
        self.button.clicked.connect(self.saveImage)
        self.outfolder = outfolder
    
    def saveImage(self):
        filename = str(QFileDialog.getSaveFileName(
            self, "Save File", self.outfolder, "Images (*.png *.bmp)"))
        if(filename):
            fileN = filename.split('.')
            if(fileN[-1] != "png" and fileN !="bmp"):
                QMessageBox.warning(self, "Output file format error:",
                                "Image file must be 'png' or 'bmp' file!")
            else:
                try:
                    self.canvas.print_figure(filename)
                    QMessageBox.information(self, "Image file saved", "Image file saved!")
                except IOError as e:
                    QMessageBox.warning(self, "Error:",
                                "Error: %s" %e.strerror, QMessageBox.Ok|QMessageBox.Default)
    
    def plotDensity(self, data):
        ax = self.figure.add_subplot(111)
        ax.hist(data, bins = 20)

        ax.set_title("Read Length Distribution")
        ax.set_xlabel("Read Length")
        ax.set_ylabel("Count")
        self.canvas.draw()
    
    def plotLengthRange(self, seqlengthinfo):
        locusnames = seqlengthinfo['name']
        data = seqlengthinfo['data']
        ax = self.figure.add_subplot(111)
        
        colnum = len(data)
        pos = list(range(colnum))
        w = min(0.15*max(colnum,1.0),0.5)
        for i in range(colnum):
            k = gaussian_kde(data[i])
            m = k.dataset.min()
            M = k.dataset.max()
            x = arange(m,M,(M-m)/100.)
            v = k.evaluate(x)
            v = v/v.max()*w
            ax.fill_betweenx(x,i,v+i,facecolor='y',alpha=0.3)
            ax.fill_betweenx(x,i,-v+i,facecolor='y',alpha=0.3)
        ax.boxplot(data,notch=1,positions=pos,vert=1)
        #ax.boxplot(data,'gD')
        ax.set_xticklabels(locusnames,rotation=40)
        ax.set_title("Length distribution of loci")
        self.canvas.draw()
    
    def plotAlignRatio(self, data):
        ax = self.figure.add_subplot(111)
        #data1 = [100,10,5]
        labels = ['Both Aligned','Barcode Only','Unaligned']
        explode=(0.05, 0, 0)
        #ax.hold(False)
        def my_autopct(pct):
            total=sum(data)
            val=int(pct*total/100.0)
            return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
        
        ax.pie(data, labels = labels, explode=explode, autopct=my_autopct, shadow=True)
        ax.set_title("Alignment ratio")
        #ax.set_title("Read alignment ratio")
        self.canvas.draw()
    

class TextViewer(QWidget):
    def __init__(self, parent=None):
        super(TextViewer, self).__init__(parent)
        self.viewer = QTextEdit()
        self.viewer.setReadOnly(True)
        layout = QHBoxLayout()
        layout.addWidget(self.viewer)
        self.setLayout(layout)
    
    def readfile(self, fn):
        fh = open(fn, 'r')
        content = fh.read()
        self.viewer.append(content)
    
    def addtext(self, text):
        self.viewer.append(text)

    def clear(self):
        self.viewer.clear()
        
    def setstyle(self, style):
        self.viewer.setStyleSheet(style)


class ConsensusViewer(QWidget):
    def __init__(self, consSeqs=None, Primers=None, Barcodes=None, outfolder = "", parent=None):
        super(ConsensusViewer, self).__init__(parent)
        self.outfolder = outfolder
        self.viewer = QTableWidget()
        self.button = QPushButton('Export table')
        buttonlayout = QHBoxLayout()
        buttonlayout.addStretch()
        buttonlayout.addWidget(self.button)
        layout = QVBoxLayout()
        layout.addWidget(self.viewer)
        layout.addLayout(buttonlayout)
        self.setLayout(layout)
        self.consSeqs = consSeqs
        self.strainIds = self._UniqueIDs(Barcodes,1)
        self.geneIds = self._UniqueIDs(Primers,2)
        #self.strainIds = ["H99"]
        #self.geneIds = ["LAC1","CAP59","SOD1","PLB1","TEF1","URA5","IGS1","GPD1"]
        self.button.clicked.connect(self.exportData)
        self.viewer.doubleClicked.connect(self.showSeqs)
    
    def showData(self):
        strainCons = {}
        for gene in self.consSeqs:
            geneSeqs = self.consSeqs[gene]
            for seq in geneSeqs:
                strainid = seq.id
                if(strainid not in strainCons):
                    tempstrain = {}
                    tempstrain[gene] = seq
                    strainCons[strainid] = tempstrain
                else:
                    tempstrain = strainCons[strainid]
                    tempstrain[gene] = seq
                    strainCons[strainid] = tempstrain
                    
        colnum = len(self.geneIds) + 1
        rownum = len(self.strainIds)
        self.viewer.setRowCount(rownum)
        self.viewer.setColumnCount(colnum)
        self.viewer.setShowGrid(True)
        self.viewer.setHorizontalHeaderItem(0, QTableWidgetItem("Strain"))
        for i in range(1,colnum):
            self.viewer.setHorizontalHeaderItem(i, QTableWidgetItem(self.geneIds[i-1]))
            
        rowcount = 0
        for sid in self.strainIds:
            colcount = 0
            isInfo = True
            strainSeqs = {}
            item = QTableWidgetItem(sid)
            if (sid in strainCons):
                strainSeqs = strainCons[sid]
                #item.setBackground(QColor(249,148,46))
            item.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEnabled)   
            self.viewer.setItem(rowcount,colcount,QTableWidgetItem(item))
            for gene in self.geneIds:
                colcount += 1
                item = None
                if(gene in strainSeqs):
                    item = QTableWidgetItem("Y")
                    #item.setBackground(QColor(249,148,46))
                if(item is None):
                    item = QTableWidgetItem("NA")
                    item.setBackground(QColor(170,170,170))
                item.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEnabled)
                self.viewer.setItem(rowcount,colcount,QTableWidgetItem(item))
            rowcount += 1
        self.viewer.resizeColumnsToContents()
        self.viewer.horizontalHeader().setStretchLastSection(True)
        self.strainCons = strainCons

    def exportData(self):
        import csv
        filename = str(QFileDialog.getSaveFileName(
            self, "Export data", self.outfolder, "CSV file (*.csv)"))
        if(filename):
            fileN = filename.split('.')
            if(fileN[-1] != "csv"):
                filename += ".csv"
            try:
                fh_outfile = open(filename,'w')
                csv_writer = csv.writer(fh_outfile, dialect='excel')
                header = self.geneIds
                header.insert(0,'Strain')
                csv_writer.writerow(header)
                for i in range(self.viewer.rowCount()):
                    content = []
                    for j in range(self.viewer.columnCount()):
                        content.append(self.viewer.item(i,j).text())
                    csv_writer.writerow(content)
                
                fh_outfile.close()
                QMessageBox.information(self, "Data Exported", "Data exported!")
            except IOError as e:
                QMessageBox.warning(self, "Error:",
                                    "Error: %s" %e.strerror, QMessageBox.Ok|QMessageBox.Default)

    def showSeqs(self,index):
        rownum = index.row()
        colnum = index.column()
        strain = self.strainIds[rownum]
        outseq = ""
        windowtitle = ""
        if(colnum == 0):
            if(strain in self.strainCons):
                strainSeqs = self.strainCons[strain]
                windowtitle = "Consensus sequences of " + strain
                for gene in strainSeqs:
                    seq = strainSeqs[gene]
                    outseq += seq.format("fasta")
            else:
                return
        else:
            gene = self.geneIds[colnum - 1]
            windowtitle = "Consensus sequence of " + strain + " " + gene
            if(strain in self.strainCons):
                strainSeqs = self.strainCons[strain]
                if(gene in strainSeqs):
                    seq = strainSeqs[gene]
                    outseq = seq.format("fasta")
                else:
                    return
            else:
                return

        if(len(outseq) > 0):
            seqview = SeqViewDlg.SeqViewDlg(1,outseq)
            seqview.setWindowTitle(windowtitle)
            seqview.exec_()

    def _UniqueIDs(self, refseqs, mode):
    ##mode 1 for barcode, mode 2 for primer
        ids = []
        for refinfo in refseqs:
            if(mode == 1): ids.append(refinfo.des)
            elif(mode == 2): ids.append(refinfo.id)
        uniqueids = sorted(set(ids))
        return list(uniqueids)
            
    
class HetViewer(QWidget):
    def __init__(self, outfolder = "", hetSeqs=None, parent=None):
        super(HetViewer, self).__init__(parent)
        self.viewer = QTableWidget()
        self.button = QPushButton('Export table')
        self.hetSeqs = hetSeqs
        buttonlayout = QHBoxLayout()
        buttonlayout.addStretch()
        buttonlayout.addWidget(self.button)
        layout = QVBoxLayout()
        layout.addWidget(self.viewer)
        layout.addLayout(buttonlayout)
        self.outfolder = outfolder
        self.setLayout(layout)
        self.button.clicked.connect(self.exportData)
        self.viewer.doubleClicked.connect(self.showSeqs)
    
    def showData(self, data=None):
        self.data = data
        header = data[0]
        self.colname = header
        self.rowname = []
        rownum = len(data) - 1
        colnum = len(header)
        self.viewer.setRowCount(rownum)
        self.viewer.setColumnCount(colnum)
        self.viewer.setShowGrid(True)
        for i in range(colnum):
            self.viewer.setHorizontalHeaderItem(i, QTableWidgetItem(header[i]))
        for i in range(rownum):
            for j in range(colnum):
                datainfo = str(data[i+1][j])
                if(j == 0):
                    self.rowname.append(datainfo)
                item = QTableWidgetItem(datainfo)
                if(datainfo == 'Yes'):
                    item.setBackground(QColor(249,148,46))
                    item.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEnabled)
                self.viewer.setItem(i,j,QTableWidgetItem(item))
        self.viewer.resizeColumnsToContents()
        self.viewer.horizontalHeader().setStretchLastSection(True)
    
    def showSeqs(self, index):
        rownum = index.row()
        colnum = index.column()
        strain = self.rowname[rownum]
        gene = self.colname[colnum]
        windowtitle = "Alleles of "+ strain + " " + gene
        try:
            seq1 = self.hetSeqs[strain][gene]['seq1']
            seq2 = self.hetSeqs[strain][gene]['seq2']
            outseq = seq1.format('fasta')
            outseq += seq2.format('fasta')
            seqview = SeqViewDlg.SeqViewDlg(1,outseq)
            seqview.setWindowTitle(windowtitle)
            seqview.exec_()
        except:
            return
    
    def exportData(self):
        import csv
        filename = str(QFileDialog.getSaveFileName(
            self, "Export data", self.outfolder, "CSV file (*.csv)"))
        if(filename):
            fileN = filename.split('.')
            if(fileN[-1] != "csv"):
                filename += ".csv"
            try:
                fh_outfile = open(filename,'w')
                csv_writer = csv.writer(fh_outfile, dialect='excel')
                for i in range(len(self.data)):
                    csv_writer.writerow(self.data[i])
                fh_outfile.close()
                QMessageBox.information(self, "Data Exported", "Data exported!")
            except IOError as e:
                QMessageBox.warning(self, "Error:",
                                    "Error: %s" %e.strerror, QMessageBox.Ok|QMessageBox.Default)
                
class TableViewer(QWidget):
    def __init__(self, outfolder = "", alignedSeqs=None, locusLengthRange=None,parent=None):
        super(TableViewer, self).__init__(parent)
        self.viewer = QTableWidget()
        self.button = QPushButton('Export table')
        self.alignedSeqs = alignedSeqs
        self.locusLengthRange = locusLengthRange
        buttonlayout = QHBoxLayout()
        buttonlayout.addStretch()
        buttonlayout.addWidget(self.button)
        layout = QVBoxLayout()
        layout.addWidget(self.viewer)
        layout.addLayout(buttonlayout)
        self.outfolder = outfolder
        self.setLayout(layout)
        self.button.clicked.connect(self.exportData)
        self.viewer.doubleClicked.connect(self.showseqs)
    
    def showData(self, data=None, minReadNum=0):
        #import numpy as np
        #self.data = np.random.rand(96,8).tolist()
        #locusName = ["LAC1","CAP59","SOD1","PLB1","TEF1","URA5","IGS1","GPD1"]
        self.data = data
        header = data[0]
        self.colname = header
        self.rowname = []
        rownum = len(data) - 1
        colnum = len(header)
        self.viewer.setRowCount(rownum)
        self.viewer.setColumnCount(colnum)
        self.viewer.setShowGrid(True)
        for i in range(colnum):
            if(header[i] == "unmapped"):
                header[i] = "UnMap"
            self.viewer.setHorizontalHeaderItem(i, QTableWidgetItem(header[i]))
        for i in range(rownum):
            for j in range(colnum):
                datainfo = str(data[i+1][j])
                if(j == 0):
                    self.rowname.append(datainfo)
                item = QTableWidgetItem(datainfo)
                item.setFlags(Qt.ItemIsSelectable|Qt.ItemIsEnabled)
                if(i < rownum-1):
                    if(j > 0 and j< (colnum - 2)):
                        nums = datainfo.split('(')
                        if(int(nums[0]) < minReadNum):
                            item.setBackground(QColor(170,170,170))
                self.viewer.setItem(i,j,QTableWidgetItem(item))
        self.viewer.resizeColumnsToContents()
        self.viewer.horizontalHeader().setStretchLastSection(True)
        
    def exportData(self):
        import csv
        filename = str(QFileDialog.getSaveFileName(
            self, "Export data", self.outfolder, "CSV file (*.csv)"))
        if(filename):
            fileN = filename.split('.')
            if(fileN[-1] != "csv"):
                filename += ".csv"
            try:
                fh_outfile = open(filename,'w')
                csv_writer = csv.writer(fh_outfile, dialect='excel')
                for i in range(len(self.data)):
                    csv_writer.writerow(self.data[i])
                fh_outfile.close()
                QMessageBox.information(self, "Data Exported", "Data exported!")
            except IOError as e:
                QMessageBox.warning(self, "Error:",
                                    "Error: %s" %e.strerror, QMessageBox.Ok|QMessageBox.Default)

    def showseqs(self, index):
        rownum = index.row()
        colnum = index.column()
        strain = self.rowname[rownum]
        gene = self.colname[colnum]
        windowtitle = "Reads of "+ strain + " " + gene
        lengthRange = self.locusLengthRange[gene]
        (isget, trimmedseqs, untrimmedSeqs) = self.__extractSeqs(strain,gene,lengthRange)
        if(not isget):
            return
        seqview = SeqViewDlg.SeqViewDlg(0, trimmedseqs, untrimmedSeqs)
        seqview.setWindowTitle(windowtitle)
        seqview.exec_()
    
    def __extractSeqs(self, strain, gene, lengthRange):
        #def ExtractSeqs(alignedSeqs,strain,gene,outformat,outmode):
        #stderr.write ('\nExtracting sequences...\n')
        selSeqs = []
        trimmedSeqs = ""
        untrimmedSeqs = ""
        for seq in self.alignedSeqs:
            if(seq.strain == strain and seq.gene == gene):
               selSeqs.append(seq)
        
        if(len(selSeqs) == 0):
            return(False,None,None)

        for seq in selSeqs:
            trimmedseq = seq.TrimPrimer()
            untrimmedseq = seq.seq
            untrimmedSeqs += untrimmedseq.format('fasta')
            if(len(trimmedseq) >= lengthRange['s1'] and len(trimmedseq) <= lengthRange['s2']):
                trimmedSeqs += trimmedseq.format('fasta')
        return (True, trimmedSeqs, untrimmedSeqs)

class Form(QDialog):
    def __init__(self, parent=None):
        super(Form, self).__init__(parent)
        
        #self.viewer = ImageViewer("")
        #self.button = QPushButton("draw")
        self.viewer = ConsensusViewer()
        layout = QVBoxLayout()
        layout.addWidget(self.viewer)
        self.setLayout(layout)

if(__name__ == '__main__'):
    app = QApplication(sys.argv)
    form = Form()
    form.show()
    #data = data*100
    #data = data.tolist()
    #info['data'] = data
    #info['name'] = locusName
    #form.viewer.plotLengthRange(info)
    form.viewer.showData()
    app.exec_()

    