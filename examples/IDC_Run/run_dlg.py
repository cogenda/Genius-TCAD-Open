from run import run

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import os

class run_dlg(QDialog):

    def __init__(self, parent=None):
        super(run_dlg, self).__init__(parent)

        fvalidate = QDoubleValidator()
        fvalidate.setNotation(QDoubleValidator.StandardNotation|
                              QDoubleValidator.ScientificNotation)
        fvalidate.setBottom(0)

        self.setWindowTitle(self.tr("Run"))

        vLayout = QVBoxLayout()
        mLayout = QGridLayout()

        mLayout.addWidget(QLabel(self.tr("Source file")), 0, 0)
        self.leSrcFile = QLineEdit()
        self.btSrcFile = QPushButton("...")
        mLayout.addWidget(self.leSrcFile, 0, 1)
        mLayout.addWidget(self.btSrcFile, 0, 2)

        mLayout.addWidget(QLabel(self.tr("Track file")), 1, 0)
        self.leTrkFile = QLineEdit()
        self.btTrkFile = QPushButton("...")
        mLayout.addWidget(self.leTrkFile, 1, 1)
        mLayout.addWidget(self.btTrkFile, 1, 2)

        mLayout.addWidget(QLabel(self.tr("Particle file")), 2, 0)
        self.lePtkFile = QLineEdit()
        self.btPtkFile = QPushButton("...")
        mLayout.addWidget(self.lePtkFile, 2, 1)
        mLayout.addWidget(self.btPtkFile, 2, 2)

        mLayout.addWidget(QLabel(self.tr("Particle weight")), 3, 0)
        self.lePtkWeight = QLineEdit("1.0")
        self.lePtkWeight.setValidator(fvalidate)
        self.lePtkWeight.setMaximumWidth(100)
        mLayout.addWidget(self.lePtkWeight, 3, 1)

        mLayout.addWidget(QLabel(self.tr("Grid resolution(mm)")), 4, 0)
        self.leGridRslt = QLineEdit("0.05")
        self.leGridRslt.setValidator(fvalidate)
        self.leGridRslt.setMaximumWidth(100)
        mLayout.addWidget(self.leGridRslt, 4, 1)

        vLayout.addLayout(mLayout)

        hLayout = QHBoxLayout()
        self.btRun = QPushButton(self.tr("Run"))
        self.btRun.setMaximumWidth(120)
        hLayout.addWidget(self.btRun)
        self.btCancle = QPushButton(self.tr("Cancle"))
        self.btCancle.setMaximumWidth(120)
        hLayout.addWidget(self.btCancle)
        vLayout.addLayout(hLayout)

        self.setLayout(vLayout)

        self.resize(600, 200)

        self.connect(self.btSrcFile, SIGNAL('clicked()'), self.onBtSrcFileClicked)
        self.connect(self.btTrkFile, SIGNAL('clicked()'), self.onBtTrkFileClicked)
        self.connect(self.btPtkFile, SIGNAL('clicked()'), self.onBtPtkFileClicked)
        self.connect(self.btRun, SIGNAL('clicked()'), self.onBtRunClicked)
        self.connect(self.btCancle, SIGNAL('clicked()'), self.reject)

    def onBtSrcFileClicked(self):
        homepath = QDir.homePath()
        path = str( QFileDialog.getOpenFileName(self, self.tr("Open File"),\
                                    homepath, self.tr("Model file(*.unv *.tif3d)")) )
        if len(path)==0: return
        self.leSrcFile.setText(path)

    def onBtTrkFileClicked(self):
        homepath = QDir.homePath()
        path = str( QFileDialog.getOpenFileName(self, self.tr("Open File"),\
                                    homepath, self.tr("Track file(*.*)")) )
        if len(path)==0: return
        self.leTrkFile.setText(path)

    def onBtPtkFileClicked(self):
        homepath = QDir.homePath()
        path = str( QFileDialog.getOpenFileName(self, self.tr("Open File"),\
                                    homepath, self.tr("Particle file(*.*)")) )
        if len(path)==0: return
        self.lePtkFile.setText(path)

    def onBtRunClicked(self):
        srcFile = str(self.leSrcFile.text())
        if not os.path.isfile(srcFile):
            msg = self.tr("Source File '%s' does not exist"%srcFile)
            QMessageBox.warning(None, self.tr("File Error!"), msg)
            return

        trkFile = str(self.leTrkFile.text())
        if not os.path.isfile(trkFile):
            msg = self.tr("Track File '%s' does not exist"%trkFile)
            QMessageBox.warning(None, self.tr("File Error!"), msg)
            return

        ptkFile = str(self.lePtkFile.text())
        if not os.path.isfile(ptkFile):
            msg = self.tr("Particle File '%s' does not exist"%ptkFile)
            QMessageBox.warning(None, self.tr("File Error!"), msg)
            return

        ptkWeight = str(self.lePtkWeight.text())
        if len(ptkWeight)==0:
            msg = self.tr("Particle weight is not set")
            return
        ptkWeight = float(ptkWeight)

        gridRslt = str(self.leGridRslt.text())
        if len(gridRslt)==0:
            msg = self.tr("Grid resolution is not set")
            return
        gridRslt = float(gridRslt)

        run(srcFile, trkFile, ptkFile, ptkWeight, gridRslt)

    def reject(self):
        super(run_dlg, self).reject()



if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    dlg = run_dlg()

    dlg.show()
    app.exec_()
