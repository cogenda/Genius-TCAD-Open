from gdmlgen import gdmlgen

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import os

class gdmlgen_dlg(QDialog):

    def __init__(self, parent=None):
        super(gdmlgen_dlg, self).__init__(parent)

        self.setWindowTitle(self.tr("gdml generator"))

        vLayout = QVBoxLayout()
        hLayout = QHBoxLayout()
        hLayout.addWidget(QLabel(self.tr("3D Model File: ")))
        self.leFPath = QLineEdit()
        hLayout.addWidget(self.leFPath)
        self.btFPath = QPushButton("...")
        self.btFPath.setMaximumWidth(50)
        self.btFPath.setMaximumHeight(25)
        hLayout.addWidget(self.btFPath)
        vLayout.addLayout(hLayout)

        hLayout = QHBoxLayout()
        self.btGdmlGen = QPushButton(self.tr("Generate GDML"))
        self.btGdmlGen.setMaximumWidth(120)
        hLayout.addWidget(self.btGdmlGen)
        self.btCancle = QPushButton(self.tr("Cancle"))
        self.btCancle.setMaximumWidth(120)
        hLayout.addWidget(self.btCancle)
        vLayout.addLayout(hLayout)

        self.setLayout(vLayout)

        self.resize(500, 100)

        self.connect(self.btFPath, SIGNAL('clicked()'), self.onBtFPathClicked)
        self.connect(self.btGdmlGen, SIGNAL('clicked()'), self.onbtGdmlGenClicked)
        self.connect(self.btCancle, SIGNAL('clicked()'), self.reject)

    def onBtFPathClicked(self):
        homepath = QDir.homePath()
        path = str( QFileDialog.getOpenFileName(self, self.tr("Open File"),\
                                                homepath, self.tr("Model file(*.unv *.tif3d)")) )
        if len(path)==0: return

        self.leFPath.setText(path)

    def onbtGdmlGenClicked(self):
        path = str(self.leFPath.text())
        if not os.path.isfile(path):
            msg = self.tr("File '%s' does not exist"%path)
            QMessageBox.warning(None, self.tr("File Error!"), msg)
            return

        gdmlgen(path)

    def reject(self):
        super(gdmlgen_dlg, self).reject()



if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    dlg = gdmlgen_dlg()

    dlg.show()
    app.exec_()