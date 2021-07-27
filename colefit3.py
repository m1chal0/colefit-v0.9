import sys
from PyQt5.QtWidgets import QApplication, QWidget, QHBoxLayout, QVBoxLayout, QLabel, QPushButton, QMessageBox, QFileDialog, QGridLayout
from PyQt5 import QtCore, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize.minpack import curve_fit
from mpl_toolkits import mplot3d
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Cursor
import matplotlib.patches as patches
import os
import sys

#pyinstaller --onefile --clean --windowed --icon=icon.ico  colefit3.py

class LoadApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle('ColeFIT 1.0')
        self.setWindowIcon(QtGui.QIcon('res/icon.ico'))
        self.setMinimumSize(600, 480)
        self.setMaximumSize(600, 480)
        self.setWindowFlag(QtCore.Qt.FramelessWindowHint)

        loadlay = QVBoxLayout()
        self.labelT = QLabel("Tato aplikace TEST")
        loadlay.addWidget(self.labelT)
        self.setLayout(loadlay)



class MyApp(QWidget):
    file_path = ""
    tpreal = []
    tpimag = []
    tf = []
    tt = []
    teplota_sel = 0
    pi=np.pi
    sin=np.sin
    cos=np.cos
    arctan = np.arctan
    REALNAE = []
    IMAGINARNIE = []
    FREKVENCE = []
    limitfromplot = 0
    cparams = ""

    def __init__(self):
        global file_path
        super().__init__()
        self.setWindowTitle('ColeFIT 0.9')
        self.setWindowIcon(QtGui.QIcon('res/icon.ico'))
        self.setMinimumSize(800, 800)
        self.setMaximumSize(800, 800)
        self.w = None
        layout = QVBoxLayout()

        self.button = QPushButton("1. Vybrat data", self)
        self.button.clicked.connect(self.openFileNamesDialog)

        self.buttonP = QPushButton("4. Uložit parametry", self)
        self.buttonP.clicked.connect(self.SaveWindowBtn)
        self.buttonP.setEnabled(False) 

        self.button2 = QPushButton("2. Zpracovat data", self)
        self.button2.clicked.connect(self.workWithData)
        self.button2.setEnabled(False)  

        self.button3 = QPushButton("3. Provézt Fit", self)
        self.button3.clicked.connect(self.fitdata)
        self.button3.setEnabled(False) 

        self.buttonH = QPushButton("Návod", self)
        self.buttonH.clicked.connect(self.openHelp)

        self.buttonR = QPushButton("Restart", self)
        self.buttonR.clicked.connect(self.restartApp)

        self.canvas = FigureCanvas(plt.Figure(figsize=(10, 10), dpi=60))
        self.canvas.mpl_connect('button_press_event', self.data_cursor_sub)

        self.canvas2 = FigureCanvas(plt.Figure(figsize=(10, 10), dpi=60))

        self.canvas3 = FigureCanvas(plt.Figure(figsize=(10, 10), dpi=60))
        

        self.toolbar1 = NavigationToolbar(self.canvas, self)
        self.toolbar2 = NavigationToolbar(self.canvas2, self)
        self.toolbar3 = NavigationToolbar(self.canvas3, self)

        self.toolbar1.setMinimumWidth(self.canvas.width())
        self.toolbar2.setMinimumWidth(self.canvas2.width())
        self.toolbar3.setMinimumWidth(self.canvas3.width())

        self.coleparam = QLabel()
        self.coleparam.setText("")

        layout = QGridLayout()

        layout.addWidget(self.button,0,0)
        layout.addWidget(self.button2,0,1)
        layout.addWidget(self.button3,0,2)

        layout.addWidget(self.toolbar1,1,0,1,2)
        layout.addWidget(self.canvas,2,0,1,2)

        layout.addWidget(self.toolbar2,1,2,1,2)
        layout.addWidget(self.canvas2,2,2,1,2)

        layout.addWidget(self.toolbar3,3,0,1,4)
        layout.addWidget(self.canvas3,4,0,1,4)

        layout.addWidget(self.buttonP,0,3)
        layout.addWidget(self.buttonH,5,0,1,1)
        layout.addWidget(self.buttonR,5,3,1,1)

        self.insert_ax()

        self.setLayout(layout)
        
    def restartApp(self):
        os.execl(sys.executable, os.path.abspath(__file__), *sys.argv)

    def openHelp(self):
        self.w = HelpW()
        self.w.show()

    def SaveWindowBtn(self):
        global cparams
        global file_path
        with open(file_path + "_coleparam.txt", "w", encoding="utf-8") as f:
            f.write(cparams)
            f.close()
            QMessageBox.about(self, "Zpráva", "Uloženo do "+file_path + "_coleparam.txt")

    def find_nearest(self, array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def data_cursor_sub(self, event):
        global limitfromplot
        global FREKVENCE
        global dlimit
        global hlimit
        limitfromplot = 0
        
        while True:
            try: 
                dlimit
            except NameError:
                limitfromplot = self.find_nearest(FREKVENCE, event.xdata)
                rect = patches.Rectangle((self.ax.get_xlim()[0], self.ax.get_ylim()[0]), limitfromplot, self.ax.get_ylim()[1], alpha=0.2, facecolor='grey')
                self.ax.add_patch(rect)
                dlimit = limitfromplot
                break
            else:
                try:
                    hlimit
                except NameError:
                    limitfromplot = self.find_nearest(FREKVENCE, event.xdata)
                    rect = patches.Rectangle((limitfromplot, self.ax.get_ylim()[0]), (self.ax.get_xlim()[1] - limitfromplot), self.ax.get_ylim()[1], alpha=0.2, facecolor='grey')
                    self.ax.add_patch(rect)
                    hlimit = limitfromplot
                    break
                else:
                    self.ax.patches = []
                    del dlimit
                    del hlimit
                    continue

        

    def openFileNamesDialog(self):
        global file_path
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getOpenFileName(self,"Vybrat soubor s daty", "","Soubor z Novocontrol ALFA *txt", options=options)
        self.button2.setEnabled(True)
    
    def insert_ax(self):
        font = {
            'weight': 'normal',
            'size': 14
        }
        plt.rc('font', **font)
        self.ax = self.canvas.figure.subplots()
        self.ax2 = self.canvas2.figure.subplots()
        self.ax3 = self.canvas3.figure.subplots()          

    def update_chart(self, FREKVENCE, IMAGINARNIE, REALNAE):
        color = "b"
        color_alpha = 1
        self.ax.clear()
        self.ax2.clear()
        self.ax3.clear()

        self.ax.set_xlabel('Frequency [Hz]', fontdict=dict(weight='bold'))
        self.ax.set_ylabel('Imag. permittivity', fontdict=dict(weight='bold'))
        self.canvas.figure.tight_layout(pad=2)
        self.ax.grid(True)

        self.ax2.set_xlabel('Frequency [Hz]', fontdict=dict(weight='bold'))
        self.ax2.set_ylabel('Real. permittivity', fontdict=dict(weight='bold'))
        self.canvas2.figure.tight_layout(pad=2)
        self.ax2.grid(True)

        self.ax3.set_xlabel('Real. permittivity', fontdict=dict(weight='bold'))
        self.ax3.set_ylabel('Imag. permittivity', fontdict=dict(weight='bold'))
        self.canvas3.figure.tight_layout(pad=2)
        self.ax3.grid(True)

        self.ax.semilogx(FREKVENCE, IMAGINARNIE, color=color, alpha=color_alpha, linestyle='solid')
        cursor = Cursor(self.ax, horizOn = False, vertOn=True, color='red', linewidth=1)
        self.canvas.draw()

        self.ax2.semilogx(FREKVENCE, REALNAE, color=color, alpha=color_alpha, linestyle='solid')
        self.canvas2.draw()

        self.ax3.plot(REALNAE, IMAGINARNIE, color=color, alpha=color_alpha, linestyle='solid')
        #self.canvas3.figure.ginput(2)

        self.canvas3.draw()
        self.button3.setEnabled(True)
        

    def update_chart_fit(self, FREKVENCE, IMAGINARNIE_FIT, REALNAE_FIT, REALCALC, IMAGCALC, PERMSTAT, PERMNEK, TAU, ALFA, BETA):
        color = "b"
        color_alpha = 1
        global cparams
        self.ax.patches = []
        self.ax.semilogx(FREKVENCE, IMAGINARNIE_FIT, alpha=color_alpha, linestyle='dashed')
        self.canvas.draw()
       # self.canvas.figure.ginput(4)

        self.ax2.semilogx(FREKVENCE, REALNAE_FIT, alpha=color_alpha, linestyle='dashed')
        self.canvas2.draw()

        line, = self.ax3.plot(REALCALC, IMAGCALC, alpha=color_alpha, linestyle='dashed')
        scient_tau = "{:.2e}".format(TAU)
        maximum = np.argmax(IMAGCALC)
        self.ax3.text(self.ax3.get_xlim()[0]+0.01*self.ax3.get_xlim()[0], IMAGCALC[maximum], "ε_stat: "+ str(np.round(PERMSTAT,3)) + " ε_nek: "+ str(np.round(PERMNEK,3)) + " λ: "+ str(scient_tau) + " α: "+ str(np.round(ALFA,3)) + " β: "+ str(np.round(BETA,3)), color = line.get_color(), fontweight = "black", fontsize = "medium")
        self.canvas3.draw()

        try:
            cparams
        except NameError:
            cparams = ""

        params = "ε_stat "+ str(np.round(PERMSTAT,3)) + " ε_nek "+ str(np.round(PERMNEK,3)) + " λ "+ str(scient_tau) + " α "+ str(np.round(ALFA,3)) + " β "+ str(np.round(BETA,3))+"\n"
        cparams = cparams + params
        self.buttonP.setEnabled(True) 

    def workWithData(self):
        global file_path, tpreal, tpimag, teplota_sel, tf, tt, pi, sin, cos, arctan, REALNAE, IMAGINARNIE, FREKVENCE
        self.button2.setEnabled(False)
        self.button.setEnabled(False)
        teplota_sel = -50.0
        tpreal = []
        tpimag = []
        tf = []
        tt = []
        pi=np.pi
        sin=np.sin
        cos=np.cos
        arctan = np.arctan
        REALNAE = []
        IMAGINARNIE = []
        FREKVENCE = []

        print(teplota_sel)
        if file_path:
            f=open(file_path,"r")
            lines=f.readlines()
            f.close()
            x=0
            route = file_path.split("/")
            name_of_file = route[len(route) - 1]
            print(name_of_file)
            for x in range(len(lines)):
                line = lines[x].split()
                if x >=3 :
                    if float(line[1]) == teplota_sel:
                        tf.append(line[0])
                        tt.append(line[1])
                        tpreal.append(line[2])
                        tpimag.append(line[3])
            tpreal = list(map(float, tpreal))
            tpimag = list(map(float, tpimag))
            tf = list(map(float, tf))
            tt = list(map(float, tt))
            n = int(len(tpreal)/1)
            print(n)
            def divide_chunks(l, n): 
                for i in range(0, len(l), n):  
                    yield l[i:i + n] 
            tpreal = list(divide_chunks(tpreal, n))
            tpimag = list(divide_chunks(tpimag, n))
            tf = list(divide_chunks(tf, n))
            tt = list(divide_chunks(tt, n))

            REALNAE = np.array(tpreal[0])
            IMAGINARNIE = np.array(tpimag[0])
            FREKVENCE = np.array(tf[0])
            TEPLOTA = np.array(tt[0])
        self.update_chart(FREKVENCE, IMAGINARNIE, REALNAE)
            

    def fitdata(self):
        global pi, sin, cos, arctan, REALNAE, IMAGINARNIE, FREKVENCE
        global dlimit, hlimit
        windex = [0,0]

        try:
            hlimit
        except NameError:
            windex = [len(FREKVENCE), 0]
        else:
            windex[0] = np.where(FREKVENCE == dlimit)[0]
            windex[1] = np.where(FREKVENCE == hlimit)[0]
            print(windex[0])
        working_FREKVENCE = []
        working_REALNAE = []
        working_IMAGINARNIE = []

        for i in range(len(FREKVENCE)):
            if windex[1] < i <= windex[0]:
                working_FREKVENCE.append(FREKVENCE[i])
                working_REALNAE.append(REALNAE[i])
                working_IMAGINARNIE.append(IMAGINARNIE[i])

        def func(frequency, PERMNEK, PERMSTAT, ALFA, BETA, TAU):
            fi = arctan(((2*pi*frequency)**(ALFA))*(TAU**(ALFA))*sin((ALFA*pi)/2)/(1+((2*pi*frequency)**(ALFA))*(TAU**(ALFA))*cos((ALFA*pi)/2)))
            return (PERMNEK-PERMSTAT)*(sin(BETA*fi))/((1+2*((2*pi*frequency)**(ALFA))*(TAU**(ALFA))*cos((ALFA*pi)/2)+((2*pi*frequency)**(2*ALFA))*(TAU**(2*ALFA)))**(BETA/2))
        popt, pcov = curve_fit(func, working_FREKVENCE, working_IMAGINARNIE, method="trf", bounds=(0, [10., 10., 1., 1., 10.]), maxfev=10000000)
        IMAGINARNIE_FIT = func(FREKVENCE, *popt)

        FITTED_PARAMETERS = tuple(popt)
        ALFA = FITTED_PARAMETERS[2]
        BETA = FITTED_PARAMETERS[3]
        TAU = FITTED_PARAMETERS[4]
        #maximum = np.argmax(IMAGINARNIE_FIT)
        #TAUC = 1/(2 * np.pi * FREKVENCE[maximum])
        #print(TAU)
        #print(maximum)
        #print(FREKVENCE[maximum])
        #print(TAUF)

        def func2(frequency,PERMNEK,PERMSTAT):                      
            fi = arctan(((2*pi*frequency)**(ALFA))*(TAU**(ALFA))*sin((ALFA*pi)/2)/(1+((2*pi*frequency)**(ALFA))*(TAU**(ALFA))*cos((ALFA*pi)/2)))
            return PERMNEK+(PERMSTAT-PERMNEK)*cos(BETA*fi)/((1+2*((2*pi*frequency)**(ALFA))*(TAU**(ALFA))*cos((ALFA*pi)/2)+((2*pi*frequency)**(2*ALFA))*(TAU**(2*ALFA)))**(BETA/2))
        popt, pcov = curve_fit(func2, working_FREKVENCE, working_REALNAE, method="trf", bounds=(0, [10., 10.]), maxfev=10000000)
        REALNAE_FIT = func2(FREKVENCE, *popt)

        FITTED_PARAMETERS = tuple(popt)
        PERMNEK = FITTED_PARAMETERS[0]
        PERMSTAT = FITTED_PARAMETERS[1]

        REALCALC = []
        IMAGCALC = []
        flist = np.logspace(-10, 100, num=100000)
        for fr in flist:
            frequency = fr
            fi = np.arctan(((2*pi*frequency)**(ALFA))*(TAU**(ALFA))*np.sin((ALFA*pi)/2)/(1+((2*pi*frequency)**(ALFA))*(TAU**(ALFA))*np.cos((ALFA*pi)/2)))
            REALCALCitem = PERMNEK+(PERMSTAT-PERMNEK)*cos(BETA*fi)/((1+2*(2*pi*frequency)**(ALFA)*TAU**(ALFA)*cos(ALFA*pi/2)+(2*pi*frequency)**(2*ALFA)*TAU**(2*ALFA))**(BETA/2))
            IMAGCALCitem = - (PERMNEK-PERMSTAT)*np.sin(BETA*np.arctan((2*pi*frequency)**(ALFA)*TAU**(ALFA)*np.sin(ALFA*pi/2)/(1+(2*pi*frequency)**(ALFA)*TAU**(ALFA)*np.cos(ALFA*pi/2))))/((1+2*(2*pi*frequency)**(ALFA)*TAU**(ALFA)*np.cos(ALFA*pi/2)+(2*pi*frequency)**(2*ALFA)*TAU**(2*ALFA))**(BETA/2))
            REALCALC.append(REALCALCitem)
            IMAGCALC.append(IMAGCALCitem)
        self.update_chart_fit(FREKVENCE, IMAGINARNIE_FIT, REALNAE_FIT, REALCALC, IMAGCALC, PERMSTAT, PERMNEK, TAU, ALFA, BETA)

class HelpW(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Návod k aplikaci ColeFIT")
        self.setWindowIcon(QtGui.QIcon('res/iconH.ico'))
        layoutH = QVBoxLayout()
        self.labelH = QLabel('ColeFIT ver. 0.9', self)
        self.label = QLabel('Tato aplikace slouží pro analýzu parametrů Cole-Cole diagramu. Konkrétně je použito vyjádření Harviliak-Negami, které slouží pro popis složitějších amorfních struktur.', self)
        self.label2 = QLabel('1. Klikněte na tlačítko "Vybrat Data". Zvolte soubor s daty ze spektroskopu NOVOCONTROL ALPHA.', self)
        self.label3 = QLabel('2. Klikněte na tlačítko "Zpracovat data". Dojde ke zpracování surových dat a zobrazení grafů.', self)
        self.label4 = QLabel('3. Pokud chcete vybrat rozsah frekvencí, pro který má být udělán fit, zvolte jej pomocí kurzoru na prvním grafu. Jinak bude využit celý rozsah frekvencí.', self)
        self.label5 = QLabel('4. Klikněte na tlačítko "Provézt Fit". Dojde k získání parametrů diagramu na základě teoretických rovnic a naměřených dat. Ty lze následně uložit do .txt souboru pomocí tlačítka "Uložit Parametry".', self)
        
        layoutH.addWidget(self.labelH)
        layoutH.addWidget(self.label)
        layoutH.addWidget(self.label2)
        layoutH.addWidget(self.label3)
        layoutH.addWidget(self.label4)
        layoutH.addWidget(self.label5)

        self.setLayout(layoutH)

if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    app.setStyleSheet('''
        QWidget {
            font-size: 16px;
        }
        QLabel {
            font-size: 12px;
        }
    ''')
    loadApp = LoadApp()
    loadApp.show()
    myApp = MyApp()
    myApp.show()
    loadApp.hide()

    try:
        sys.exit(app.exec_())
    except SystemExit:
        print('Closing Window...')