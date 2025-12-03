"""
The QT interface has been taken from a public domain software from
Eli Bendersky (eliben@gmail.com), updated by Ondrej Holesovsky.
"""
import sys, os
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import numpy as np
from scipy import signal as sp_signal


class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Effects of Arrivals and Services on the performances')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        self.generateRnd = True

        self.on_draw()
        
        
    def runModel(self, arrv, arrCVv, arrCorv, srvv, srvCVv, srvCorv):
 #       print('CV: ', arrCVv, srvCVv)
        rSamples = 25000
        corrEps = 0.0001
        
        if self.generateRnd:
            self.baserand = np.random.randn(rSamples, 2)
            self.generateRnd = False
        
        asqrt1pcv2 = np.sqrt(1.0 + arrCVv * arrCVv)
        ssqrt1pcv2 = np.sqrt(1.0 + srvCVv * srvCVv)
        asqrtlog1pcv2 = np.sqrt(np.log(1.0 + arrCVv * arrCVv))
        ssqrtlog1pcv2 = np.sqrt(np.log(1.0 + srvCVv * srvCVv))
        
        randv = self.baserand.copy()
        ac = arrCorv
        if ac == 1.0:
            ac = ac - corrEps
        elif ac == -1.0:
            ac = ac + corrEps
        AaFilter = np.array([1, -ac])
        AbFilter = np.sqrt(1-ac*ac)
        randv[0,0] = randv[0,0] / AbFilter
        randv[:,0] = sp_signal.lfilter(AbFilter, AaFilter, randv[:,0]) 

        sc = srvCorv
        if sc == 1.0:
            sc = sc - corrEps
        elif sc == -1.0:
            sc = sc + corrEps
        SaFilter = np.array([1, -sc])
        SbFilter = np.sqrt(1-sc*sc)
        randv[0,1] = randv[0,1] / SbFilter
        randv[:,1] = sp_signal.lfilter(SbFilter, SaFilter, randv[:,1]) 

        
#        print('a1osqrt1pcv2: ', a1osqrt1pcv2, ', s1osqrt1pcv2', s1osqrt1pcv2)
        samp = np.exp(randv * np.array([asqrtlog1pcv2, ssqrtlog1pcv2])) * np.array([arrv / asqrt1pcv2, srvv / ssqrt1pcv2])

        act = np.zeros((rSamples, 3))
        act[:,0] = samp[:,0].cumsum()
        act[0,1] = act[0,0] + samp[0, 1]
        act[0,2] = samp[0, 1]
        for i in range(1, rSamples):
            act[i, 1] = max(act[i, 0], act[i-1, 1]) + samp[i, 1]
            act[i, 2] = act[i, 1] - act[i, 0]

        avgiat = np.mean(samp[:,0])
        stdiat = np.std(samp[:,0])
        arate = 1.0 / avgiat
        avgst  = np.mean(samp[:,1])
        stdst = np.std(samp[:,1])
        U = min(arate * avgst, 1.0)
        X = min(arate, 1.0 / avgst)
        R = np.mean(act[:, 2])
        N = X * R
        
        perfidx = {'avgiat': avgiat, 'stdiat': stdiat, 'arate': arate, 'avgst': avgst, 'stdst': stdst, 'X': X, 'U': U, 'R': R, 'N': N}
        
        return {'samp': samp, 'perfidx': perfidx, 'act': act}
   
    def on_about(self):
        msg = """ Arrival and Services:
        
         * See the effect of average interarrival rates and service times
         * See the effect of the variance to the performances
         * See the effect of correlation
        """
        QMessageBox.about(self, "About this tool", msg.strip())
        
    def on_draw(self):
        """ Redraws the figure
        """
        
        # clear the axes and redraw the plot anew
        #
        self.axes.clear()        
        
        arrv    = self.Arrivals.getVal()
        arrCVv  = self.ArrivalsCV.getVal()
        arrCorv = self.ArrivalsCor.getVal()
        srvv    = self.Services.getVal()
        srvCVv  = self.ServicesCV.getVal()
        srvCorv = self.ServicesCor.getVal()
        
        res = self.runModel(arrv, arrCVv, arrCorv, srvv, srvCVv, srvCorv)
        #print(np.mean(res,0))
        avgiat = res['perfidx']['avgiat']
        stdiat = res['perfidx']['stdiat']
        arate  = res['perfidx']['arate']
        avgst  = res['perfidx']['avgst']
        stdst  = res['perfidx']['stdst']
        X      = res['perfidx']['X']
        U      = res['perfidx']['U']
        R      = res['perfidx']['R']
        N      = res['perfidx']['N']

        self.TextOutput.setHtml(f"<H1>Arrivals</H1><BR><B>Average Inter-arriavl Time:</B> <I>{avgiat}</I><BR><B>Standard Deviation</B> <I>{stdiat}</I><BR><B>Arriavl Rate:</B> <I>{arate}</I><BR> <H1>Services</H1><BR><B>Average Service Time:</B> <I>{avgst}</I><BR><B>Standard Deviation</B> <I>{stdst}</I><BR> <H1>Performance Indices</H1><BR><B>Utilization:</B> <I>{U}</I><BR><B>Throughput:</B> <I>{X}</I><BR><B>Average Response Time:</B> <I>{R}</I><BR><B>Average Number of Jobs:</B> <I>{N}</I><BR>")
        
        vSamp = 500
        self.axes.plot(res['samp'][0:vSamp,0], range(0,vSamp), "+", res['samp'][0:vSamp,1], range(0,vSamp), "+")
        self.axes.set_xlim([0, 25])
        self.axes.legend(['Arrivals', 'Services'], loc='upper left')
        self.canvas.draw()
        
    class sliderAndText(object):
        slider = None
        
        def __init__(self, rmin, rmax, rval, scale, callback):
            self.scale = scale
            frval = rval / scale
            self.label = QLabel(f'{frval}')
            self.slider = QSlider(Qt.Horizontal)
            self.slider.setRange(rmin, rmax)
            self.slider.setValue(rval)
            self.slider.setTracking(True)
            self.slider.setTickPosition(QSlider.TicksBothSides)
            self.slider.valueChanged.connect(callback)
        
        def getVal(self):
            val = self.slider.value()
            frval = val / self.scale
            self.label.setText(f'{frval}')
            return frval
    
    def create_main_frame(self):
        self.main_frame = QWidget()
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((10.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(111)
        
        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        
        # Other GUI controls
        # 
        L00 = QLabel('')
        L01 = QLabel('Average')
        L02 = QLabel('Coef. of Var.')
        L03 = QLabel('Correlation')
        L10 = QLabel('Inter-Arrivals')
        L20 = QLabel('')
        L30 = QLabel('Services')
        L40 = QLabel('')
        
        self.Arrivals = self.sliderAndText(1, 100, 16, 10, self.on_draw)
        self.Services = self.sliderAndText(1, 100,  5, 10, self.on_draw)
        self.ArrivalsCV = self.sliderAndText(0, 200, 10, 20, self.on_draw)
        self.ServicesCV = self.sliderAndText(0, 200, 10, 20, self.on_draw)
        self.ArrivalsCor = self.sliderAndText(-100, 100, 0, 100, self.on_draw)
        self.ServicesCor = self.sliderAndText(-100, 100, 0, 100, self.on_draw)
        
        grid = QGridLayout()
        
        i = 0
        cpr = 4
        for w in [L00, L01, L02, L03,
                  L10, self.Arrivals.slider, self.ArrivalsCV.slider, self.ArrivalsCor.slider,
                  L20, self.Arrivals.label,  self.ArrivalsCV.label,  self.ArrivalsCor.label,
                  L30, self.Services.slider, self.ServicesCV.slider, self.ServicesCor.slider,
                  L20, self.Services.label,  self.ServicesCV.label,  self.ServicesCor.label
                 ]:
            grid.addWidget(w, i // cpr, i % cpr)
            grid.setAlignment(w, Qt.AlignVCenter)
            i = i + 1

        self.TextOutput = QTextEdit()
        self.TextOutput.setReadOnly(True)

        vbox = QVBoxLayout()
        chbox = QHBoxLayout()
        chbox.addWidget(self.canvas)
        chbox.addWidget(self.TextOutput)
        vbox.addLayout(chbox)
        vbox.addWidget(self.mpl_toolbar)
        vbox.addLayout(grid)
        
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QLabel("Running")
        self.statusBar().addWidget(self.status_text, 1)
        
    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("&File")
        
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (quit_action,))
        
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", 
            shortcut='F1', slot=self.on_about, 
            tip='About the demo')
        
        self.add_actions(self.help_menu, (about_action,))

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False):
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


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
