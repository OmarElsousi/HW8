#region imports
import sys
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from Rankine_GUI import Ui_Form
from Rankine_Classes_MVC import rankineController
from UnitConversions import UnitConverter as UC
#these imports are necessary for drawing a matplot lib graph on my GUI
#no simple widget for this exists in QT Designer, so I have to add the widget in code.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
#endregion

#region class definitions
class MainWindow(qtw.QWidget, Ui_Form):
    def __init__(self):
        """
        MainWindow constructor
        """
        super().__init__()
        self.setupUi(self)
        self.AssignSlots()
        self.MakeCanvas()

        # Create lists of input and display widgets
        self.input_widgets = [
            self.rb_SI, self.le_PHigh, self.le_PLow,
            self.le_TurbineInletCondition, self.rdo_Quality,
            self.le_TurbineEff, self.cmb_XAxis, self.cmb_YAxis,
            self.chk_logX, self.chk_logY
        ]
        self.display_widgets = [
            self.lbl_PHigh, self.lbl_PLow, self.lbl_SatPropLow, self.lbl_SatPropHigh,
            self.lbl_TurbineInletCondition, self.lbl_H1, self.lbl_H1Units,
            self.lbl_H2, self.lbl_H2Units, self.lbl_H3, self.lbl_H3Units,
            self.lbl_H4, self.lbl_H4Units, self.lbl_TurbineWork, self.lbl_TurbineWorkUnits,
            self.lbl_PumpWork, self.lbl_PumpWorkUnits, self.lbl_HeatAdded,
            self.lbl_HeatAddedUnits, self.lbl_ThermalEfficiency,
            self.canvas, self.figure, self.ax
        ]
        # Instantiate the controller
        self.RC = rankineController(self.input_widgets, self.display_widgets)

        self.setNewPHigh()
        self.setNewPLow()

        self.Calculate()  # initial calculation

        # a place to store coordinates from last position on graph
        self.oldXData = 0.0
        self.oldYData = 0.0
        self.show()

    def AssignSlots(self):
        self.btn_Calculate.clicked.connect(self.Calculate)
        self.rdo_Quality.clicked.connect(self.SelectQualityOrTHigh)
        self.rdo_THigh.clicked.connect(self.SelectQualityOrTHigh)
        # Note that we connect editingFinished to setNewPHigh / setNewPLow
        self.le_PHigh.editingFinished.connect(self.setNewPHigh)
        self.le_PLow.editingFinished.connect(self.setNewPLow)
        # SI/English radio buttons
        self.rb_SI.clicked.connect(self.SetUnits)
        self.rb_English.clicked.connect(self.SetUnits)
        # Plot controls
        self.cmb_XAxis.currentIndexChanged.connect(self.SetPlotVariables)
        self.cmb_YAxis.currentIndexChanged.connect(self.SetPlotVariables)
        self.chk_logX.toggled.connect(self.SetPlotVariables)
        self.chk_logY.toggled.connect(self.SetPlotVariables)

    def MakeCanvas(self):
        self.figure = Figure(figsize=(1,1), tight_layout=True, frameon=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot()
        self.Layout_Plot.addWidget(NavigationToolbar2QT(self.canvas, self))
        self.Layout_Plot.addWidget(self.canvas)
        self.canvas.mpl_connect("motion_notify_event", self.mouseMoveEvent_Canvas)

    def mouseMoveEvent_Canvas(self, event):
        self.oldXData = event.xdata if event.xdata is not None else self.oldXData
        self.oldYData = event.ydata if event.ydata is not None else self.oldYData
        sUnit = 'kJ/(kg*K)' if self.rb_SI.isChecked() else 'BTU/(lb*R)'
        TUnit = 'C' if self.rb_SI.isChecked() else 'F'
        self.setWindowTitle(f's: {self.oldXData:0.2f} {sUnit}, '
                            f'T: {self.oldYData:0.2f} {TUnit}')

    def Calculate(self):
        self.RC.updateModel()

    def SelectQualityOrTHigh(self):
        self.RC.selectQualityOrTHigh()

    def SetPlotVariables(self):
        self.RC.updatePlot()

    # -------------------------------
    # CHANGES FOR OBJECTIVES 2 & 3:
    # -------------------------------
    def SetUnits(self):
        """
        Called when the user clicks SI or English radio button.
        We immediately update all labels/values/plots in the GUI
        via the controller's updateUnits().
        """
        self.RC.updateUnits()

    def setNewPHigh(self):
        """
        Called when user finishes editing P High (editingFinished signal).
        Let the controller update the saturations and any T High if needed,
        and then immediately refresh model results.
        """
        self.RC.setNewPHigh()
        self.RC.updateModel()  # so user sees new saturation properties right away

    def setNewPLow(self):
        """
        Similar to setNewPHigh, but for P Low.
        """
        self.RC.setNewPLow()
        self.RC.updateModel()

#endregion

if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    mw = MainWindow()
    mw.setWindowTitle('Rankine calculator')
    sys.exit(app.exec())
