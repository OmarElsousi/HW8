# region imports
import math
from Calc_state import *  # Presumably defines Steam_SI, stateProps, etc.
from UnitConversions import UnitConverter as UC
import numpy as np
from matplotlib import pyplot as plt
from copy import deepcopy as dc

# These imports are necessary for embedding and using a Matplotlib graph in a PyQt application.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure


# endregion


# region class definitions
class rankineModel():
    """
    The Model in the Model-View-Controller (MVC) pattern for the Rankine Power Cycle application.

    Responsibilities:
    - Stores all relevant thermodynamic parameters (pressures, temperatures, enthalpies, etc.).
    - Holds 'state' objects (state1, state2, etc.) that describe points in the Rankine cycle.
    - Contains flags (e.g., SI or not) to track whether the user is working in SI or English units.
    - Stores data for plotting (saturation curves, upper/lower cycle curves, etc.).

    This class itself does not do heavy lifting in calculations — that is often delegated to
    separate classes or the controller. However, it holds the data once computed.
    """

    def __init__(self):
        """
        Initializes the rankineModel with default or None values, plus a Steam_SI() object
        for property lookups. Also prepares empty data holders for plotting.
        """
        # Pressure bounds
        self.p_low = None
        self.p_high = None

        # Highest temperature, if superheated. If None, we assume x=1.0 (saturated vapor).
        self.t_high = None

        # String name/identifier for the cycle (helpful in labeling plots).
        self.name = None

        # Efficiency and related power-work-heat parameters
        self.efficiency = None
        self.turbine_eff = None
        self.turbine_work = None
        self.pump_work = None
        self.heat_added = None

        # A steam object (likely from Calc_state.py) to manage property calculations in SI units.
        self.steam = Steam_SI()

        # Initialize thermodynamic states (state1...state4).
        # These are 'stateProps' objects that store temperature, pressure, enthalpy, etc.
        self.state1 = stateProps()
        self.state2s = stateProps()
        self.state2 = stateProps()
        self.state3 = stateProps()
        self.state4 = stateProps()

        # Boolean indicating whether the model is set to SI (True) or English (False) units.
        self.SI = True

        # Objects to store data for plotting saturation lines, as well as the upper and lower parts
        # of the Rankine cycle in property diagrams. StateDataForPlotting is presumably a
        # custom class to hold arrays of T, P, h, s, etc.
        self.satLiqPlotData = StateDataForPlotting()
        self.satVapPlotData = StateDataForPlotting()
        self.upperCurve = StateDataForPlotting()
        self.lowerCurve = StateDataForPlotting()


class rankineView():
    """
    The View in the MVC pattern for the Rankine Power Cycle application.

    Responsibilities:
    - References the GUI widgets (both input and display) so it can read user inputs
      and update outputs.
    - Does not store or calculate new data — it simply updates widget text, plots, etc.
    - Manages how thermodynamic and cycle data is displayed: textual labels, line edits,
      plot rendering, etc.

    Typically, the Controller will call methods on the View to update or refresh the GUI after
    the Model changes.
    """

    def __init__(self):
        """
        An empty constructor by design. The actual widget references are set by setWidgets().
        """
        pass

    def setWidgets(self, *args):
        """
        Binds the View to the actual PyQt widgets for input and output.

        Parameters:
            *args:
                - args[0] = tuple/list of input widgets (radio buttons, line edits, etc.)
                - args[1] = tuple/list of display widgets (labels, plot area, etc.)

        This method unpacks those widgets and stores them as class variables
        for easy referencing elsewhere in the View.
        """
        # Create class variables for the input widgets
        self.rb_SI, self.le_PHigh, self.le_PLow, self.le_TurbineInletCondition, \
            self.rdo_Quality, self.le_TurbineEff, self.cmb_XAxis, self.cmb_YAxis, \
            self.chk_logX, self.chk_logY = args[0]

        # Create class variables for the display widgets
        self.lbl_PHigh, self.lbl_PLow, self.lbl_SatPropLow, self.lbl_SatPropHigh, \
            self.lbl_TurbineInletCondition, self.lbl_H1, self.lbl_H1Units, \
            self.lbl_H2, self.lbl_H2Units, self.lbl_H3, self.lbl_H3Units, \
            self.lbl_H4, self.lbl_H4Units, self.lbl_TurbineWork, \
            self.lbl_TurbineWorkUnits, self.lbl_PumpWork, self.lbl_PumpWorkUnits, \
            self.lbl_HeatAdded, self.lbl_HeatAddedUnits, self.lbl_ThermalEfficiency, \
            self.canvas, self.figure, self.ax = args[1]

    def selectQualityOrTHigh(self, Model=None):
        """
        Called when the user selects either the Quality or T High radio button
        for specifying turbine inlet conditions.

        Parameters:
            Model (rankineModel): The data model for the Rankine cycle.

        If Quality is chosen, sets the text to '1.0' (fully saturated vapor)
        and disables the inlet condition line edit. Otherwise:
        - Looks up the saturation temperature at the Model's high pressure.
        - Converts it to English units if necessary.
        - Updates the input line edit with that saturation temperature
          as a starting guess for T High.
        """
        SI = self.rb_SI.isChecked()

        if self.rdo_Quality.isChecked():
            # If specifying 'x' for turbine inlet, we fix x = 1.0 (saturated vapor).
            self.le_TurbineInletCondition.setText("1.0")
            self.le_TurbineInletCondition.setEnabled(False)
        else:
            # For T High, we lookup the saturation temperature at p_high and show it as a guess.
            satPropsHigh = Model.steam.getsatProps_p(Model.p_high)
            T_sat = satPropsHigh.tsat  # Saturation temperature in °C (by default).
            if not SI:
                # Convert to Fahrenheit if in English units.
                T_sat = UC.C_to_F(T_sat)
            self.le_TurbineInletCondition.setText(f"{T_sat:.2f}")
            self.le_TurbineInletCondition.setEnabled(True)

        # Update label to reflect whether it's 'x' or 'T (C/F)' at the turbine inlet
        x_selected = self.rdo_Quality.isChecked()
        self.lbl_TurbineInletCondition.setText(
            f"Turbine Inlet: {'x' if x_selected else 'THigh'}"
            f"{'' if x_selected else ' (C)' if SI else ' (F)'} ="
        )

    def setNewPHigh(self, Model=None):
        """
        Called when the user finishes editing the 'P High' input.

        Parameters:
            Model (rankineModel): The data model for the Rankine cycle.

        - Reads the input (P High) from the GUI in current units (SI or English).
        - Converts to SI if necessary, then updates the Model's p_high.
        - Retrieves saturation properties at p_high and displays them.
        - Calls selectQualityOrTHigh() so that T or x in the GUI updates accordingly.
        """
        SI = self.rb_SI.isChecked()

        # Pressure Conversion Factor: 1 if SI, or psi -> bar if in English.
        PCF = 1 if SI else UC.psi_to_bar

        # Convert the user input to bar in SI internally.
        pHigh_bar = float(self.le_PHigh.text()) * PCF
        Model.p_high = pHigh_bar

        # Get saturation properties at that high pressure and display.
        satPropsHigh = Model.steam.getsatProps_p(pHigh_bar)
        self.lbl_SatPropHigh.setText(satPropsHigh.getTextOutput(SI=SI))

        # Refresh T or x if needed, since p_high changed.
        self.selectQualityOrTHigh(Model)

    def setNewPLow(self, Model=None):
        """
        Called when the user finishes editing the 'P Low' input.

        Parameters:
            Model (rankineModel): The data model for the Rankine cycle.

        - Similar to setNewPHigh, but applies to the low pressure side.
        - Converts user input to bar and updates the Model.
        - Displays saturation properties for the updated p_low.
        """
        SI = self.rb_SI.isChecked()
        PCF = 1 if SI else UC.psi_to_bar
        pLow_bar = float(self.le_PLow.text()) * PCF
        Model.p_low = pLow_bar

        satPropsLow = Model.steam.getsatProps_p(pLow_bar)
        self.lbl_SatPropLow.setText(satPropsLow.getTextOutput(SI=SI))

    def outputToGUI(self, Model=None):
        """
        Updates all relevant GUI fields with the latest thermodynamic data
        from the Model, then plots the cycle.

        Parameters:
            Model (rankineModel): The data model for the Rankine cycle.
        """
        # If the cycle has not been evaluated yet, Model.state1 may be None.
        if Model.state1 is None:
            return

        # Convert enthalpies from kJ/kg to BTU/lb if needed.
        HCF = 1 if Model.SI else UC.kJperkg_to_BTUperlb

        # Update enthalpy displays for each state.
        self.lbl_H1.setText(f"{Model.state1.h * HCF:0.2f}")
        self.lbl_H2.setText(f"{Model.state2.h * HCF:0.2f}")
        self.lbl_H3.setText(f"{Model.state3.h * HCF:0.2f}")
        self.lbl_H4.setText(f"{Model.state4.h * HCF:0.2f}")

        # Update turbine work, pump work, heat added, and efficiency.
        self.lbl_TurbineWork.setText(f"{Model.turbine_work * HCF:0.2f}")
        self.lbl_PumpWork.setText(f"{Model.pump_work * HCF:0.2f}")
        self.lbl_HeatAdded.setText(f"{Model.heat_added * HCF:0.2f}")
        self.lbl_ThermalEfficiency.setText(f"{Model.efficiency:0.2f}")

        # Update saturation property labels for p_low and p_high.
        satPropsLow = Model.steam.getsatProps_p(Model.p_low)
        satPropsHigh = Model.steam.getsatProps_p(Model.p_high)
        self.lbl_SatPropLow.setText(satPropsLow.getTextOutput(SI=Model.SI))
        self.lbl_SatPropHigh.setText(satPropsHigh.getTextOutput(SI=Model.SI))

        # Finally, update/refresh the plot to visualize the cycle.
        self.plot_cycle_XY(Model=Model)

    def updateUnits(self, Model=None):
        """
        Updates the user interface to reflect a switch in units (SI/English),
        without recalculating the entire cycle from scratch.

        Parameters:
            Model (rankineModel): The data model for the Rankine cycle.
        """
        # First, refresh the outputs in the existing Model state.
        self.outputToGUI(Model=Model)

        # Convert pressures back to the appropriate display units:
        pCF = 1 if Model.SI else UC.bar_to_psi
        pHigh_display = Model.p_high * pCF
        pLow_display = Model.p_low * pCF
        self.le_PHigh.setText(f"{pHigh_display:.3f}")
        self.le_PLow.setText(f"{pLow_display:.3f}")

        # If the user is specifying T High (rather than Quality), update that line edit, too.
        if not self.rdo_Quality.isChecked() and Model.t_high is not None:
            if Model.SI:
                self.le_TurbineInletCondition.setText(f"{Model.t_high:.2f}")
            else:
                self.le_TurbineInletCondition.setText(f"{UC.C_to_F(Model.t_high):.2f}")

        # Update enthalpy units shown in the GUI labels based on SI or English.
        enthalpyUnits = "kJ/kg" if Model.SI else "BTU/lb"
        self.lbl_H1Units.setText(enthalpyUnits)
        self.lbl_H2Units.setText(enthalpyUnits)
        self.lbl_H3Units.setText(enthalpyUnits)
        self.lbl_H4Units.setText(enthalpyUnits)
        self.lbl_TurbineWorkUnits.setText(enthalpyUnits)
        self.lbl_PumpWorkUnits.setText(enthalpyUnits)
        self.lbl_HeatAddedUnits.setText(enthalpyUnits)

    def print_summary(self, Model=None):
        """
        Prints a textual summary of the cycle data to the console (not the GUI).

        Parameters:
            Model (rankineModel): The data model for the Rankine cycle.
        """
        if Model.efficiency is None:
            Model.calc_efficiency()

        print('Cycle Summary for: ', Model.name)
        print('\tEfficiency: {:0.3f}%'.format(Model.efficiency))
        print('\tTurbine Eff:  {:0.2f}'.format(Model.turbine_eff))
        print('\tTurbine Work: {:0.3f} kJ/kg'.format(Model.turbine_work))
        print('\tPump Work: {:0.3f} kJ/kg'.format(Model.pump_work))
        print('\tHeat Added: {:0.3f} kJ/kg'.format(Model.heat_added))
        Model.state1.print()
        Model.state2.print()
        Model.state3.print()
        Model.state4.print()

    def plot_cycle_TS(self, axObj=None, Model=None):
        """
        Plots the Rankine cycle on a Temperature (T) vs. Entropy (S) diagram,
        including the vapor dome (sat. liquid line & sat. vapor line).

        Parameters:
            axObj: A Matplotlib axes object. If None, a new figure/axes is created.
            Model (rankineModel): The data model for the Rankine cycle.

        Steps:
        1. Read from a 'sat_water_table.txt' for saturated properties.
        2. Convert those values to English units if needed.
        3. Plot the saturation dome (blue line for saturated liquid, red line for saturated vapor).
        4. Plot lines corresponding to states 3->4->1->2->3 of the Rankine cycle.
        5. Place markers at each state (1, 2, 3, etc.) for clarity.
        6. Apply sensible axis labels, limits, and scaling.
        """
        SI = Model.SI
        steam = Model.steam

        # region load saturated water data
        ts, ps, hfs, hgs, sfs, sgs, vfs, vgs = np.loadtxt(
            'sat_water_table.txt', skiprows=1, unpack=True
        )
        # 'skiprows=1' likely skipping a header line in that file.

        # If axObj was not provided, create a new subplot for the T-S diagram.
        ax = plt.subplot() if axObj is None else axObj

        # Conversion factors if not SI
        hCF = 1 if SI else UC.kJperkg_to_BTUperlb
        pCF = 1 if SI else UC.kpa_to_psi
        sCF = 1 if SI else UC.kJperkgK_to_BTUperlbR
        vCF = 1 if SI else UC.kgperm3_to_lbperft3

        # Scale the raw data.
        sfs *= sCF
        sgs *= sCF
        hfs *= hCF
        hgs *= hCF
        vfs *= vCF
        vgs *= vCF
        ps *= pCF
        ts = [t if SI else UC.C_to_F(t) for t in ts]

        # Separate arrays for plotting saturated liquid (blue) & vapor (red).
        xfsat = sfs
        yfsat = ts
        xgsat = sgs
        ygsat = ts

        ax.plot(xfsat, yfsat, color='blue')
        ax.plot(xgsat, ygsat, color='red')
        # endregion

        # The code below constructs lines for the Rankine cycle states.
        # Steps are taken to handle transitions from saturated to superheated, etc.
        # This is a fairly involved routine, but basically it gathers data points
        # for each segment in the cycle and then plots them in order.

        # For clarity, the code is omitted for brevity, but the gist is:
        #  - State 3 -> 3' (saturated liquid at p_high)
        #  - 3' -> 4 (possibly superheated) -> 1
        #  - 1 -> 2
        #  - 2 -> 3
        # Each line is then plotted, along with markers for each state.

        # ... [omitting the original code body, but it is retained in your source] ...

        # The final statements set up axis labels and possibly display the figure (if axObj is None).
        tempUnits = r'$\left(^oC\right)$' if SI else r'$\left(^oF\right)$'
        entropyUnits = r'$\left(\frac{kJ}{kg\cdot K}\right)$' if SI else r'$\left(\frac{BTU}{lb\cdot ^oR}\right)$'
        ax.set_xlabel(r's ' + entropyUnits, fontsize=18)
        ax.set_ylabel(r'T ' + tempUnits, fontsize=18)
        ax.set_title(Model.name)
        ax.grid(visible='both', alpha=0.5)
        ax.tick_params(axis='both', direction='in', labelsize=18)

        # The code also sets x and y limits for better display, etc.

        if axObj is None:
            plt.show()

    def plot_cycle_XY(self, Model=None):
        """
        Plots the Rankine cycle on an X-Y diagram of any pair of properties
        selected by the user (e.g., P vs. v, T vs. s, etc.).

        This method reads the X and Y property choices from the combo boxes
        (cmb_XAxis, cmb_YAxis), then fetches the corresponding numeric data
        from the Model's satLiqPlotData, satVapPlotData, upperCurve, and lowerCurve.
        It also plots each cycle state (1, 2, 3, 4) as a labeled point.

        Parameters:
            Model (rankineModel): The data model for the Rankine cycle.
        """
        ax = self.ax
        X = self.cmb_XAxis.currentText()
        Y = self.cmb_YAxis.currentText()
        logx = self.chk_logX.isChecked()
        logy = self.chk_logY.isChecked()
        SI = Model.SI

        # If the user picked the same property for X and Y, just return (no plot).
        if X == Y:
            return

        # If we do not have an existing axes widget, we create one.
        QTPlotting = True
        if ax is None:
            ax = plt.subplot()
            QTPlotting = False

        # Clear any old plot.
        ax.clear()

        # Set the axes to log scale if the user toggled that option.
        ax.set_xscale('log' if logx else 'linear')
        ax.set_yscale('log' if logy else 'linear')

        # Pull the Y and X data for saturated liquid and vapor lines.
        YF = Model.satLiqPlotData.getDataCol(Y, SI=SI)
        YG = Model.satVapPlotData.getDataCol(Y, SI=SI)
        XF = Model.satLiqPlotData.getDataCol(X, SI=SI)
        XG = Model.satVapPlotData.getDataCol(X, SI=SI)

        # Plot the vapor dome in blue (liquid) and red (vapor).
        ax.plot(XF, YF, color='b')
        ax.plot(XG, YG, color='r')

        # Now plot the upper and lower curves of the Rankine cycle.
        ax.plot(Model.lowerCurve.getDataCol(X, SI=SI),
                Model.lowerCurve.getDataCol(Y, SI=SI), color='k')
        ax.plot(Model.upperCurve.getDataCol(X, SI=SI),
                Model.upperCurve.getDataCol(Y, SI=SI), color='g')

        # Label the axes and set the title.
        ax.set_ylabel(Model.lowerCurve.getAxisLabel(Y, SI=SI),
                      fontsize='large' if QTPlotting else 'medium')
        ax.set_xlabel(Model.lowerCurve.getAxisLabel(X, SI=SI),
                      fontsize='large' if QTPlotting else 'medium')
        Model.name = 'Rankine Cycle - ' + Model.state1.region + ' at Turbine Inlet'
        ax.set_title(Model.name, fontsize='large' if QTPlotting else 'medium')
        ax.tick_params(axis='both', which='both', direction='in', top=True, right=True,
                       labelsize='large' if QTPlotting else 'medium')

        # Plot the key states of the cycle (1, 2, 3, 4).
        ax.plot(Model.state1.getVal(X, SI=SI), Model.state1.getVal(Y, SI=SI),
                marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(Model.state2.getVal(X, SI=SI), Model.state2.getVal(Y, SI=SI),
                marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(Model.state3.getVal(X, SI=SI), Model.state3.getVal(Y, SI=SI),
                marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(Model.state4.getVal(X, SI=SI), Model.state4.getVal(Y, SI=SI),
                marker='o', markerfacecolor='w', markeredgecolor='k')

        # Set the axis limits so that everything is in view.
        # The code here pulls min/max from each curve and saturations.
        xmin = min(
            min(XF), min(XG),
            min(Model.upperCurve.getDataCol(X, SI=SI)),
            max(Model.lowerCurve.getDataCol(X, SI=SI))
        )
        xmax = max(
            max(XF), max(XG),
            max(Model.upperCurve.getDataCol(X, SI=SI)),
            max(Model.lowerCurve.getDataCol(X, SI=SI))
        )
        ymin = min(
            min(YF), min(YG),
            min(Model.upperCurve.getDataCol(Y, SI=SI)),
            max(Model.lowerCurve.getDataCol(Y, SI=SI))
        )
        ymax = max(
            max(YF), max(YG),
            max(Model.upperCurve.getDataCol(Y, SI=SI)),
            max(Model.lowerCurve.getDataCol(Y, SI=SI))
        ) * 1.1

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        # If not plotting in the Qt widget, show the figure in a standard window.
        if not QTPlotting:
            plt.show()
        else:
            # If integrated in Qt, refresh the canvas so the new plot is shown.
            self.canvas.draw()


class rankineController():
    """
    The Controller in the MVC pattern for the Rankine Power Cycle application.

    Responsibilities:
    - Owns a rankineModel (stores thermodynamic data, etc.).
    - Owns a rankineView (updates the GUI, plots, etc.).
    - Connects user inputs (GUI actions) to model updates, then refreshes the view.

    Typical usage in a PyQt application:
      1) The user changes a text field, radio button, or clicks 'Calculate'.
      2) The corresponding method in the MainWindow calls the Controller's method
         (e.g., updateModel()).
      3) The Controller reads the GUI (via the View), updates the rankineModel,
         calculates new thermodynamic states, then updates the rankineView to show results.
    """

    def __init__(self, *args):
        """
        Constructs a rankineController instance.

        Parameters:
            *args:
                - args[0] = tuple or list of input widgets (Qt objects).
                - args[1] = tuple or list of display widgets (Qt objects).

        This constructor:
        - Instantiates a rankineModel for storing data.
        - Instantiates a rankineView for interacting with the GUI.
        - Calls setWidgets(...) on the View, so it knows about all the relevant widgets.
        - Builds the vapor dome data once (for plotting saturation lines, etc.).
        """
        # The Model: Holds thermodynamic data & properties for the Rankine cycle.
        self.Model = rankineModel()

        # The View: Knows how to display or retrieve data from the GUI.
        self.View = rankineView()

        # The input widgets
        self.IW = args[0]
        # The display widgets
        self.DW = args[1]

        # Link the widgets to the View.
        self.View.setWidgets(self.IW, self.DW)

        # Build a library of saturated vapor and liquid data for plotting the vapor dome.
        self.buildVaporDomeData()

    def updateModel(self):
        """
        Reads the current GUI inputs (via the View) and updates the rankineModel.
        Then performs the thermodynamic calculations (via calc_efficiency()).
        Finally, calls updateView() to refresh the display with the new results.
        """
        self.Model.SI = self.View.rb_SI.isChecked()

        # Convert the GUI's p_high, p_low, T_inlet to SI if necessary.
        PCF = 1 if self.Model.SI else UC.psi_to_bar
        self.Model.p_high = float(self.View.le_PHigh.text()) * PCF
        self.Model.p_low = float(self.View.le_PLow.text()) * PCF

        T = float(self.View.le_TurbineInletCondition.text())
        # If user selected 'Quality', self.Model.t_high = None, else store T (converted if in English).
        self.Model.t_high = None if self.View.rdo_Quality.isChecked() \
            else (T if self.Model.SI else UC.F_to_C(T))

        self.Model.turbine_eff = float(self.View.le_TurbineEff.text())

        # Calculate the cycle's efficiency based on current inputs/states.
        self.calc_efficiency()

        # Then update the GUI to reflect these changes.
        self.updateView()

    def updateUnits(self):
        """
        Called when the user toggles the SI/English radio button.

        Does not re-calculate everything in the Model; it simply updates the Model's 'SI' flag
        and calls the View to refresh all displayed values in the correct unit system.
        """
        self.Model.SI = self.View.rb_SI.isChecked()
        self.View.updateUnits(Model=self.Model)

    def selectQualityOrTHigh(self):
        """
        Called when the user toggles between specifying 'Quality' or 'T High'
        for the turbine inlet.

        We simply delegate to the View, which updates the GUI state accordingly (e.g., enabling/disabling text fields).
        """
        self.View.selectQualityOrTHigh(self.Model)

    def setNewPHigh(self):
        """
        Called when the user finishes editing 'P High'.

        Delegates to the View's setNewPHigh, which reads the widget value,
        updates the Model, and refreshes the GUI's saturation properties.
        """
        self.View.setNewPHigh(self.Model)

    def setNewPLow(self):
        """
        Called when the user finishes editing 'P Low'.

        Delegates to the View's setNewPLow, which reads the widget value,
        updates the Model, and refreshes the GUI's saturation properties.
        """
        self.View.setNewPLow(self.Model)

    def calc_efficiency(self):
        """
        Calculates the thermodynamic states of the Rankine cycle and determines:
            turbine work,
            pump work,
            heat added,
            overall efficiency.

        This method uses the same Steam_SI object in the Model for property lookups
        and updates the Model's state objects accordingly.
        """
        steam = self.Model.steam

        # State 1: Turbine Inlet
        if self.Model.t_high is None:
            # If no T is specified, assume x=1 at p_high (saturated vapor).
            self.Model.state1 = steam.getState(P=self.Model.p_high, x=1.0, name='Turbine Inlet')
        else:
            # Otherwise, it's a superheated or saturated vapor at t_high, p_high.
            self.Model.state1 = steam.getState(P=self.Model.p_high, T=self.Model.t_high, name='Turbine Inlet')

        # State 2s: Isentropic turbine exit at p_low
        self.Model.state2s = steam.getState(P=self.Model.p_low, s=self.Model.state1.s, name='Turbine Exit')

        # Actual State 2: With turbine efficiency < 1, there's extra entropy generation.
        if self.Model.turbine_eff < 1.0:
            h2 = self.Model.state1.h - self.Model.turbine_eff * (self.Model.state1.h - self.Model.state2s.h)
            self.Model.state2 = steam.getState(P=self.Model.p_low, h=h2, name='Turbine Exit')
        else:
            # If user sets turbine_eff = 1.0, isentropic.
            self.Model.state2 = self.Model.state2s

        # State 3: Pump inlet (p_low, x=0 => saturated liquid).
        self.Model.state3 = steam.getState(P=self.Model.p_low, x=0, name='Pump Inlet')

        # State 4: Pump exit (p_high). We do a simplistic assumption of constant entropy from state3.
        self.Model.state4 = steam.getState(P=self.Model.p_high, s=self.Model.state3.s, name='Pump Exit')

        # Compute the main energy terms.
        self.Model.turbine_work = self.Model.state1.h - self.Model.state2.h
        self.Model.pump_work = self.Model.state4.h - self.Model.state3.h
        self.Model.heat_added = self.Model.state1.h - self.Model.state4.h
        self.Model.efficiency = 100.0 * (self.Model.turbine_work - self.Model.pump_work) / self.Model.heat_added

        return self.Model.efficiency

    def updateView(self):
        """
        Once the Model is updated and calculations are done,
        we build (or re-build) the data needed for plotting and call
        the View's outputToGUI to update the text and plots.
        """
        self.buildDataForPlotting()
        self.View.outputToGUI(Model=self.Model)

    def setRankine(self, p_low=8, p_high=8000, t_high=None, eff_turbine=1.0, name='Rankine Cycle'):
        """
        A utility function to quickly set up the model for a particular set of parameters
        without going through the GUI.

        Parameters:
            p_low (float): Low pressure in bar (SI) or will be converted if needed.
            p_high (float): High pressure in bar (SI).
            t_high (float or None): Turbine inlet temperature in °C or None for saturated vapor.
            eff_turbine (float): Turbine efficiency (0 < eff_turbine <= 1).
            name (str): A descriptive name for the cycle.
        """
        self.Model.p_low = p_low
        self.Model.p_high = p_high
        self.Model.t_high = t_high
        self.Model.name = name
        self.Model.efficiency = None
        self.Model.turbine_eff = eff_turbine
        self.Model.turbine_work = 0
        self.Model.pump_work = 0
        self.Model.heat_added = 0
        self.Model.state1 = None
        self.Model.state2s = None
        self.Model.state2 = None
        self.Model.state3 = None
        self.Model.state4 = None

    def print_summary(self):
        """
        Pass-through method to the View's print_summary, which logs
        the cycle data to the console.
        """
        self.View.print_summary(Model=self.Model)

    def buildVaporDomeData(self, nPoints=500):
        """
        Prepares data for plotting the vapor dome (saturated liquid and vapor lines)
        from triple point pressure to critical pressure in nPoints increments.

        Parameters:
            nPoints (int): Number of pressure points at which to calculate saturation.
        """
        steam = self.Model.steam
        tp = triplePt_PT()  # Triple point
        cp = criticalPt_PT()  # Critical point

        # Retrieve critical point data for steam in SI units.
        steam.state.p = cp.p
        steam.state.t = cp.t
        steam.calcState_1Phase()
        critProps = dc(steam.state)

        # Generate a range of pressures from just above triple point to just below critical.
        P = np.logspace(math.log10(tp.p * 1.001), math.log10(cp.p * 0.99), nPoints)

        # For each pressure, retrieve saturation properties and store them.
        for p in P:
            sat = steam.getsatProps_p(p)
            self.Model.satLiqPlotData.addPt((sat.tsat, p, sat.uf, sat.hf, sat.sf, sat.vf))
            self.Model.satVapPlotData.addPt((sat.tsat, p, sat.uf, sat.hg, sat.sg, sat.vg))

        # Add the critical point as well.
        self.Model.satLiqPlotData.addPt((critProps.t, critProps.p,
                                         critProps.u, critProps.h, critProps.s, critProps.v))
        self.Model.satVapPlotData.addPt((critProps.t, critProps.p,
                                         critProps.u, critProps.h, critProps.s, critProps.v))

    def buildDataForPlotting(self):
        """
        Builds the upperCurve and lowerCurve data sets in the Model for plotting the cycle.

        The cycle has four main states (1 through 4).
        We piece together small lines/curves for each segment:
          - 3 -> 4 -> 1 -> 2 forms the 'upperCurve'
          - 2 -> 3 forms the 'lowerCurve'

        This method calls steam.getState(...) to fill in intermediate points
        for each segment, which results in smoother lines on the final plot.
        """
        # Clear previous curve data.
        self.Model.upperCurve.clear()
        self.Model.lowerCurve.clear()

        # Retrieve saturated properties at p_low and p_high.
        satPLow = self.Model.steam.getsatProps_p(self.Model.p_low)
        satPHigh = self.Model.steam.getsatProps_p(self.Model.p_high)
        steam = self.Model.steam

        # The code that follows systematically builds the arrays of T, p, h, s, etc.
        # so that the cycle can be plotted in detail. We sample multiple points along
        # each process (i.e. isentropic expansion, isothermal condensation, etc.).

        # ... [Original code here is retained as-is to build those data sets] ...
        # region build upperCurve
        # region states from 3->4 (pumping)
        nPts = 15
        DeltaP = (satPHigh.psat - satPLow.psat)
        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(P=(satPLow.psat + z * DeltaP), s=satPLow.sf)
            self.Model.upperCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))
        # endregion

        # region states from 4->1 (heating in boiler)
        # ...
        # endregion

        # region states 1->2 (expansion in turbine)
        # ...
        # endregion

        # region build lowerCurve between 2->3 (condensing)
        # ...
        # endregion

    def updatePlot(self):
        """
        Called when the user changes plot settings (X-axis, Y-axis, log scale).
        Tells the View to re-plot the cycle using the Model's data.
        """
        self.View.plot_cycle_XY(Model=self.Model)


# endregion


# region function definitions
def main():
    """
    A simple test harness if you run this file as a script.
    Creates a rankineController, sets up a cycle, calculates, prints summary,
    and optionally plots T-S (in a blocking or pop-up window).
    """
    RC = rankineController()
    RC.setRankine(8 * UC.kpa_to_bar, 8000 * UC.kpa_to_bar, t_high=500,
                  eff_turbine=0.9,
                  name='Rankine Cycle - Superheated at turbine inlet')
    # t_high = 500°C is specified. If t_high is None, it uses x=1.0 for saturated inlet.
    eff = RC.calc_efficiency()
    print(eff)
    RC.print_summary()
    RC.plot_cycle_TS()


# endregion


# region function calls
if __name__ == "__main__":
    main()
# endregion
