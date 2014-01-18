# -*- coding: utf-8 -*-
"""
This module contains functions for running CACTUS and manipulating output data

"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import os
import json
import geom

# Translation dictionary for quantities printed in CSV file
quantities = {"t" : " Normalized Time (-)",
              "theta" : "Theta (rad)",
              "rev" : "Rev",
              "ct" : "Torque Coeff. (-)",
              "cp" : "Power Coeff. (-)",
              "cd" : "Fx Coeff. (-)"}
# Default configuration parameters
config_defaults = {"GPFlag" : 0,
                   "nr" : 10,
                   "nti" : 16,
                   "convrg" : 0.0001,
                   "iut" : 0,
                   "ifc" : 0,
                   "ixterm" : 0,
                   "ntif" : 16,
                   "iutf" : 1,
                   "nric" : 9,
                   "convrgf" : 0.0001,
                   "DiagOutFlag" : 1,
                   "Output_ELFlag" : 1}
# Default case inputs 
case_defaults = {"jbtitle" : "TestHAWT",
                 "rho" : 0.002378,
                 "vis" : 0.3739E-6,
                 "tempr" : 60.0,
                 "hBLRef" : 56.57,
                 "slex" : 0.0,
                 "hAG" : 30.0,
                 "RPM" : 70,
                 "Ut" : 4,
                 "GeomFilePath" : "../TestGeom/TestHAWT.geom",
                 "nSect" : 1,
                 "AFDPath" : "../../Airfoil_Section_Data/NACA_0015.dat"} 

def read_time_data(casename):
    """Load csv time data."""
    with open(casename + "_TimeData.csv", "r") as f:
        data = pd.read_csv(f)
    return data
    
def read_rev_data(casename):
    """Load csv time data."""
    with open(casename + "_RevData.csv", "r") as f:
        data = pd.read_csv(f)
    return data
    
def plot_time_data(casename, quantity):
    data = read_time_data(casename)
    t = data[quantities["t"]]
    q = data[quantities[quantity]]
    plt.figure()
    plt.plot(t, q)
    plt.show()
    
def get_mean_cp(casename):
    data = read_rev_data(casename)
    cp = data["Power Coeff. (-)"]
    meancp = np.mean(cp[len(cp)/2:])
    return meancp
    
def write_input_file(name, casedict, configdict):
    casedict["jbtitle"] = name
    with open(name+".in", "w") as f:
        f.write("&ConfigInputs\n")
        for key, value in configdict.iteritems():
            f.write("\t" + key + " = " + str(value) + "\n")
        f.write("/End\n")
        f.write("&CaseInputs\n")
        for key, value in casedict.iteritems():
            if "title" in key or "Path" in key:
                f.write("\t" + key + " = '" + str(value) + "'\n")
            else:
                f.write("\t" + key + " = " + str(value) + "\n")
        f.write("/End\n")
        
def runcase(name):
    """Run a case. Currently assumes that if the user is running Windows, 
    MinGW is installed in the default configuration. If running Linux,
    it is assumed CACTUS is on the path."""
    if os.name == "nt":
        cm = 'C:/MinGW/msys/1.0/bin/bash.exe --login -c "cd \''
        cm += os.getcwd().replace('\\', '/') 
        cm += '\' ; cactus ' + name + ".in\""
    elif os.name == "posix":
        cm = "cactus " + name + ".in"
    subprocess.call(cm)
    

class Case(object):
    """
    A CACTUS case object for running a single simulation.
    
    Default geometry files will be same as case name and in same directory.
    
    If an input file with the case name exists in the working directory it
    will not be loaded unless `loadconfig=True`.
    """
    def __init__(self, name, loadconfig=False):
        self.name = name
        # Default numerical configuration parameters
        self.GPFlag = 0
        self.nr = 10
        self.nti = 16
        self.convrg = 0.0001
        self.iut = 0
        self.ifc = 0
        self.ixterm = 0
        self.ntif = 16
        self.iutf = 1
        self.nric = 9
        self.convrgf = 0.0001
        self.diag_out_flag = 1
        self.output_el_flag = 1
        # Default case configuration parameters
        self.rho = 0.002378
        self.viscosity = 0.3739E-6
        self.tempr = 60.0
        self.hBLRef = 56.57
        self.slex = 0.0
        self.hAG = 30.0
        self.rpm = 70.0
        self.tsr = 4.0
        self.nSect = 1
        self.geomfile = name + ".geom"
        self.foil_data_file = "../../Airfoil_Section_Data/NACA_0015.dat"
        self.setconfig()
        if loadconfig: 
            try: 
                self.loadconfig()
            except IOError:
                pass
    def creategeom(self, n_blade, n_belem, n_strut, s_selem, ref_r, rot_n,
                   rot_p, ref_ar, turb_type=None, **kwargs):
        """Create new geometry file."""
        turbgeom = geom.Turbine(n_blade, n_belem, n_strut, s_selem, ref_r, 
                                rot_n, rot_p, ref_ar, turb_type=turb_type, **kwargs)
        turbgeom.writefile(self.name)
    def loadconfig(self):
        """loads a configuration file with the same name as the case, 
        inside the same directory."""
        with open(self.name + ".in", "r") as f:
            for line in f.readlines():
                if "=" in line:
                    i = line.split()[0]
                    v = line.split()[-1]
                    if i == "GPFlag": 
                        self.GPFlag = int(v)
                    elif i == "nr":
                        self.nr = int(v)
                    elif i == "nti":
                        self.nti = int(v)
                    elif i == "convrg":
                        self.convrg = float(v)
                    elif i == "iut":
                        self.iut = int(v)
                    elif i == "ifc":
                        self.ifc = int(v)
                    elif i == "ixterm":
                        self.ixterm = int(v)
                    elif i == "ntif":
                        self.ntif = int(v)
                    elif i == "iutf":
                        self.iutf = int(v)
                    elif i == "nric":
                        self.nric = int(v)
                    elif i == "convrgf":
                        self.convrgf = float(v)
                    elif i == "DiagOutFlag":
                        self.diag_out_flag = int(v)
                    elif i == "Output_ELFlag":
                        self.output_el_flag = int(v)
                    elif i == "rho":
                        self.rho = float(v)
                    elif i == "vis":
                        self.viscosity = float(v)
                    elif i == "tempr":
                        self.tempr = float(v)
                    elif i == "hBLRef":
                        self.hBLRef = float(v)
                    elif i == "slex":
                        self.slex = float(v)
                    elif i == "hAG":
                        self.hAG = float(v)
                    elif i == "RPM":
                        self.rpm = float(v)
                    elif i == "Ut":
                        self.tsr = float(v)
                    elif i == "GeomFilePath":
                        self.geomfile = str(v)
                    elif i == "nSect":
                        self.nSect = int(v)
                    elif i == "AFDPath":
                        self.foil_data_file = str(v)
        self.setconfig()
    def setconfig(self):
        """Sets configuration dictionaries and writes input file."""
        self.numconfig = {"GPFlag" : self.GPFlag,
                          "nr" : self.nr,
                          "nti" : self.nti,
                          "convrg" : self.convrg,
                          "iut" : self.iut,
                          "ifc" : self.ifc,
                          "ixterm" : self.ixterm,
                          "ntif" : self.ntif,
                          "iutf" : self.iutf,
                          "nric" : self.nric,
                          "convrgf" : self.convrgf,
                          "DiagOutFlag" : self.diag_out_flag,
                          "Output_ELFlag" : self.output_el_flag}
        self.caseconfig = {"jbtitle" : self.name,
                           "rho" : self.rho,
                           "vis" : self.viscosity,
                           "tempr" : self.tempr,
                           "hBLRef" : self.hBLRef,
                           "slex" : self.slex,
                           "hAG" : self.hAG,
                           "RPM" : self.rpm,
                           "Ut" : self.tsr,
                           "GeomFilePath" : self.geomfile,
                           "nSect" : 1,
                           "AFDPath" : self.foil_data_file}
    def writeconfig(self):
        self.setconfig()
        with open(self.name+".in", "w") as f:
            f.write("&ConfigInputs\n")
            for key, value in self.numconfig.iteritems():
                f.write("\t" + key + " = " + str(value) + "\n")
            f.write("/End\n")
            f.write("&CaseInputs\n")
            for key, value in self.caseconfig.iteritems():
                if "title" in key or "Path" in key:
                    f.write("\t" + key + " = '" + value + "'\n")
                else:
                    f.write("\t" + key + " = " + str(value) + "\n")
            f.write("/End\n")
    def calc_cp(self):
        """Returns average power coefficient."""
        self.cp = get_mean_cp(self.name)
        return self.cp
    def run(self):
        self.writeconfig()
        if os.name == "nt":
            cm = 'C:/MinGW/msys/1.0/bin/bash.exe --login -c "cd \''
            cm += os.getcwd().replace('\\', '/') 
            cm += '\' ; cactus ' + self.name + ".in\""
        elif os.name == "posix":
            cm = "cactus " + self.name + ".in"
        subprocess.call(cm)
    

class PerfCurve(object):
    """Object that represents a performance curve, or batch runs of cases.
    Settings can be changed through the properties of `self.basecase`."""
    def __init__(self, name):
        self.name = name
        self.basecase = Case(self.name)
    def run(self, tsr_start, tsr_stop, tsr_step):
        # Create empty power coefficient list
        self.cp = []
        self.tsr = np.arange(tsr_start, tsr_stop+tsr_step, tsr_step)
        for tsr in self.tsr:
            case = self.basecase
            case.tsr = tsr
            case.run()
            self.cp.append(case.calc_cp())
        self.lastcase = case
        self.save()
        os.remove(self.name+"_ElementData.csv")
        os.remove(self.name+"_Param.csv")
        os.remove(self.name+"_RevData.csv")
        os.remove(self.name+"_TimeData.csv")
        os.remove(self.name+".in")
    def save(self):
        with open(self.name+"_PerfCurve.json", "w") as f:
            data = {"tsr" : self.tsr.tolist(), "cp" : self.cp,
                    "numconfig_lastcase" : self.lastcase.numconfig, 
                    "caseconfig_lastcase" : self.lastcase.caseconfig}
            f.write(json.dumps(data, indent=2))
    def clean(self):
        """Deletes all files associated with this performance curve."""
        pass
    def plot(self):
        with open(self.name+"_PerfCurve.json") as f:
            data = json.load(f)
        plt.figure()
        plt.plot(data["tsr"], data["cp"], '-o')
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$C_P$")
        plt.show()
    
    
if __name__ == "__main__":
#    case = Case("test")
#    case.geomfile = "../../Test/TestGeom/TestHAWT.geom"
#    case.nr = 4
#    case.creategeom(2, 5, 0, 5, 3.52, [0,0,1], [0,0,0], 31.5, turb_type="HAWT",
#                    **geom.hawt_defaults)
#    case.run()
    pc = PerfCurve("test")
    pc.basecase.geomfile = "../../Test/TestGeom/TestHAWT.geom"
    pc.basecase.nr = 4
    pc.run(4, 8, 0.25)
    pc.plot()