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
    
def write_input_file(name, casedict=case_defaults, configdict=config_defaults):
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


class PerfCurve(object):
    def __init__(self, name):
        self.name = name
        self.numconfig = {}
        self.caseconfig = {}
        # Numerical configuration parameters
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
        self.DiagOutFlag = 1
        self.Output_ELFlag = 1
        # Case configuration parameters
        self.rho = 0.002378
        self.vis = 0.3739E-6
        self.tempr = 60.0
        self.hBLRef = 56.57
        self.slex = 0.0
        self.hAG = 30.0
        self.RPM = 70
        self.GeomFilePath = "../TestGeom/TestHAWT.geom"
        self.nSect = 1
        self.AFDPath = "../../Airfoil_Section_Data/NACA_0015.dat"
    def setconfig(self):
        self.numconfig["GPFlag"] = self.GPFlag
        self.numconfig["nr"] = self.nr
        self.numconfig["nti"] = self.nti
        self.numconfig["convrg"] = self.convrg
        self.numconfig["iut"] = self.iut
        self.numconfig["ifc"] = self.ifc
        self.numconfig["ixterm"] = self.ixterm
        self.numconfig["ntif"] = self.ntif
        self.numconfig["iutf"] = self.iutf
        self.numconfig["nric"] = self.nric
        self.numconfig["convrgf"] = self.convrgf
        self.numconfig["DiagOutFlag"] = self.DiagOutFlag
        self.numconfig["Output_ELFlag"] = self.Output_ELFlag
        self.caseconfig["jbtitle"] = self.name
        self.caseconfig["rho"] = self.rho
        self.caseconfig["vis"] = self.vis
        self.caseconfig["tempr"] = self.tempr
        self.caseconfig["hBLRef"] = self.hBLRef
        self.caseconfig["slex"] = self.slex
        self.caseconfig["hAG"] = self.hAG
        self.caseconfig["RPM"] = self.RPM
        self.caseconfig["GeomFilePath"] = self.GeomFilePath
        self.caseconfig["nSect"] = self.nSect
        self.caseconfig["AFDPath"] = self.AFDPath
    def run(self, tsr_start, tsr_stop, tsr_step):
        # Create empty power coefficient list
        self.cp = []
        self.setconfig()
        self.tsr = np.arange(tsr_start, tsr_stop+tsr_step, tsr_step)
        for tsr in self.tsr:
            self.caseconfig["Ut"] = tsr
            write_input_file(self.name, self.caseconfig, self.numconfig)
            runcase(self.name)
            self.cp.append(get_mean_cp(self.name))
        self.save()
    def save(self):
        with open(self.name+"_PerfCurve.json", "w") as f:
            data = {"tsr" : self.tsr.tolist(), "cp" : self.cp,
                    "numconfig" : self.numconfig, 
                    "caseconfig" : self.caseconfig}
            f.write(json.dumps(data, indent=2))
    def plot(self):
        with open(self.name+"_PerfCurve.json") as f:
            data = json.load(f)
        plt.figure()
        plt.plot(data["tsr"], data["cp"], '-o')
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$C_P$")
        plt.show()
    
    
if __name__ == "__main__":
    os.chdir("../../Test/TestCase2")
    plt.close("all")
    pc = PerfCurve("TestVAWT2")
    pc.GeomFilePath = "../TestGeom/TestVAWT.geom"
    pc.run(2.0, 4.0, 0.5)
    pc.plot()
