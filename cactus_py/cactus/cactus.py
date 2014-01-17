# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 13:42:15 2014
@author: Pete Bachant

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
    def __init__(self, name, rho=0.002378, vis=1e-6, nr=10):
        self.name = name
        self.numconfig = config_defaults
        self.numconfig["nr"] = nr
        self.caseconfig = case_defaults
        self.caseconfig["rho"] = rho
        self.caseconfig["vis"] = vis
        self.cp = []
    def run(self, tsr_start, tsr_stop, tsr_step):
        self.tsr = np.arange(tsr_start, tsr_stop+tsr_step, tsr_step)
        for tsr in self.tsr:
            self.caseconfig["Ut"] = tsr
            write_input_file(self.name, self.caseconfig, self.numconfig)
            runcase(self.name)
            self.cp.append(get_mean_cp(self.name))
        self.save()
    def save(self):
        with open(self.name+"_PerfCurve.json", "w") as f:
            data = {"tsr" : self.tsr.tolist(), "cp" : self.cp}
            f.write(json.dumps(data, indent=2))
    def plot(self):
        with open(self.name+"_PerfCurve.json") as f:
            data = json.load(f)
        plt.figure()
        plt.plot(data["tsr"], data["cp"], '-o')
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$C_P$")
    
    
if __name__ == "__main__":
    os.chdir("../Test/TestCase1")
    plt.close("all")
    pc = PerfCurve("TestHAWT2")
#    pc.run(2.0, 6.0, 0.5)
    pc.plot()