#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script calculates the performance curve of the test VAWT
"""

from cactus import PerfCurve

pc = PerfCurve("TestVAWT2")
pc.basecase.geomfile = "../TestGeom/TestVAWT.geom"
pc.basecase.foil_data_file = "../../Airfoil_Section_Data/NACA_0015.dat"
pc.basecase.rho = 1000.0
pc.basecase.viscosity = 1.0e-6
pc.run(3.0, 7.0, 0.5) # arguments are tsr_start, tsr_stop, tsr_step
pc.plot()
