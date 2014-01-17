# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 19:19:55 2014

@author: Pete Bachant
"""

from cactus import PerfCurve

pc = PerfCurve("TestVAWT2")
pc.GeomFilePath = "../TestGeom/TestVAWT.geom"
pc.rho = 1000.0
pc.vis = 1.0e-6
pc.run(3.0, 7.0, 0.5)
pc.plot()