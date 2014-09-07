#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script calculates the performance curve of the test VAWT
"""

from cactus import PerfCurve

pc = PerfCurve("TestVAWT2")
pc.basecase.geomfile = "../TestGeom/unh-rvat.geom"
pc.basecase.rho = 1000.0
pc.basecase.viscosity = 1.0e-6
pc.run(1.0, 2.5, 0.25)
pc.plot()
