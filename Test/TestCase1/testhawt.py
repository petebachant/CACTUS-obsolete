#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This sample script runs a performance curve for the test HAWT case.
"""

from cactus import PerfCurve

pc = PerfCurve("TestHAWT2")
pc.basecase.geomfile = "../TestGeom/TestHAWT.geom"
pc.run(3.0, 7.0, 0.5)
pc.plot()