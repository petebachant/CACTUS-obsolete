# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 19:19:55 2014

@author: Pete Bachant
"""

from cactus import PerfCurve

pc = PerfCurve("TestHAWT2")
pc.run(3.0, 7.0, 0.5)
pc.plot()