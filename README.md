CACTUS
======
CACTUS: Code for Axial- and Cross-flow TUrbine Simulation

Copyright (c) 2013, Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. 


Building and installing on Linux
================================
Also see `install.txt`

Ubuntu
------
  * Install gfortran if not already installed: 
    * From a terminal `sudo apt-get update`
    * `sudo apt-get install gfortran`
  * Install LAPACK and BLAS:
	* `sudo apt-get install liblapack-dev`
	* `sudo apt-get install libblas-dev`
  * Inside this directory run `make -f makefile.gfortran`
  * Locate the cactus executable in the directory above, and place inside a directory on the
    system path.
    * `gedit ~/.bashrc` to add a folder to the system path.


Building and installing on Windows
==================================
  * Download and install MinGW (32-bit) and MSYS (http://mingw.org)
  * Open MinGW bash in this directory and execute `make`
  * Execute `make install` (Note: this assumes default MinGW path)


Usage
=====
Linux
-----
  * cd into project director

Windows
-------
  * Open MinGW bash in project directory (where the `*.in` file is located)
  * Run `cactus project_name.in`
  
  
Python interface
================
An experimental Python 2.7 interface is provided in the `cactus_py` directory. 
This package provides functions for reading CACTUS data and a class for running a batch 
of CACTUS runs to create a performance curve. An example script is located in
`Test/TestCase1/testhawt.py`

Installation
------------
Inside the `cactus_py` directory, run `python setup.py install`

Dependencies
------------
  * Numpy
  * matplotlib
  * pandas
  
Usage
-----
See `Test/TestCase1/testhawt.py` for a simple example script.


3rd party libraries: LAPACK and BLAS
====================================
Included in the `lib` directory. Prebuilt for MinGW-Win32.
http://icl.cs.utk.edu/lapack-for-windows/lapack/#libraries_mingw