#!/usr/bin/env python3
import sys
sys.path.insert(0, '/dtOO-install/tools')
import dtOOPythonSWIG as dtOO

# Initialize
dtOO.dtDefaults.init()
print("dtOO init OK")

# Build the T1_9 geometry
bC, cV, aF, aG, bV, dC, dP = dtOO.tistos.build()
print("Build OK")
print("Models:", list(bV.keys()))

# Check if gmsh model exists
import os
os.system("gmsh -info 2>&1 | head -5")
