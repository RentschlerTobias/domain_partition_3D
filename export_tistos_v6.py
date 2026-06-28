#!/usr/bin/env python3
import sys
sys.path.insert(0, '/dtOO-install/tools')
from dtOOPythonSWIG import *

# Initialize
logMe.initLog('build.log')
dtXmlParser.init("machine.xml", "T2_7461.xml")
parser = dtXmlParser.reference()
parser.parse()

# Create containers
bC = baseContainer()
cV = labeledVectorHandlingConstValue()
aF = labeledVectorHandlingAnalyticFunction()
aG = labeledVectorHandlingAnalyticGeometry()
bV = labeledVectorHandlingBoundedVolume()
dC = labeledVectorHandlingDtCase()
dP = labeledVectorHandlingDtPlugin()

# Build
parser.createConstValue(cV)
parser.loadStateToConst("T2_7461", cV)
parser.destroyAndCreate(bC, cV, aF, aG, bV, dC, dP)

# Run case
dC.get("tistos_ru_of_n").runCurrentState()

# Make grid
bV.get("ruWithRounding_mechMesh").makeGrid()

# Check available methods
mesh = bV.get("ruWithRounding_mechMesh")
print("Mesh methods:")
for attr in dir(mesh):
    if not attr.startswith('_'):
        print(f"  {attr}")

# Try to get mesh as MSH
print("\nTrying to get mesh...")
try:
    # Try to write to file using system gmsh
    import os
    os.system("echo 'test' > /work/test.txt")
    print("System test OK")
except Exception as e:
    print(f"Error: {e}")
