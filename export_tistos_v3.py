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

# Check available models
print("Available bounded volumes:")
for ii in range(bV.size()):
    label = bV.label(ii)
    print(f"  {label}")

# Try to export mesh using dtOO's built-in functions
mesh = bV.get("ruWithRounding_mechMesh")
print("\nMesh label:", mesh.getLabel())
print("Mesh type:", type(mesh))

# Try to export
mesh.dumpGrid("/work/tistos_output.msh")
print("Mesh dumped to /work/tistos_output.msh")
