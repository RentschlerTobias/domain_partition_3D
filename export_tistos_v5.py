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
labels = bV.labels()
print(labels)

# Try to export mesh
mesh = bV.get("ruWithRounding_mechMesh")
print("\nMesh label:", mesh.getLabel())

# Try to export mesh
print("\nTrying to export mesh...")
try:
    mesh.dumpGrid("/work/tistos_output.msh")
    print("Mesh dumped successfully")
except Exception as e:
    print(f"Error dumping mesh: {e}")
