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

# Save mesh
import gmsh
model_name = gmsh.model.list()[0]
print(f"Model name: {model_name}")
gmsh.write('/work/tistos_output.msh')
print("Mesh saved to /work/tistos_output.msh")

# List surfaces
surfaces = gmsh.model.getEntities(2)
print(f"Surfaces: {len(surfaces)}")
for dim, tag in surfaces:
    name = gmsh.model.getEntityName(dim, tag)
    print(f"  Surface {tag}: {name}")
