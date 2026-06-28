#!/usr/bin/env python3
import sys
sys.path.insert(0, '/dtOO-install/tools')

from dtOOPythonSWIG import *

# Set up OpenFOAM environment
import os
os.environ['LD_LIBRARY_PATH'] = '/dtOO-install/lib:/usr/lib/openfoam/openfoam2406/platforms/linux64GccDPInt32Opt/lib/sys-openmpi:/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib/hpc/gnu7/mpi/openmpi/4.1.6/lib64'

# Build the T2_7461 case
logMe.initLog('build.log')
dtXmlParser.init("machine.xml", "T2_7461.xml")
parser = dtXmlParser.reference()
parser.parse()
bC = baseContainer()
cV = labeledVectorHandlingConstValue()
aF = labeledVectorHandlingAnalyticFunction()
aG = labeledVectorHandlingAnalyticGeometry()
bV = labeledVectorHandlingBoundedVolume()
dC = labeledVectorHandlingDtCase()
dP = labeledVectorHandlingDtPlugin()
parser.createConstValue(cV)
parser.loadStateToConst("T2_7461", cV)
parser.destroyAndCreate(bC, cV, aF, aG, bV, dC, dP)
dC.get("tistos_ru_of_n").runCurrentState()
bV.get("ruWithRounding_mechMesh").makeGrid()

# Export the mesh
model = bV.get("ruWithRounding_mechMesh").getModel()
print("Exporting mesh...")
model.writeMSH("/work/T2_7461_ru_gridGmsh.msh")
print("Mesh exported to /work/T2_7461_ru_gridGmsh.msh")

# Also export the grid model
bV.get("ru_gridGmsh").makeGrid()
model_grid = bV.get("ru_gridGmsh").getModel()
model_grid.writeMSH("/work/T2_7461_ru_gridGmsh_grid.msh")
print("Grid mesh exported to /work/T2_7461_ru_gridGmsh_grid.msh")
