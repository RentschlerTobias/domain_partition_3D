#!/usr/bin/env python3
"""
Extract coarse block structure from T2_7461 hub/shroud surfaces
using dtOO Docker container + Gmsh quasi-structured quad algorithm.
"""
import sys
import os
sys.path.insert(0, '/dtOO-install/tools')

import gmsh
import numpy as np
from dtOOPythonSWIG import *

# Set environment
os.environ['LD_LIBRARY_PATH'] = '/dtOO-install/lib:/usr/lib/openfoam/openfoam2406/platforms/linux64GccDPInt32Opt/lib/sys-openmpi:/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib/hpc/gnu7/mpi/openmpi/4.1.6/lib64'

def main():
    # Build T2_7461
    print("Building T2_7461...")
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
    
    # Get grid mesh
    print("Generating grid mesh...")
    bV.get("ru_gridGmsh").makeGrid()
    model = bV.get("ru_gridGmsh").getModel()
    
    # Export the grid mesh
    print("Exporting grid mesh...")
    model.writeMSH("/work/T2_7461_ru_gridGmsh_export.msh")
    
    # Now use Gmsh to process it
    print("Processing with Gmsh...")
    gmsh.initialize()
    gmsh.open("/work/T2_7461_ru_gridGmsh_export.msh")
    
    # Get hub and shroud surfaces
    surfaces = gmsh.model.getEntities(2)
    hub_tags = []
    shroud_tags = []
    
    for dim, tag in surfaces:
        phys_names = gmsh.model.getPhysicalGroupsForEntity(dim, tag)
        for phys in phys_names:
            name = gmsh.model.getPhysicalName(dim, phys)
            if name == "RU_HUB":
                hub_tags.append(tag)
            elif name == "RU_SHROUD":
                shroud_tags.append(tag)
    
    print(f"Hub surfaces: {hub_tags}")
    print(f"Shroud surfaces: {shroud_tags}")
    
    # Set quasi-structured quad algorithm
    gmsh.option.setNumber("Mesh.Algorithm", 11)
    
    # Set element size for hub and shroud
    for tag in hub_tags + shroud_tags:
        gmsh.model.mesh.setSize([(2, tag)], 0.5)
    
    # Remesh
    gmsh.model.mesh.clear()
    gmsh.model.mesh.generate(2)
    
    # Export remeshed surfaces
    print("Exporting remeshed surfaces...")
    gmsh.write("/work/T2_7461_remeshed.msh")
    
    # Extract coarse block structure
    # Get all quads
    _, quad_tags, _ = gmsh.model.mesh.getElementsByType(3)
    print(f"Number of quads: {len(quad_tags)}")
    
    gmsh.finalize()
    print("Done!")

if __name__ == '__main__':
    main()
