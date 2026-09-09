"""Build ONLY the block-structured Gmsh mesh of a dtOO case.

The demo `build.py` scripts run the whole chain -- geometry, mesh, decompose,
two simpleFoam runs, reconstruct. This pipeline needs the mesh and nothing
else, and the mesh is one `boundedVolume` of type `map3dTo3dGmsh`: the same
generator that produced `T1_9_ru_gridGmsh.msh` for the tistos case, which is
why its geometric entity ids are the ones `tet_prep_v5.classify` reads.

Run inside the dtOO image, from the case directory:

    docker run --rm -v /root/repos/dtOO:/dtOO -w /dtOO/demo/canada \\
        atismer/dtoo-opensuse:stable bash -lc \\
        'export LD_LIBRARY_PATH=/dtOO-install/lib:$LD_LIBRARY_PATH; \\
         export PYTHONPATH=/dtOO-install/tools:$PYTHONPATH; \\
         python3 dtoo_mesh_only.py machine.xml E1_12685.xml \\
                 ingvrudtout_blkGridMesh canada.msh'

The library path has to be set by hand: the image ships libTKFeat and the rest
of OpenCASCADE in /dtOO-install/lib but puts nothing on LD_LIBRARY_PATH, so a
plain `import dtOOPythonSWIG` fails on a missing shared object.
"""

import sys

from dtOOPythonSWIG import (baseContainer, bVOWriteMSH, dtXmlParser,
                            jsonPrimitive, labeledVectorHandlingAnalyticFunction,
                            labeledVectorHandlingAnalyticGeometry,
                            labeledVectorHandlingBoundedVolume,
                            labeledVectorHandlingConstValue,
                            labeledVectorHandlingDtCase,
                            labeledVectorHandlingDtPlugin, logMe)


def main(machine, state_file, bv_label, out_msh, state=None, plugins=()):
    logMe.initLog("mesh_only.log")
    dtXmlParser.init(machine, state_file)
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
    parser.loadStateToConst(state or state_file.split(".")[0], cV)
    parser.destroyAndCreate(bC, cV, aF, aG, bV, dC, dP)
    # some cases adjust the domain through a plugin and must be rebuilt after
    for p in plugins:
        dP.get(p).apply()
        parser.destroyAndCreate(bC, cV, aF, aG, bV, dC, dP)

    print(f"[dtoo] bounded volumes: {[b.getLabel() for b in bV]}", flush=True)
    ref = bV.get(bv_label)
    print(f"[dtoo] meshing {bv_label} ...", flush=True)
    ref.makeGrid()

    ob = bVOWriteMSH()
    ob.jInit(jsonPrimitive('{"_filename" : "%s", "_saveAll" : true}' % out_msh),
             None, None, None, None, None, ref)
    ob.postUpdate()
    print(f"[dtoo] wrote {out_msh}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) < 5:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
         state=sys.argv[5] if len(sys.argv) > 5 else None)
