"""Build one bounded volume of the canadaLight case and write it as .msh.

Like datagen/gen_mesh.py (tistos), but without dtCase (the OpenFOAM builder is
missing in this dtOO build) and without ru_adjustDomain (canada has no such plugin).
Run inside a copy of /dtOO/test/canadaLight:
    python gen_mesh_canada.py <state> <bv_label> <out.msh> [params.json]
"""
import json, os, sys, time
from dtOOPythonSWIG import *

state, bv_label, out = sys.argv[1:4]
logMe.initLog("build.log")
parser = dtXmlParser.init("machine.xml", "machineSave.xml").reference()
parser.parse()
cV = labeledVectorHandlingConstValue()
parser.createConstValue(cV)
parser.loadStateToConst(state, cV)
if len(sys.argv) > 4:
    for k, v in json.load(open(sys.argv[4])).items():
        cV.get(k).setValue(float(v))

t0 = time.time()
bC = baseContainer(); aF = labeledVectorHandlingAnalyticFunction()
aG = labeledVectorHandlingAnalyticGeometry(); bV = labeledVectorHandlingBoundedVolume()
parser.createAnalyticFunction(bC, cV, aF)
parser.createAnalyticGeometry(bC, cV, aF, aG)
parser.createBoundedVolume(bC, cV, aF, aG, bV)
print(f"[dtoo] bounded volumes: {[b.getLabel() for b in bV]}", flush=True)
ref = bV.get(bv_label); ref.makeGrid()
ob = bVOWriteMSH()
ob.jInit(jsonPrimitive('{"_filename" : "%s", "_saveAll" : %s}' % (out, os.environ.get("SAVE_ALL", "true"))), None, None, None, None, None, ref)
ob.postUpdate()
print(f"[time] dtoo {time.time()-t0:.1f}s -> {out}", flush=True)
