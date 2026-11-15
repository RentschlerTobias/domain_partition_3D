"""Build the canada geometry once from state E1_12685 + params.json and write several bounded volumes.
Run inside a copy of /dtOO/test/canadaLight:  python gen_mesh_multi.py ../params.json ru_gridGmsh:../ru.msh gv_gridGmsh:../gv.msh
"""
import json, sys, time
from dtOOPythonSWIG import *

params = json.load(open(sys.argv[1]))
logMe.initLog("build.log")
parser = dtXmlParser.init("machine.xml", "machineSave.xml").reference()
parser.parse()
cV = labeledVectorHandlingConstValue()
parser.createConstValue(cV)
parser.loadStateToConst("E1_12685", cV)
for k, val in params.items():
    cV.get(k).setValue(float(val))
t0 = time.time()
bC = baseContainer(); aF = labeledVectorHandlingAnalyticFunction()
aG = labeledVectorHandlingAnalyticGeometry(); bV = labeledVectorHandlingBoundedVolume()
parser.createAnalyticFunction(bC, cV, aF)
parser.createAnalyticGeometry(bC, cV, aF, aG)
parser.createBoundedVolume(bC, cV, aF, aG, bV)
for spec in sys.argv[2:]:
    label, out = spec.split(":")
    ref = bV.get(label); ref.makeGrid()
    ob = bVOWriteMSH()
    ob.jInit(jsonPrimitive('{"_filename" : "%s", "_saveAll" : true}' % out), None, None, None, None, None, ref)
    ob.postUpdate()
    print(f"[time] {label} {time.time()-t0:.1f}s -> {out}", flush=True)
