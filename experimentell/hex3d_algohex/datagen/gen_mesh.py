"""Build the tistos geometry from params.json and write ru_gridGmsh as mesh.msh.

Same order as the dataset export (eigenfrequencies dtoo/export.py): load
templateState, set the parameters, build, ru_adjustDomain, build again, mesh.
Without dtCase (the OpenFOAM builder is missing in this dtOO build).
Run inside a copy of /dtOO/demo/tistos:  python gen_mesh.py ../params.json
"""
import json, sys, time
from dtOOPythonSWIG import *

params = json.load(open(sys.argv[1]))
logMe.initLog("build.log")
parser = dtXmlParser.init("machine.xml", "machineSave.xml").reference()
parser.parse()
cV = labeledVectorHandlingConstValue()
parser.createConstValue(cV)
parser.loadStateToConst("templateState", cV)
for k, v in params.items():
    cV.get(k).setValue(float(v))


def build():
    bC = baseContainer(); aF = labeledVectorHandlingAnalyticFunction()
    aG = labeledVectorHandlingAnalyticGeometry(); bV = labeledVectorHandlingBoundedVolume()
    dC = labeledVectorHandlingDtCase(); dP = labeledVectorHandlingDtPlugin()
    parser.createAnalyticFunction(bC, cV, aF)
    parser.createAnalyticGeometry(bC, cV, aF, aG)
    parser.createBoundedVolume(bC, cV, aF, aG, bV)
    parser.createPlugin(bC, cV, aF, aG, bV, dC, dP)
    return bC, aF, aG, bV, dC, dP


t0 = time.time()
c = build()
c[5].get("ru_adjustDomain").apply()
bC, aF, aG, bV, dC, dP = build()
ref = bV.get("ru_gridGmsh"); ref.makeGrid()
ob = bVOWriteMSH()
ob.jInit(jsonPrimitive('{"_filename" : "mesh.msh", "_saveAll" : true}'), None, None, None, None, None, ref)
ob.postUpdate()
print(f"[time] dtoo {time.time()-t0:.1f}s", flush=True)
