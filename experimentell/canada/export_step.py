"""Build canadaLight state E1_12685 and run the cad_writeStep plugin (all blk boundary surfaces as STEP)."""
import time
from dtOOPythonSWIG import *
logMe.initLog("build_step.log")
parser = dtXmlParser.init("machine.xml", "machineSave.xml").reference()
parser.parse()
cV = labeledVectorHandlingConstValue(); parser.createConstValue(cV); parser.loadStateToConst("E1_12685", cV)
bC = baseContainer(); aF = labeledVectorHandlingAnalyticFunction(); aG = labeledVectorHandlingAnalyticGeometry()
bV = labeledVectorHandlingBoundedVolume(); dC = labeledVectorHandlingDtCase(); dP = labeledVectorHandlingDtPlugin()
parser.createAnalyticFunction(bC, cV, aF); parser.createAnalyticGeometry(bC, cV, aF, aG)
parser.createBoundedVolume(bC, cV, aF, aG, bV); parser.createPlugin(bC, cV, aF, aG, bV, dC, dP)
t0 = time.time(); dP.get("cad_writeStep").apply(); print(f"[step] written in {time.time()-t0:.1f}s", flush=True)
