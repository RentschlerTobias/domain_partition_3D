#------------------------------------------------------------------------------
#  dtOO < design tool Object-Oriented >
#    
#    Copyright (C) 2024 A. Tismer.
#------------------------------------------------------------------------------
#License
#    This file is part of dtOO.
#
#    dtOO is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the LICENSE.txt file in the
#    dtOO root directory for more details.
#
#    You should have received a copy of the License along with dtOO.
#
#------------------------------------------------------------------------------

from dtOOPythonSWIG import *

logMe.initLog("build.log")
dtXmlParser.init("machine.xml", "E1_12685.xml")
parser = dtXmlParser.reference()
parser.parse()
container = dtBundle()
bC = container.cptr_bC()
cV = container.cptr_cV()
aF = container.cptr_aF() 
aG = container.cptr_aG() 
bV = container.cptr_bV() 
dC = container.cptr_dC() 
dP = container.cptr_dP() 

parser.createConstValue(cV)
parser.loadStateToConst("E1_12685", cV)
parser.destroyAndCreate(bC, cV, aF, aG, bV, dC, dP)

dC.get("ingvrudtout_coupled_of").runCurrentState()
