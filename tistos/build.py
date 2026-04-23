

def build():

    import logging
    import os
    os.environ["OSLO_LOCK_PATH"] = "/tmp/oslo_locks/"
    import sys
    sys.path.insert(0, '/mnt/opt.net/src/dtOO/ThirdParty~opensuse-leap~15.6/lib64')
    import gmsh

    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    os.chdir(BASE_DIR)
    logging.basicConfig(
        format='%(asctime)s %(levelname)s : %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.INFO
    )
    from dtOOPythonSWIG import (
        logMe,
        dtXmlParser,
        baseContainer,
        labeledVectorHandlingConstValue,
        labeledVectorHandlingAnalyticFunction,
        labeledVectorHandlingAnalyticGeometry,
        labeledVectorHandlingBoundedVolume,
        labeledVectorHandlingDtCase,
        labeledVectorHandlingDtPlugin,
        lVHOstateHandler,
    )
    from pyDtOO import (
        dtScalarDeveloping,
        dtForceDeveloping,
        dtDeveloping
    )
    from pyDtOO import dtClusteredSingletonState as stateCounter
    import foamlib
    import numpy as np
    import sys
    import os
    clone = "git clone https://github.com/ihs-ustutt/axial_turbine_database.git"
    if not os.path.isdir("./axial_turbine_database"):
            logging.info("Clone repository.")
            ret = os.system(clone)
    stateCounter.PREFIX = 'T1'
    stateCounter.CASE = 'tistos_ru_of'
    stateCounter.DATADIR = './axial_turbine_database/runData'
    stateCounter.ADDDATA = [
        'P',
        'dH',
        'eta',
        'VCav',
        'history',
        'islandID',
    ]
    stateCounter.ADDDATADEF = [
        {"tl": 0, "n": 0, "vl": 0},  # P
        {"tl": 0, "n": 0, "vl": 0},  # dH
        {"tl": 0, "n": 0, "vl": 0},  # eta
        {"tl": 0, "n": 0, "vl": 0},  # VCav
        {},  # history
        -1,  # islandID
    ]
    cVArr = [
        {'label': 'cV_ru_alpha_1_ex_0.0', 'min': -0.155, 'max': 0.025},
        {'label': 'cV_ru_alpha_1_ex_0.5', 'min': -0.19, 'max': -0.01},
        {'label': 'cV_ru_alpha_1_ex_1.0', 'min': -0.19, 'max': -0.01},
        {'label': 'cV_ru_alpha_2_ex_0.0', 'min': -0.08, 'max': 0.1},
        {'label': 'cV_ru_alpha_2_ex_0.5', 'min': -0.08, 'max': 0.1},
        {'label': 'cV_ru_alpha_2_ex_1.0', 'min': -0.08, 'max': 0.07},
        {'label': 'cV_ru_offsetM_ex_0.0', 'min': 1.0, 'max': 1.5},
        {'label': 'cV_ru_offsetM_ex_0.5', 'min': 1.0, 'max': 1.5},
        {'label': 'cV_ru_offsetM_ex_1.0', 'min': 1.0, 'max': 1.5},
        {'label': 'cV_ru_ratio_0.0', 'min': 0.4, 'max': 0.6},
        {'label': 'cV_ru_ratio_0.5', 'min': 0.4, 'max': 0.6},
        {'label': 'cV_ru_ratio_1.0', 'min': 0.4, 'max': 0.6},
        {'label': 'cV_ru_offsetPhiR_ex_0.0', 'min': -0.15, 'max': 0.15},
        {'label': 'cV_ru_offsetPhiR_ex_0.5', 'min': -0.15, 'max': 0.15},
        {'label': 'cV_ru_offsetPhiR_ex_1.0', 'min': -0.15, 'max': 0.15},
        {'label': 'cV_ru_bladeLength_0.0', 'min': 0.4, 'max': 0.8},
        {'label': 'cV_ru_bladeLength_0.5', 'min': 0.6, 'max': 1.0},
        {'label': 'cV_ru_bladeLength_1.0', 'min': 0.8, 'max': 1.3},
        {'label': 'cV_ru_t_le_a_0', 'min': 0.005, 'max': 0.06},
        {'label': 'cV_ru_t_le_a_0.5', 'min': 0.005, 'max': 0.06},
        {'label': 'cV_ru_t_le_a_1', 'min': 0.005, 'max': 0.06},
        {'label': 'cV_ru_t_mid_a_0', 'min': 0.005, 'max': 0.06},
        {'label': 'cV_ru_t_mid_a_0.5', 'min': 0.005, 'max': 0.06},
        {'label': 'cV_ru_t_mid_a_1', 'min': 0.005, 'max': 0.06},
        {'label': 'cV_ru_t_te_a_0', 'min': 0.005, 'max': 0.06},
        {'label': 'cV_ru_t_te_a_0.5', 'min': 0.005, 'max': 0.06},
        {'label': 'cV_ru_t_te_a_1', 'min': 0.005, 'max': 0.06},
        {'label': 'cV_ru_u_mid_a_0', 'min': 0.4, 'max': 0.6},
        {'label': 'cV_ru_u_mid_a_0.5', 'min': 0.4, 'max': 0.6},
        {'label': 'cV_ru_u_mid_a_1', 'min': 0.4, 'max': 0.6},
    ]
    sc = stateCounter(21260)
    logMe.initLog('build-'+sc.state()+'.log')
    parser = dtXmlParser.init("machine.xml", "machineSave.xml").reference()
    parser.parse()
    bC = baseContainer()
    cV = labeledVectorHandlingConstValue()
    aF = labeledVectorHandlingAnalyticFunction()
    aG = labeledVectorHandlingAnalyticGeometry()
    bV = labeledVectorHandlingBoundedVolume()
    dC = labeledVectorHandlingDtCase()
    dP = labeledVectorHandlingDtPlugin()
    parser.createConstValue(cV)
    parser.loadStateToConst("templateState", cV)
    cc = 0
    for anObj in sc.objective():
        cV[cVArr[cc]['label']].setValue(anObj)
        cc = cc + 1
    parser.destroyAndCreate(bC, cV, aF, aG, bV, dC, dP)
    dP.get('ru_adjustDomain').apply()
    parser.destroyAndCreate(bC, cV, aF, aG, bV, dC, dP)
    lVHOstateHandler().makeState(sc.state())

    return bC, cV, aF, aG, bV, dC, dP 


