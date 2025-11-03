import awkward as ak
from pocket_coffea.lib.cut_definition import Cut

def vlt_presel(events, params, **kwargs):
    
    passed_nPhotonGood = events["nPhotonGood"] > 0
    passed_nPhotonCRB = events["nPhotonCRB"] > 0
    passed_nPhotonCRC = events["nPhotonCRC"] > 0
    passed_nPhotonCRD = events["nPhotonCRD"] > 0
    passed_nPhotonPLJ = events["nPhotonPLJ"] > 0
    passed_nPhoton = passed_nPhotonGood | passed_nPhotonCRB | passed_nPhotonCRC | passed_nPhotonCRD | passed_nPhotonPLJ
    
    passed_nLeptonGood = events["nLeptonGood"] == 1
    if events.flavor[0] == "Electron":
        passed_nMuonLoose = events["nMuonLoose"] == 0
        passed_nElectronVeto = events["nElectronVeto"] == events["nLeptonGood"]
    elif events.flavor[0] == "Muon":
        passed_nMuonLoose = events["nMuonLoose"] == events["nLeptonGood"]
        passed_nElectronVeto = events["nElectronVeto"] == 0
    else:
        raise Excepion("No Leton")
    passed_met_pt = events.MET.pt > 30
    
    mask = passed_nPhoton & passed_nLeptonGood & passed_nMuonLoose & passed_nElectronVeto & passed_met_pt

    return mask

vlt_presel = Cut(
    name="vlt_presel",
    params={
    },
    function=vlt_presel
)

def SR_selection(events, params, **kwargs):
    passed_nPhotonGood = events["nPhotonGood"] == 1
    # passed_nPhotonCRB = events["nPhotonCRB"] == 0
    # passed_nPhotonCRC = events["nPhotonCRC"] == 0
    # passed_nPhotonCRD = events["nPhotonCRD"] == 0
    # passed_nPhotonPLJ = events["nPhotonPLJ"] == 0
    passed_nBJetGood = events["nBJetGood"] == params["nb"]
    # passed_nJetGood = events["nJetGood"] == 1

    return passed_nPhotonGood & passed_nBJetGood

SR_cut = Cut(
    name="sr",
    params={
        "nb": 1
    },
    function=SR_selection
)
b0_SR_cut = Cut(
    name="b0_sr",
    params={
        "nb": 0
    },
    function=SR_selection
)

def PLJ_selection(events, params, **kwargs):
    passed_nPhotonGood = events["nPhotonGood"] == 0
    passed_nPhotonPLJ = events["nPhotonPLJ"] == 1
    passed_nBJetGood = events["nBJetGoodPLJ"] == params["nb"]
    # passed_nJetGood = events["nJetGoodPLJ"] == 1

    return passed_nPhotonPLJ & passed_nPhotonGood & passed_nBJetGood

PLJ_cut = Cut(
    name="plj",
    params={
        "nb": 1
    },
    function=PLJ_selection
)
b0_PLJ_cut = Cut(
    name="b0_plj",
    params={
        "nb": 0
    },
    function=PLJ_selection
)

def CRB_selection(events, params, **kwargs):
    passed_nPhotonGood = events["nPhotonGood"] == 0
    passed_nPhotonCRB = events["nPhotonCRB"] == 1
    passed_nPhotonCRC = events["nPhotonCRC"] == 0
    passed_nPhotonCRD = events["nPhotonCRD"] == 0
    passed_nBJetGood = events["nBJetGoodCRB"] == params["nb"]
    # passed_nJetGood = events["nJetGoodCRB"] == 1
    
    return passed_nPhotonGood & passed_nPhotonCRB & passed_nPhotonCRC & passed_nPhotonCRD & passed_nBJetGood

CRB_cut = Cut(
    name="crb",
    params={
        "nb": 1
    },
    function=CRB_selection
)
b0_CRB_cut = Cut(
    name="b0_crb",
    params={
        "nb": 0
    },
    function=CRB_selection
)

def CRC_selection(events, params, **kwargs):
    passed_nPhotonGood = events["nPhotonGood"] == 0
    passed_nPhotonCRB = events["nPhotonCRB"] == 0
    passed_nPhotonCRC = events["nPhotonCRC"] == 1
    passed_nPhotonCRD = events["nPhotonCRD"] == 0
    passed_nBJetGood = events["nBJetGoodCRC"] == params["nb"]
    # passed_nJetGood = events["nJetGoodCRC"] == 1

    return passed_nPhotonGood & passed_nPhotonCRB & passed_nPhotonCRC & passed_nPhotonCRD & passed_nBJetGood

CRC_cut = Cut(
    name="crc",
    params={
        "nb": 1
    },
    function=CRC_selection
)
b0_CRC_cut = Cut(
    name="b0_crc",
    params={
        "nb": 0
    },
    function=CRC_selection
)

def CRD_selection(events, params, **kwargs):
    passed_nPhotonGood = events["nPhotonGood"] == 0
    passed_nPhotonCRB = events["nPhotonCRB"] == 0
    passed_nPhotonCRC = events["nPhotonCRC"] == 0
    passed_nPhotonCRD = events["nPhotonCRD"] == 1
    passed_nBJetGood = events["nBJetGoodCRD"] == params["nb"]
    # passed_nJetGood = events["nJetGoodCRD"] == 1

    return passed_nPhotonGood & passed_nPhotonCRB & passed_nPhotonCRC & passed_nPhotonCRD & passed_nBJetGood

CRD_cut = Cut(
    name="crd",
    params={
        "nb": 1
    },
    function=CRD_selection
)
b0_CRD_cut = Cut(
    name="b0_crd",
    params={
        "nb": 0
    },
    function=CRD_selection
)