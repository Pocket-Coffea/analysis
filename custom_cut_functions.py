import awkward as ak
from pocket_coffea.lib.cut_definition import Cut

def vlt_had_presel(events, params, year, sample, **kwargs):
    
    passed_nPhotonSR = events["nPhotonSR"] > 0
    passed_nPhotonCRB = events["nPhotonCRB"] > 0
    passed_nPhotonCRC = events["nPhotonCRC"] > 0
    passed_nPhotonCRD = events["nPhotonCRD"] > 0
    passed_nPhotonPLJ = events["nPhotonPLJ"] > 0
    passes_photon = passed_nPhotonSR | passed_nPhotonCRB | passed_nPhotonCRC | passed_nPhotonCRD | passed_nPhotonPLJ
    
    passed_nMuonLoose = events["nMuonLoose"] == 0
    passed_nElectronVeto = events["nElectronVeto"] == 0
    passed_nJetGood_NotB = events["nJetGood_NotB"] >= 2
    passed_nBJetGood = events["nBJetGood"] >= 1
    passes_other_objects = passed_nMuonLoose & passed_nElectronVeto & passed_nJetGood_NotB & passed_nBJetGood
    
    mask = passes_photon & passes_other_objects

    return mask

vlt_had_presel = Cut(
    name="vlt_had_presel",
    params={
    },
    function=vlt_had_presel
)

def SR_selection(events, params, year, sample, **kwargs):
    passed_nPhotonSR = events["nPhotonSR"] == 1
    # passed_nPhotonCRB = events["nPhotonCRB"] == 0
    # passed_nPhotonCRC = events["nPhotonCRC"] == 0
    # passed_nPhotonCRD = events["nPhotonCRD"] == 0
    # passed_nPhotonPLJ = events["nPhotonPLJ"] == 0

    return passed_nPhotonSR

SR_cut = Cut(
    name="sr_cut",
    params={
    },
    function=SR_selection
)

def PLJ_selection(events, params, year, sample, **kwargs):
    # passed_nPhotonSR = events["nPhotonSR"] == 0
    passed_nPhotonPLJ = events["nPhotonPLJ"] == 1

    return passed_nPhotonPLJ

PLJ_cut = Cut(
    name="plj_cut",
    params={
    },
    function=PLJ_selection
)

def CRB_selection(events, params, year, sample, **kwargs):
    # passed_nPhotonSR = events["nPhotonSR"] == 0
    passed_nPhotonCRB = events["nPhotonCRB"] == 1
    
    return passed_nPhotonCRB

CRB_cut = Cut(
    name="crb_cut",
    params={
    },
    function=CRB_selection
)

def CRC_selection(events, params, year, sample, **kwargs):
    # passed_nPhotonSR = events["nPhotonSR"] == 1
    passed_nPhotonCRC = events["nPhotonCRC"] == 1

    return passed_nPhotonCRC

CRC_cut = Cut(
    name="crc_cut",
    params={
    },
    function=CRC_selection
)

def CRD_selection(events, params, year, sample, **kwargs):
    # passed_nPhotonSR = events["nPhotonSR"] == 1
    passed_nPhotonCRD = events["nPhotonCRD"] == 1

    return passed_nPhotonCRD

CRD_cut = Cut(
    name="crd_cut",
    params={
    },
    function=CRD_selection
)