import numpy as np
import awkward as ak

def photon_selection(events, photon, params, region, leptons_collection=""):

    photons = events["Photon"]
    cuts = params.object_preselection[photon]
    # Requirements on pT and eta
    passes_eta = abs(photons.eta) < cuts["eta"]
    passes_transition = np.invert(( abs(photons.eta) >= 1.4442) & (abs(photons.eta) <= 1.5660))
    passes_pt = photons.pt > cuts["pt"]
    passes_pixelseed = ~photons.pixelSeed

    if leptons_collection != "":
        dR_photons_lep = photons.metric_table(events[leptons_collection])
        mask_lepton_cleaning = ak.prod(dR_photons_lep > cuts["dr_lepton"], axis=2) == 1
    else:
        mask_lepton_cleaning = True
        
    bitMap = photons.vidNestedWPBitmap
    cutbased_ids_medium = parse_photon_vid_cuts(bitMap, 2)
    cutbased_ids_loose = parse_photon_vid_cuts(bitMap, 1)
    
    if region == "SR":
        good_photons = (
            passes_eta & passes_pt & passes_pixelseed & passes_transition & mask_lepton_cleaning & 
            cutbased_ids_medium["passed_HNP_id"] &
            cutbased_ids_medium["passSIEIE"] &
            cutbased_ids_medium["passed_chIso"]
        )
    if region == "CRB":
        good_photons = (
            passes_eta & passes_pt & passes_pixelseed & passes_transition & mask_lepton_cleaning & 
            cutbased_ids_medium["passed_HNP_id"] &
            ~cutbased_ids_loose["passSIEIE"] &
            cutbased_ids_medium["passed_chIso"]
        )
    if region == "CRC":
        good_photons = (
            passes_eta & passes_pt & passes_pixelseed & passes_transition & mask_lepton_cleaning & 
            cutbased_ids_medium["passed_HNP_id"] &
            cutbased_ids_medium["passSIEIE"] &
            ~cutbased_ids_loose["passed_chIso"]
        )
    if region == "CRD":
        good_photons = (
            passes_eta & passes_pt & passes_pixelseed & passes_transition & mask_lepton_cleaning & 
            cutbased_ids_medium["passed_HNP_id"] &
            ~cutbased_ids_loose["passSIEIE"] &
            ~cutbased_ids_loose["passed_chIso"]
        )
    if region == "PLJ":
        passes
        good_photons = (
            (passes_eta & passes_pt & passes_pixelseed & passes_transition & mask_lepton_cleaning & 
            cutbased_ids_medium["passed_hoe"]) &
            ((~cutbased_ids_loose["passSIEIE"] & cutbased_ids_medium["passed_chIso"] & cutbased_ids_medium["passed_neuIso"] & cutbased_ids_medium["passed_phoIso"]) |
            (cutbased_ids_medium["passSIEIE"] & ~cutbased_ids_loose["passed_chIso"] & cutbased_ids_medium["passed_neuIso"] & cutbased_ids_medium["passed_phoIso"]) |
            (cutbased_ids_medium["passSIEIE"] & cutbased_ids_medium["passed_chIso"] & ~cutbased_ids_loose["passed_neuIso"] & cutbased_ids_medium["passed_phoIso"]) |
            (cutbased_ids_medium["passSIEIE"] & cutbased_ids_medium["passed_chIso"] & cutbased_ids_medium["passed_neuIso"] & ~cutbased_ids_loose["passed_phoIso"]))
        )

    return photons[good_photons]

def parse_photon_vid_cuts(bitMap, cutLevel):
    cutbased_ids = {}
    # Create a list of cuts as per the C++ function logic
    passHoverE = ((bitMap >> 4) & 3) >= cutLevel
    passSIEIE = ((bitMap >> 6) & 3) >= cutLevel
    passChIso = ((bitMap >> 8) & 3) >= cutLevel
    passNeuIso = ((bitMap >> 10) & 3) >= cutLevel
    passPhoIso = ((bitMap >> 12) & 3) >= cutLevel

    cutbased_ids["passID"] = passHoverE & passSIEIE & passChIso & passNeuIso & passPhoIso
    cutbased_ids["passed_hoe"] = passHoverE
    cutbased_ids["passed_neuIso"] = passNeuIso
    cutbased_ids["passed_phoIso"] = passPhoIso
    cutbased_ids["passed_HNP_id"] = passHoverE & passNeuIso & passPhoIso
    cutbased_ids["passed_chIso"] = passChIso
    cutbased_ids["passSIEIE"] = passSIEIE

    return cutbased_ids

def lepton_selection(events, lepton, params, id):

    leptons = events[lepton]
    cuts = params.object_preselection[lepton]
    passes_eta = abs(leptons.eta) < cuts["eta"]
    passes_pt = leptons.pt >= cuts[id]["pt"]

    if lepton == "Electron":
        passes_transition = np.invert((abs(leptons.eta) >= 1.4442) & (abs(leptons.eta) <= 1.5660))
        passes_iso = leptons.pfRelIso03_all < cuts["iso"]
        passes_id = leptons.cutBased >= cuts[id]["cutBased"]

        good_leptons = passes_eta & passes_pt & passes_transition & passes_iso & passes_id

    elif lepton == "Muon":
        passes_iso = leptons.pfRelIso04_all < cuts[id]["iso"]
        passes_id = leptons[cuts[id]["id"]]

        good_leptons = passes_eta & passes_pt & passes_iso & passes_id

    return leptons[good_leptons]

def btagging(Jet, btag, params, veto=False):
    cuts = params.object_preselection["Jet"]["btag"]
    if veto:
        return Jet[(Jet[btag["btagging_algorithm"]] < btag["btagging_WP"][cuts["wp"]]) & (abs(Jet.eta < cuts["eta"]))]
    else:
        return Jet[(Jet[btag["btagging_algorithm"]] > btag["btagging_WP"][cuts["wp"]]) & (abs(Jet.eta < cuts["eta"]))]