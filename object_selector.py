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

    # if region == "CRA":
    #     good_photons = (
    #         passes_eta & passes_pt & passes_pixelseed & passes_transition & mask_lepton_cleaning & 
    #         cutbased_ids_medium["passed_HNP_id"] &
    #         ~cutbased_ids_loose["passSIEIE"] &
    #         cutbased_ids_medium["passed_chIso"]
            
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
        good_photons = (
            (passes_eta & passes_pt & passes_pixelseed & passes_transition & mask_lepton_cleaning & 
            cutbased_ids_medium["passed_hoe"]) &
            ((~cutbased_ids_loose["passSIEIE"] & cutbased_ids_medium["passed_chIso"] & cutbased_ids_medium["passed_neuIso"] & cutbased_ids_medium["passed_phoIso"]) |
            (cutbased_ids_medium["passSIEIE"] & ~cutbased_ids_loose["passed_chIso"] & cutbased_ids_medium["passed_neuIso"] & cutbased_ids_medium["passed_phoIso"]) |
            (cutbased_ids_medium["passSIEIE"] & cutbased_ids_medium["passed_chIso"] & ~cutbased_ids_loose["passed_neuIso"] & cutbased_ids_medium["passed_phoIso"]) |
            (cutbased_ids_medium["passSIEIE"] & cutbased_ids_medium["passed_chIso"] & cutbased_ids_medium["passed_neuIso"] & ~cutbased_ids_loose["passed_phoIso"]) |
            (~cutbased_ids_loose["passSIEIE"] & ~cutbased_ids_loose["passed_chIso"] & cutbased_ids_medium["passed_neuIso"] & cutbased_ids_medium["passed_phoIso"]))
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

    electron_etaSC = events.Electron.eta + events.Electron.deltaEtaSC
    events["Electron"] = ak.with_field(
        events.Electron, electron_etaSC, "etaSC"
    )
    leptons = events[lepton]
    cuts = params.object_preselection[lepton]
    passes_eta = abs(leptons.eta) < cuts["eta"]
    passes_pt = leptons.pt >= cuts[id]["pt"]

    if lepton == "Electron":
        passes_transition = np.invert((abs(leptons.etaSC) >= 1.4442) & (abs(leptons.etaSC) <= 1.5660))
        passes_iso = leptons.pfRelIso03_all < cuts["iso"]
        passes_id = leptons.cutBased >= cuts[id]["cutBased"]

        good_leptons = passes_eta & passes_pt & passes_transition & passes_iso & passes_id

    elif lepton == "Muon":
        passes_iso = leptons.pfRelIso04_all < cuts[id]["iso"]
        passes_id = leptons[cuts[id]["id"]]

        good_leptons = passes_eta & passes_pt & passes_iso & passes_id

    return leptons[good_leptons]

def btagging(Jet, btag, params, veto=False):
    if veto:
        return Jet[(Jet[btag["btagging_algorithm"]] < btag["btagging_WP"][params["wp"]])]
    else:
        return Jet[(Jet[btag["btagging_algorithm"]] > btag["btagging_WP"][params["wp"]]) & (abs(Jet.eta) < params["eta"])]

def calculateNu4vec(lepton, MET):
    # MET components
    MET_pt = MET.pt 
    MET_phi = MET.phi
    MET_px = MET_pt * np.cos(MET_phi)
    MET_py = MET_pt * np.sin(MET_phi) 
    # Lepton components
    lep_m = lepton.mass
    lep_eta = lepton.eta
    lep_pt = lepton.pt
    lep_phi = lepton.phi
    lep_py = lep_pt * np.sin(lep_phi)
    lep_px = lep_pt * np.cos(lep_phi)
    lep_pz = lep_pt * np.sinh(lep_eta)
    lep_E = np.sqrt(lep_px**2 + lep_py**2 + lep_pz**2 + lep_m**2)
    # Constants
    MW = 80.38  # W boson mass in GeV


    #Discriminant
    A = pow(lep_pz,2)-pow(lep_E,2)
    alpha = pow(MW,2)-pow(lep_m,2) + 2 * (lep_px * MET_px + lep_py * MET_py)
    B = alpha * lep_pz
    C = (- pow(lep_E,2) * pow(MET_pt,2) ) + np.divide(pow(alpha,2),4)
    dis = (pow(B,2) - (4*A*C))  # b2-4AC

    condition = dis >= 0
    
    root = np.sqrt(ak.where(condition, dis, ak.zeros_like(dis)))
    root1 = np.divide(- B - root, 2*A)
    root2 = np.divide(- B + root, 2*A)
    pz_nu = ak.where(np.abs(root1) < np.abs(root2), root1, root2)
    E_nu = np.sqrt(MET_pt + pz_nu**2)  

    real_root = ak.where(condition, ak.zeros_like(dis), -B/(2*A)) 
    pz_nu = ak.where(condition, pz_nu, real_root)
    E_nu = np.sqrt(MET_pt + pz_nu**2)  

    
    pt = MET_pt
    phi = MET_phi
    theta = np.arctan2(pt, pz_nu)
    eta = -np.log(np.tan(theta / 2))
    m = np.sqrt(np.maximum(E_nu**2 - (MET_px**2 + MET_py**2 + pz_nu**2), 0))
    nu_p4 = ak.zip({"pt": pt, "eta": eta, "phi": phi, "mass": m},with_name="PtEtaPhiMCandidate")

    return nu_p4