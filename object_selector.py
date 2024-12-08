import numpy as np

# def photon_selection(events, photon, params, region):

#     photons = events["Photons"]
#     cuts = params.object_preselection[photon]
#     # Requirements on pT and eta
#     passes_eta = abs(photons.eta) < cuts["eta"]
#     passes_transition = np.invert(( abs(photons.eta) >= 1.4442) & (abs(photons.eta) <= 1.5660))
#     passes_pt = photons.pt > cuts["pt"]
#     passes_pixelseed = not photons.pixelSeed
#     if region == "SR":
#         if abs(photons.eta) < 1.4442:
#             passes_hoe = photons.hoe < cuts["barrel"]["SR"]["hoe"]
#             passes_sieie = photons.sieie < cuts["barrel"]["SR"]["sieie"]
#             passes_chiso = photons.pfRelIso03_chg < cuts["barrel"]["SR"]["chiso"]
#         if abs(photons.eta) > 1.5660 and abs(photons.eta) < cuts["eta"]:
#             passes_hoe = photons.hoe < cuts["endcap"]["SR"]["hoe"]
#             passes_sieie = photons.sieie < cuts["endcap"]["SR"]["sieie"]
#             passes_chiso = photons.pfRelIso03_chg < cuts["endcap"]["SR"]["chiso"]
           
        
    

#     good_photons = passes_eta & passes_pt & passes_pixelseed & passes_id

#     return photons[good_photons]

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

    elif lepton_flavour == "Muon":
        passes_iso = leptons.pfRelIso04_all < cuts[id]["iso"]
        passes_id = leptons[cuts[id]["id"]]

        good_leptons = passes_eta & passes_pt & passes_iso & passes_id

    return leptons[good_leptons]


def jet_selection(events, jet_type, params, year, leptons_collection="", jet_tagger=""):

    jets = events[jet_type]
    cuts = params.object_preselection[jet_type]
    # Only jets that are more distant than dr to ALL leptons are tagged as good jets
    # Mask for  jets not passing the preselection
    mask_presel = (
        (jets.pt > cuts["pt"])
        & (np.abs(jets.eta) < cuts["eta"])
        & (jets.jetId > cuts["tightLeptonVetoId"])
    )
    # Lepton cleaning
    if leptons_collection != "":
        dR_jets_lep = jets.metric_table(events[leptons_collection])
        mask_lepton_cleaning = ak.prod(dR_jets_lep > cuts["dr_lepton"], axis=2) == 1
    else:
        mask_lepton_cleaning = True

    if jet_type == "Jet":
        # Selection on PUid. Only available in Run2 UL, thus we need to determine which sample we run over;
        if year in ['2016_PreVFP', '2016_PostVFP','2017','2018']:
            mask_jetpuid = (jets.puId >= params.jet_scale_factors.jet_puId[year]["working_point"][cuts["puId"]["wp"]]) | (
                jets.pt >= cuts["puId"]["maxpt"]
            )
        else:
            mask_jetpuid = True
  
        mask_good_jets = mask_presel & mask_lepton_cleaning & mask_jetpuid