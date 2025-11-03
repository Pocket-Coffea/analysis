from hist import Hist
import awkward as ak


def get_correlation_hist(events, sample):
    
    
    if sample == "WJets":
        photons = events.Photon[ak.fill_none(events.Photon.matched_gen.pdgId == None, True)]
    else:
        photon_prompt = events.LHEPart[events.LHEPart.pdgId == 22]
        photons = events.Photon[ak.any(events.Photon.metric_table(photon_prompt) < 0.3, axis=2)]

    hist_2d = Hist.new.Reg(
        100, .006, 0.016, name="sieie", label="Sigma Ieta Ieta"  # Adjust binning and range for sieie
    ).Reg(
        50, 0, 10, name="chiso", label="Charged Hadron Isolation"  # Adjust binning and range for chiso
    ).Double()

    hist_2d.fill(sieie=ak.to_numpy(ak.flatten(photons.sieie), allow_missing=False), chiso=ak.to_numpy(ak.flatten(photons.chIso), allow_missing=False))
    
    return hist_2d