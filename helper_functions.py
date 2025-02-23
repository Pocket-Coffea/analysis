from hist import Hist
import awkward as ak


def get_correlation_hist(events, sample):
    
    mask = ak.fill_none(events.Photon.matched_gen.pdgId == 22, False)
    if sample == "WJets":
        photons = events.Photon[~mask]
    else:
        photons = events.Photon[mask]

    hist_2d = Hist.new.Reg(
        100, .006, 0.016, name="sieie", label="Sigma Ieta Ieta"  # Adjust binning and range for sieie
    ).Reg(
        50, 0, 10, name="chiso", label="Charged Hadron Isolation"  # Adjust binning and range for chiso
    ).Double()

    hist_2d.fill(sieie=ak.to_numpy(ak.flatten(photons.sieie), allow_missing=False), chiso=ak.to_numpy(ak.flatten(photons.chIso), allow_missing=False))
    
    return hist_2d