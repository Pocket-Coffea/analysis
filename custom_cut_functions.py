import awkward as ak
from pocket_coffea.lib.cut_definition import Cut

def dimuon(events, params, year, sample, **kwargs):

    # Masks for same-flavor (SF) and opposite-sign (OS)
    mask = (events.nJetGood >= 1)

    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

dimuon_presel = Cut(
    name="dilepton",
    params={
        "pt_leading_muon": 25,
        "mll": {'low': 25, 'high': 2000},
    },
    function=dimuon,
)
