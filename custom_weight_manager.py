from coffea.analysis_tools import Weights
from custom_cut_functions import *
from pocket_coffea.lib.categorization import StandardSelection
from pocket_coffea.lib.weights import WeightData #,WeightWrapper, WeightLambda, WeightDataMultiVariation
import numpy as np
import awkward as ak

from pocket_coffea.lib.weights.common.common import (
    genWeight,
    signOfGenWeight,
    lumi,
    XS,
    pileup,
    SF_ele_reco,
    SF_ele_id,
    SF_mu_id,
    SF_mu_iso,
    SF_mu_trigger,
    SF_btag,
    SF_btag_calib,
    SF_ctag,
    SF_ctag_calib,
    SF_jet_puId,
    # SF_PSWeight_isr,
    # SF_PSWeight_fsr
)


class ExtrapolationFactor:
    def __init__(self, events):
        self.events = events
        self.regions = StandardSelection({"PLJ": [PLJ_cut], "CRB": [CRB_cut], "CRC": [CRC_cut], "CRD": [CRD_cut]})
        self.regions.prepare(
            events=self.events,
            processor_params={},
        )
        self.EF = {
                 "Electron": {
                     "2018":{
                         '[30, 40]': 0.1721540584018723,
                         '[40, 50]': 0.15993839045051983,
                         '[50, 70]': 0.12443778110944527,
                         '[70, 100]': 0.09195402298850576,
                         '[100, 140]': 0.06654665094046806,
                         '[140, 200]': 0.06585365853658537,
                         '[200, np.inf]': 0.041025641025641026
                     }
                 }
        }
        self.EF_weight = Weights(len(self.events[self.regions.get_mask("PLJ")]))

    def get_EF_for_pt(self, pt_down, pt_up):
        
        nevents = {}
        for region, mask in self.regions.get_masks():
            masked_events = self.events[mask]
            photon = ak.flatten(getattr(masked_events, "Photon{}".format(region)))
            selected_events = masked_events[(photon.pt >= pt_down) & (photon.pt < pt_up)]
            nevents[region] = len(selected_events)
        return ((nevents["CRB"] * nevents["CRC"])/nevents["CRD"])/nevents["PLJ"]

    def compute_EF(self, year):
        mask = self.regions.get_mask("PLJ")
        masked_events = self.events[mask]
        photon = ak.flatten(masked_events["PhotonPLJ"])
        
        extrapolation_factor = ak.where(
            (photon.pt >= 30) & (photon.pt < 40),
            self.EF[self.events.flavor[0]][year]['[30, 40]'],
            ak.where(
                (photon.pt >= 40) & (photon.pt < 50),
                self.EF[self.events.flavor[0]][year]['[40, 50]'],
                ak.where(
                    (photon.pt >= 50) & (photon.pt < 70),
                    self.EF[self.events.flavor[0]][year]['[50, 70]'],
                    ak.where(
                        (photon.pt >= 70) & (photon.pt < 100),
                        self.EF[self.events.flavor[0]][year]['[70, 100]'],
                        ak.where(
                            (photon.pt >= 100) & (photon.pt < 140),
                            self.EF[self.events.flavor[0]][year]['[100, 140]'],
                            ak.where(
                                (photon.pt >= 140) & (photon.pt < 200),
                                self.EF[self.events.flavor[0]][year]['[140, 200]'],
                                ak.where(
                                    (photon.pt >= 200) & (photon.pt < 300),
                                    self.EF[self.events.flavor[0]][year]['[200, np.inf]'],
                                    1
                                )
                            )
                        )
                    )
                )
            )
        )
        weight = WeightData(
            name = "EF",
            nominal = extrapolation_factor
        )

        self.EF_weight.add(weight.name, weight.nominal)

        return self.EF_weight.weight()


# SF_ele_trigger = WeightLambda.wrap_func(
#     name="sf_ele_trigger",
#     function=lambda params, metadata, events, size, shape_variations:
#         sf_ele(params, events, metadata["year"], 'trigger'),
#     has_variations=True
#     )


common_weights = [
    genWeight,
    signOfGenWeight,
    lumi,
    XS,
    pileup,
    SF_ele_reco,
    SF_ele_id,
    SF_mu_id,
    SF_mu_iso,
    SF_mu_trigger,
    SF_btag,
    SF_btag_calib,
    SF_ctag,
    SF_ctag_calib,
    SF_jet_puId,
    # SF_PSWeight_isr,
    # SF_PSWeight_fsr
]

