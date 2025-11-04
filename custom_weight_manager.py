from coffea.analysis_tools import Weights
from custom_cut_functions import *
from pocket_coffea.lib.categorization import StandardSelection
from pocket_coffea.lib.weights import WeightData ,WeightWrapper, WeightLambda, WeightDataMultiVariation
import numpy as np
import awkward as ak

from custom_scale_factors import *

from pocket_coffea.lib.weights.common.common import (
    genWeight,
    signOfGenWeight,
    lumi,
    XS,
    pileup,
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
        self.regions = StandardSelection({"PLJ": [PLJ_cut], "b0_PLJ": [b0_PLJ_cut]})
        self.regions.prepare(
            events=self.events,
            processor_params={},
        )
        self.EF = {
            "b1": {
                "Electron": {
                    "2018":{
                        '[30, 40]': 0.2792079273132033,
                        '[40, 50]': 0.2597070300006397,
                        '[50, 60]': 0.17540994578676128,
                        '[60, 80]': 0.18533123028391169,
                        '[80, 100]': 0.15953716690042075,
                        '[100, np.inf]': 0.138282153411281
                    }
                },
                "Muon":{
                     "2018":{
                         '[30, 40]': 0.13074266214115549,
                         '[40, 50]': 0.11580594679186229,
                         '[50, 60]': 0.09128600266114875,
                         '[60, 80]': 0.06784048242120833,
                         '[80, 100]': 0.04153560194452388,
                         '[100, np.inf]': 0.05357142857142857
                     }
                 },
                "Lepton":{
                     "2018":{
                         '[60, 80]': 0.079798,
                         '[80, 100]': 0.0706147,
                         "[100, 150]": 0.06082,
                         '[150, np.inf]': 0.06299762
                     },
                }
            },
            "b0": {
                "Electron": {
                    "2018":{
                        '[30, 40]': 0.2792079273132033,
                        '[40, 50]': 0.2597070300006397,
                        '[50, 60]': 0.17540994578676128,
                        '[60, 80]': 0.18533123028391169,
                        '[80, 100]': 0.15953716690042075,
                        '[100, np.inf]': 0.138282153411281
                    }
                },
                "Muon":{
                     "2018":{
                         '[30, 40]': 0.13074266214115549,
                         '[40, 50]': 0.11580594679186229,
                         '[50, 60]': 0.09128600266114875,
                         '[60, 80]': 0.06784048242120833,
                         '[80, 100]': 0.04153560194452388,
                         '[100, np.inf]': 0.05357142857142857
                     }
                 },
                "Lepton":{
                     "2018":{
                         '[60, 80]': 0.1041118,
                         '[80, 100]': 0.0872647,
                         "[100, 150]": 0.0738956,
                         '[150, np.inf]': 0.0573968
                     },
                },
            },
        }

    def compute_EF(self, year, cat):
        mask = self.regions.get_mask(cat)
        masked_events = self.events[mask]
        photon = ak.flatten(masked_events[f"Photon{cat[-3:]}"])

        if "b0" in cat:
            ef = self.EF["b0"]
        else:
            ef = self.EF["b1"]

        flavor = "Lepton" # self.events.flavor[0] for seprate channels
        extrapolation_factor = ak.where(
            (photon.pt >= 60) & (photon.pt < 80),
            ef[flavor][year]['[60, 80]'],
            ak.where(
                (photon.pt >= 80) & (photon.pt < 100),
                ef[flavor][year]['[80, 100]'],
                ak.where(
                    (photon.pt >= 100) & (photon.pt < 150),
                    ef[flavor][year]['[100, 150]'],
                    ak.where(
                        (photon.pt >= 150),
                        ef[flavor][year]['[150, np.inf]'],
                        1
                    )
                )
            )
        )
        weight = WeightData(
            name = "EF",
            nominal = extrapolation_factor
        )

        self.EF_weight = Weights(len(self.events[self.regions.get_mask(cat)]))
        self.EF_weight.add(weight.name, weight.nominal)

        return self.EF_weight.weight()

SF_ele_reco = WeightLambda.wrap_func(
    name="custom_sf_ele_reco",
    function=lambda params, metadata, events, size, shape_variations:
        sf_ele_reco(params, events, metadata["year"]),
    has_variations=True
    )

SF_ele_id = WeightLambda.wrap_func(
    name="custom_sf_ele_id",
    function=lambda params, metadata, events, size, shape_variations:
        sf_ele_id(params, events, metadata["year"]),
    has_variations=True
    )

SF_mu_id = WeightLambda.wrap_func(
    name="custom_sf_mu_id",
    function=lambda params, metadata, events, size, shape_variations:
        sf_mu(params, events, metadata["year"], 'id'),
    has_variations=True
    )

SF_mu_iso = WeightLambda.wrap_func(
    name="custom_sf_mu_iso",
    function=lambda params, metadata, events, size, shape_variations:
        sf_mu(params, events, metadata["year"], 'iso'),
    has_variations=True
    )

SF_pho_pxseed = WeightLambda.wrap_func(
    name="custom_sf_pho_pxseed",
    function=lambda params, metadata, events, size, shape_variations:
        sf_pho_pxseed(params, events, metadata["year"]),
    has_variations=True
    )

SF_pho_id = WeightLambda.wrap_func(
    name="custom_sf_pho_id",
    function=lambda params, metadata, events, size, shape_variations:
        sf_pho_id(params, events, metadata["year"]),
    has_variations=True
    )


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
    SF_pho_pxseed,
    SF_pho_id,
    SF_btag,
    SF_btag_calib,
    SF_ctag,
    SF_ctag_calib,
    SF_jet_puId,
    # SF_PSWeight_isr,
    # SF_PSWeight_fsr
]

