from custom_configurator import CustomConfigurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel, get_nPVgood, goldenJson, eventFlags, get_nBtagEq, get_nBtagMin
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.parameters.histograms import *
import workflow
from workflow import TopPartnerBaseProcessor

# Register custom modules in cloudpickle to propagate them to dask workers
import cloudpickle
import custom_cut_functions
import custom_hist_manager
import custom_weight_manager
from custom_weight_manager import common_weights
import object_selector
cloudpickle.register_pickle_by_value(workflow)
cloudpickle.register_pickle_by_value(custom_cut_functions)
cloudpickle.register_pickle_by_value(custom_hist_manager)
cloudpickle.register_pickle_by_value(custom_weight_manager)
cloudpickle.register_pickle_by_value(object_selector)

from custom_cut_functions import *
import os
localdir = os.path.dirname(os.path.abspath(__file__))

# Loading default parameters
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")

parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection_lep.yaml",
                                                  f"{localdir}/params/triggers_lep.yaml",
                                                  f"{localdir}/params/plotting.yaml", 
                                                  f"{localdir}/params/photon_scale_factors.yaml",
                                                  f"{localdir}/params/muon_scale_factors.yaml",
                                                  update=True)



cfg = CustomConfigurator(
    lepton = "Muon",
    parameters = parameters,
    datasets = {
        "jsons": [
            f"{localdir}/datasets/local_mu_fileset.json",
            # f"{localdir}/datasets/DATA_SingleMuon.json",    
            # f"{localdir}/datasets/DYJetsToLL_M-50.json",
            # f"{localdir}/datasets/DYJetsToLL_M10To50.json",
            # f"{localdir}/datasets/TTTo2L2Nu.json",
            # f"{localdir}/datasets/TGJets.json",
            # f"{localdir}/datasets/TTToHadronic.json",
            # f"{localdir}/datasets/TTToSemiLeptonic.json"
        ],
        "filter" : {
            "samples": ["DATA_SingleMuon", "DYJets", "TGJets", "GJets", "TTG", "WJets", "WG", "WWG", "WZG", "ZG", "ST", "TT"], 
            # "TT", "ST"
            "samples_exclude" : ["JetFakePhoton"],
            "year": ['2018']
        }
    },

    workflow = TopPartnerBaseProcessor,

    skim = [#get_nPVgood(1), eventFlags, goldenJson, # basic skims
            #get_nObj_min(1, 18., "Muon"),
            # Asking only SingleMuon triggers since we are only using SingleMuon PD data
            get_HLTsel()], 
    
    preselections = [vlt_presel],
    categories = {
        "SR": [SR_cut]
    },

    weights_classes = common_weights,
    
    weights = {
        "common": {
            "inclusive": ["genWeight","lumi","XS",
                          "pileup", 
                          "custom_sf_mu_id", "custom_sf_mu_iso",
                          "sf_pho_id", "sf_pho_pxseed", "sf_btag"
                          ],
            "bycategory" : {
            }
        },
        "bysample": {
        }
    },

    variations = {
        "weights": {
            "common": {
                "inclusive": [  "pileup",
                                "custom_sf_mu_id", "custom_sf_mu_iso",
                                "sf_pho_id", "sf_pho_pxseed", "sf_btag"
                              ],
                "bycategory" : {
                }
            },
        "bysample": {
        }    
        },
    },

    
   variables = {
       # **lepton_hists(),
       # **jet_hists(),
       "photon_pt" : HistConf( [Axis(coll="PhotonGood", field="pt", bins=[(30+i*20) for i in range(11)], label="$p_T^\gamma$")] ),
       "photon_eta" : HistConf( [Axis(coll="PhotonGood", field="eta", bins=8, start=-2, stop=2, label="$\eta_\gamma$")] ),
       "BJet_pt": HistConf( [Axis(coll="BJetGood", field="pt", bins=[(i*10) for i in range(21)], label="$p_T^b$")] ),
       "BJet_eta": HistConf( [Axis(coll="BJetGood", field="eta", bins=12, start=-3, stop=3, label="$\eta_b$")] ),
       "top_pt": HistConf( [Axis(coll="top", field="pt", bins=[i*10 for i in range(21)], label="$p_T^{top}$")] ),
       "top_mass": HistConf( [Axis(coll="top", field="mass", bins=[(i*10) for i in range(41)], overflow = False, label="$top_M$")] ),
       "VLT_pt": HistConf( [Axis(coll="VLT", field="pt", bins=[i*10 for i in range(21)], overflow=True, label="$p_T^{VLT}$")] ),
       "VLT_mass": HistConf( [Axis(coll="VLT", field="mass", bins=[(100+i*50) for i in range(21)]+[1500,1700,2200], label="$VLT_M$")] ),
       "Muon_pt" : HistConf( [Axis(coll="LeptonGood", field="pt", bins=[(30+i*10) for i in range(11)]+[140,200], label="${p_T}^{\mu}$")] ),
       "WTransverse" : HistConf( [Axis(coll="events", field="W_transMass", bins=30, start=0, stop=150, overflow = False, label="$mW_T$")] )
   }
)

