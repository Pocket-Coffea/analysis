from pocket_coffea.utils.configurator import Configurator
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
                                                  update=True)



cfg = Configurator(
    parameters = parameters,
    datasets = {
        "jsons": [
            f"{localdir}/datasets/local_fileset.json",
            # f"{localdir}/datasets/DATA_SingleMuon.json",    
            # f"{localdir}/datasets/DYJetsToLL_M-50.json",
            # f"{localdir}/datasets/DYJetsToLL_M10To50.json",
            # f"{localdir}/datasets/TTTo2L2Nu.json",
            # f"{localdir}/datasets/TGJets.json",
            # f"{localdir}/datasets/TTToHadronic.json",
            # f"{localdir}/datasets/TTToSemiLeptonic.json"
        ],
        "filter" : {
            "samples": ["DATA_EGamma", "DYJets", "TGJets", "GJets", "TTG", "WJets", "WG", "WWG", "WZG", "ZG", "ST"], 
            # "DATA_SinglePhoton", "TT", "ST"
            "samples_exclude" : [],
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
                          # "sf_ele_id", "sf_ele_reco",
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
                                # "sf_ele_id", "sf_ele_reco",
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
       "photon_pt" : HistConf( [Axis(coll="PhotonGood", field="pt", bins=20, start=0, stop=1000, label="$p_T^\gamma$")] ),
       "photon_eta" : HistConf( [Axis(coll="PhotonGood", field="eta", bins=10, start=-2.5, stop=2.5, label="$\eta_\gamma$")] ),
       "BJet_pt": HistConf( [Axis(coll="BJetGood", field="pt", bins=20, start=0, stop=400, label="$p_T^b$")] ),
       "BJet_eta": HistConf( [Axis(coll="BJetGood", field="eta", bins=16, start=-4, stop=4, label="$\eta_b$")] ),
       "top_pt": HistConf( [Axis(coll="top", field="pt", bins=20, start=0, stop=1000, label="$p_T^{top}$")] ),
       "top_mass": HistConf( [Axis(coll="top", field="mass", bins=20, start=0, stop=1000, label="$top_M$")] ),
       "VLT_pt": HistConf( [Axis(coll="VLT", field="pt", bins=20, start=0, stop=500, label="$p_T^{VLT}$")] ),
       "VLT_mass": HistConf( [Axis(coll="VLT", field="mass", bins=20, start=0, stop=2000, label="$VLT_M$")] ),
       "Lepton_pt" : HistConf( [Axis(coll="LeptonGood", field="pt", bins=20, start=0, stop=1000, label="$e_{p_T}$")] ),
       "WTransverse" : HistConf( [Axis(coll="events", field="W_transMass", bins=20, start=0, stop=200, label="$W_T$")] )
   }
)

