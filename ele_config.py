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
import custom_scale_factors
import helper_functions
import object_selector
cloudpickle.register_pickle_by_value(workflow)
cloudpickle.register_pickle_by_value(custom_cut_functions)
cloudpickle.register_pickle_by_value(custom_hist_manager)
cloudpickle.register_pickle_by_value(custom_weight_manager)
cloudpickle.register_pickle_by_value(custom_scale_factors)
cloudpickle.register_pickle_by_value(object_selector)
cloudpickle.register_pickle_by_value(helper_functions)

from custom_cut_functions import *
import os
localdir = os.path.dirname(os.path.abspath(__file__))

# Loading default parameters
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")

parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection_lep.yaml",
                                                  f"{localdir}/params/triggers_ele.yaml",
                                                  f"{localdir}/params/plotting.yaml",
                                                  f"{localdir}/params/photon_scale_factors.yaml",
                                                  update=True)



cfg = Configurator(
    parameters = parameters,
    datasets = {
        "jsons": [
            # f"{localdir}/datasets/filesets_test.json",
            f"{localdir}/datasets/fileset_data.json",
            f"{localdir}/datasets/fileset_MC.json",    
            f"{localdir}/datasets/fileset_Signal.json"            
        ],
        "filter" : {
            # "samples": ["WJets", "Signal_1000", # "DATA_EGamma", 
            "samples": ["DATA_EGamma",
                        "DYJets", "TGJets", "GJets", "ST", "TTG", "TT", "WJets", "WG", "WW", "WWG", "WZG", "WZ", "ZZG", "ZZ", "ZG",
                        # "Signal_600", "Signal_1000", "Signal_1500", "Signal_2000",
                        # "Signal_700", "Signal_800", "Signal_900", "Signal_1100", "Signal_1200", "Signal_1300", "Signal_1400", "Signal_1600", "Signal_1700", "Signal_1800","Signal_1900",
                       ],
            # "DATA_SinglePhoton",
            "samples_exclude" : [],
            "year": ['2018']
        }
    },

    workflow = TopPartnerBaseProcessor,

    skim = [get_nPVgood(1), eventFlags, goldenJson, # basic skims
            #get_nObj_min(1, 18., "Muon"),
            # Asking only SingleMuon triggers since we are only using SingleMuon PD data
            get_HLTsel()], 
    
    preselections = [vlt_presel],
    categories = {
        "SR": [SR_cut],
        "PLJ": [PLJ_cut],
        "CRB": [CRB_cut],
        "CRC": [CRC_cut],
        "CRD": [CRD_cut],
        "b0_SR": [b0_SR_cut],
        "b0_PLJ": [b0_PLJ_cut],
        "b0_CRB": [b0_CRB_cut],
        "b0_CRC": [b0_CRC_cut],
        "b0_CRD": [b0_CRD_cut]
    },

    weights_classes = common_weights,
    
    weights = {
        "common": {
            "inclusive": ["genWeight","lumi","XS",
                          "pileup",
                          "custom_sf_ele_id", "custom_sf_ele_reco",
                          "custom_sf_pho_id", "custom_sf_pho_pxseed",
                          "sf_btag"
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
                "inclusive": ["pileup",
                              "custom_sf_ele_id", "custom_sf_ele_reco",
                              "custom_sf_pho_id", "custom_sf_pho_pxseed",
                              "sf_btag"
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
       "photon_pt" : HistConf([Axis(coll="PhotonGood", field="pt", bins=[(60+i*20) for i in range(12)], label="$p_{T,\gamma}(GeV)$")]),
       "photon_sieie" : HistConf([Axis(coll="PhotonGood", field="sieie", bins=20, start=0.005, stop=0.02, label="$\sigma_{i\eta i\eta}$")]),
       "photon_chiso" : HistConf([Axis(coll="PhotonGood", field="chIso", bins=15, start=0, stop=10, label="Charged Hadron Isolation")] ),
       "photon_eta" : HistConf( [Axis(coll="PhotonGood", field="eta", bins=6, start=-1.5, stop=1.5, label="$\eta_\gamma$")]),
       "Electron_pt" : HistConf([Axis(coll="LeptonGood", field="pt", bins=[35 + (i*15) for i in range(8)]+[175, 195, 215, 235, 260], label="$p_{T,e}(GeV)$")]),
       "WTransverse" : HistConf([Axis(coll="events", field="W_transMass", bins=15, start=0, stop=240, label="$mW_T(GeV)$")]),
       **count_hist(name="nJets", coll="JetGood",bins=10, start=0, stop=10),
       "BJet_pt": HistConf([Axis(coll="BJetGood", field="pt", bins=[(30+i*30) for i in range(8)], label="$p_{T,b}(GeV)$")],
                           exclude_categories=["b0_SR", "b0_PLJ", "b0_CRB", "b0_CRC", "b0_CRD"]
                          ),
       "BJet_eta": HistConf([Axis(coll="BJetGood", field="eta", bins=10, start=-2.5, stop=2.5, label="$\eta_b$")],
                           exclude_categories=["b0_SR", "b0_PLJ", "b0_CRB", "b0_CRC", "b0_CRD"]
                          ),
       "top_pt": HistConf([Axis(coll="top", field="pt", bins=[i*20 for i in range(16)], label="$p_{T,top}(GeV)$")],
                           exclude_categories=["b0_SR", "b0_PLJ", "b0_CRB", "b0_CRC", "b0_CRD"]
                          ),
       "top_mass": HistConf([Axis(coll="top", field="mass", bins=[(60+i*30) for i in range(21)], label="$top_M(GeV)$")],
                           exclude_categories=["b0_SR", "b0_PLJ", "b0_CRB", "b0_CRC", "b0_CRD"]
                          ),
       "VLT_pt": HistConf([Axis(coll="VLT", field="pt", bins=[i*15 for i in range(18)], overflow=True, label="$p_{T,VLT}(GeV)$")],
                           exclude_categories=["b0_SR", "b0_PLJ", "b0_CRB", "b0_CRC", "b0_CRD"]
                          ),
       "VLT_mass": HistConf([Axis(coll="VLT", field="mass", bins=[(i*100) for i in range(9)]+[(1000+j*200) for j in range(6)], label="$VLT_M(GeV)$")],
                           exclude_categories=["b0_SR", "b0_PLJ", "b0_CRB", "b0_CRC", "b0_CRD"]
                          )
   }
)

