import awkward as ak
import copy
# import os
# import logging
import time
from collections import defaultdict


from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.lib.categorization import StandardSelection
from pocket_coffea.parameters.cuts import passthrough
from custom_hist_manager import CustomHistManager
from object_selector import *
from custom_cut_functions import *
from custom_configurator import CustomConfigurator
from pocket_coffea.lib.objects import (
    jet_correction,
    jet_selection,
    get_dilepton
)


class TopPartnerBaseProcessor(BaseProcessorABC):
    def __init__(self, cfg: CustomHistManager):
        super().__init__(cfg)
        self.lepton = self.cfg.lepton
        self.cfg.subsamples_reversed_map["JetFakePhoton"] = "PLJ"
        self.cfg.samples_metadata["PLJ"] = {'isMC': False}

    def load_metadata_extra(self):
        self._nEvents_primary = self.events.metadata["nevents"]

    def apply_object_preselection(self, variation):
        self.events = ak.with_field(self.events, self.lepton, 'flavor')
        self.events["LeptonGood"] = lepton_selection(
            self.events, self.lepton, self.params, id="tight"
        )
        self.events["MuonLoose"] = lepton_selection(
            self.events, "Muon", self.params, id="loose"
        )
        self.events["ElectronVeto"] = lepton_selection(
            self.events, "Electron", self.params, id="veto"
        )
        self.events["PhotonGood"] = photon_selection(
            self.events, "Photon", self.params, "SR", "LeptonGood"
        )
        self.events["PhotonCRB"] = photon_selection(
            self.events, "Photon", self.params, "CRB", "LeptonGood"
        )
        self.events["PhotonCRC"] = photon_selection(
            self.events, "Photon", self.params, "CRC", "LeptonGood"
        )
        self.events["PhotonCRD"] = photon_selection(
            self.events, "Photon", self.params, "CRD", "LeptonGood"
        )
        self.events["PhotonPLJ"] = photon_selection(
            self.events, "Photon", self.params, "PLJ", "LeptonGood"
        )
        
        LepPho = ak.with_name(
            ak.concatenate((self.events.LeptonGood, self.events.PhotonGood), axis=1),
            name='PtEtaPhiMCandidate',
        )
        self.events["LepPho"] = LepPho[ak.argsort(LepPho.pt, ascending=False)]
        self.events["JetGood"], self.jetGoodMask = jet_selection(
            self.events, "Jet", self.params, self._year, "LepPho"
        )
        self.events["BJetGood"] = btagging(
            self.events["JetGood"], self.params.btagging.working_point[self._year], self.params.object_preselection.Jet.btag
        )

        LepPhoCRB = ak.with_name(
            ak.concatenate((self.events.LeptonGood, self.events.PhotonCRB), axis=1),
            name='PtEtaPhiMCandidate',
        )
        self.events["LepPhoCRB"] = LepPhoCRB[ak.argsort(LepPhoCRB.pt, ascending=False)]
        self.events["JetGoodCRB"], self.jetGoodMask = jet_selection(
            self.events, "Jet", self.params, self._year, "LepPhoCRB"
        )
        self.events["BJetGoodCRB"] = btagging(
            self.events["JetGoodCRB"], self.params.btagging.working_point[self._year], self.params.object_preselection.Jet.btag
        )

        LepPhoCRC = ak.with_name(
            ak.concatenate((self.events.LeptonGood, self.events.PhotonCRC), axis=1),
            name='PtEtaPhiMCandidate',
        )
        self.events["LepPhoCRC"] = LepPhoCRC[ak.argsort(LepPhoCRC.pt, ascending=False)]
        self.events["JetGoodCRC"], self.jetGoodMask = jet_selection(
            self.events, "Jet", self.params, self._year, "LepPhoCRC"
        )
        self.events["BJetGoodCRC"] = btagging(
            self.events["JetGoodCRC"], self.params.btagging.working_point[self._year], self.params.object_preselection.Jet.btag
        )

        LepPhoCRD = ak.with_name(
            ak.concatenate((self.events.LeptonGood, self.events.PhotonCRD), axis=1),
            name='PtEtaPhiMCandidate',
        )
        self.events["LepPhoCRD"] = LepPhoCRD[ak.argsort(LepPhoCRD.pt, ascending=False)]
        self.events["JetGoodCRD"], self.jetGoodMask = jet_selection(
            self.events, "Jet", self.params, self._year, "LepPhoCRD"
        )
        self.events["BJetGoodCRD"] = btagging(
            self.events["JetGoodCRD"], self.params.btagging.working_point[self._year], self.params.object_preselection.Jet.btag
        )

        LepPhoPLJ = ak.with_name(
            ak.concatenate((self.events.LeptonGood, self.events.PhotonPLJ), axis=1),
            name='PtEtaPhiMCandidate',
        )
        self.events["LepPhoPLJ"] = LepPhoPLJ[ak.argsort(LepPhoPLJ.pt, ascending=False)]
        self.events["JetGoodPLJ"], self.jetGoodMask = jet_selection(
            self.events, "Jet", self.params, self._year, "LepPhoPLJ"
        )
        self.events["BJetGoodPLJ"] = btagging(
            self.events["JetGoodPLJ"], self.params.btagging.working_point[self._year], self.params.object_preselection.Jet.btag
        )

    def count_objects(self, variation):
        self.events["nPhotonGood"] = ak.num(self.events.PhotonGood)
        self.events["nPhotonCRB"] = ak.num(self.events.PhotonCRB)
        self.events["nPhotonCRC"] = ak.num(self.events.PhotonCRC)
        self.events["nPhotonCRD"] = ak.num(self.events.PhotonCRD)
        self.events["nPhotonPLJ"] = ak.num(self.events.PhotonPLJ)
        self.events["nLeptonGood"] = ak.num(self.events.LeptonGood)
        self.events["nMuonLoose"] = ak.num(self.events.MuonLoose)
        self.events["nElectronVeto"] = ak.num(self.events.ElectronVeto)
        self.events["nJetGood"] = ak.num(self.events.JetGood)
        self.events["nBJetGood"] = ak.num(self.events.BJetGood)
        self.events["nJetGoodCRB"] = ak.num(self.events.JetGoodCRB)
        self.events["nBJetGoodCRB"] = ak.num(self.events.BJetGoodCRB)
        self.events["nJetGoodCRC"] = ak.num(self.events.JetGoodCRC)
        self.events["nBJetGoodCRC"] = ak.num(self.events.BJetGoodCRC)
        self.events["nJetGoodCRD"] = ak.num(self.events.JetGoodCRD)
        self.events["nBJetGoodCRD"] = ak.num(self.events.BJetGoodCRD)
        self.events["nJetGoodPLJ"] = ak.num(self.events.JetGoodPLJ)
        self.events["nBJetGoodPLJ"] = ak.num(self.events.BJetGoodPLJ)

    # Function that defines common variables employed in analyses and save them as attributes of `events`
    def define_common_variables_before_presel(self, variation):
        pass

    def define_common_variables_after_presel(self, variation):        
        # self.events["MET_"] = ak.zip({
        #                         "pt": self.events.MET.pt,
        #                         "eta": ak.zeros_like(self.events.MET.pt),
        #                         "phi": self.events.MET.phi,
        #                         "mass": ak.zeros_like(self.events.MET.pt),
        #                         "charge": ak.zeros_like(self.events.MET.pt),
        #                         },with_name="PtEtaPhiMCandidate")
        self.events["neutrino"] = calculateNu4vec(self.events.LeptonGood, self.events.MET)
        self.events["W_transMass"] = np.sqrt(2*self.events.LeptonGood.pt*self.events.MET.pt*(1-np.cos(self.events.LeptonGood.delta_phi(self.events.MET))))

        # deltaR = self.events.BJetGood.delta_r(self.events.MuonGood)
        # min_deltaR = ak.argmin(deltaR, axis=1, keepdims=True)
        # self.events["b_jet"] = self.events.BJetGood[min_deltaR]
        
        # top reconstruction: W(mu+nu)+b
        self.events["top"] = self.events.LeptonGood
        # self.events["top_m"] = self.events.top.mass
        # self.events["top_pt"] = self.events.top.pt
        # self.events["top_eta"] = self.events.top.eta
        # self.events["top_phi"] = self.events.top.phi

        self.events["VLT"] = self.events.LeptonGood

    def process_extra_after_presel(self, variation):

        if not self._isMC:
            self.regions = StandardSelection({"PLJ": [PLJ_cut], "CRB": [CRB_cut], "CRC": [CRC_cut], "CRD": [CRD_cut]})
            self.regions.prepare(
                events=self.events,
                processor_params={},
            )

            pt_intervals = {
                "[30, 40]": [30, 40],
                "[40, 50]": [40, 50],
                "[50, 70]": [50, 70],
                "[70, 100]": [70, 100],
                "[100, 140]": [100, 140],
                "[140, 200]": [140, 200],
                "[200, np.inf]": [200, np.inf]
            }
    
            self.output["nevents"] = {}
            self.output["nevents_dataset"] = {}
            
            for region, mask in self.regions.get_masks():
                masked_events = self.events[mask]
                photon = ak.flatten(getattr(masked_events, "Photon{}".format(region)))
                self.output["nevents"][region] = {}
                self.output["nevents_dataset"][region] = {pt: {} for pt in pt_intervals.keys()}
                for pt, pt_interval in pt_intervals.items():
                    selected_events = masked_events[(photon.pt >= pt_interval[0]) & (photon.pt < pt_interval[1])]
                    self.output["nevents"][region][pt] = len(selected_events)
                    self.output["nevents_dataset"][region][pt][self._dataset] = len(selected_events)
            

    def define_categories(self, variation):
        '''
        The function saves all the cut masks internally, in order to use them later
        to define categories (groups of cuts.).

        The categorization objects takes care of the details of the caching of the mask
        and expose a common interface.

        Moreover it computes the cut masks defining the subsamples for the current
        chunks and store them in the `self.subsamples` attribute for later use.
        '''

        # We make sure that for each category the list of cuts is unique in the Configurator validation
        self._categories.prepare(
            events=self.events,
            processor_params=self.params,
            year=self._year,
            sample=self._sample,
            isMC=self._isMC,
        )

        self._subsamples[self._sample].prepare(
            events=self.events,
            processor_params=self.params,
            year=self._year,
            sample=self._sample,
            isMC=self._isMC,
        )
        
        self.PLJ_category.prepare(
            events=self.events,
            processor_params=self.params,
            year=self._year,
            sample=self._sample,
            isMC=self._isMC,
        )

        self.PLJ_subsample["PLJ"].prepare(
            events=self.events,
            processor_params=self.params,
            year=self._year,
            sample=self._sample,
            isMC=self._isMC,
        )

    def define_histograms(self):
        '''
        Initialize the HistManager.
        Redefine to customize completely the creation of the histManager.
        Only one HistManager is created for all the subsamples.
        The subsamples masks are passed to `fill_histogram` and used internally.
        '''
        self.hists_manager = CustomHistManager(
            self.cfg.variables,
            self._year,
            self._sample,
            self._subsamples[self._sample].keys(),
            self._categories,
            variations_config=self.cfg.variations_config[self._sample] if self._isMC else None,
            processor_params=self.params,
            weights_manager=self.weights_manager if self._isMC else None,
            custom_axes=self.custom_axes,
            isMC=self._isMC,
        )            
            
    def define_histograms_extra(self):
        '''
        Function that get called after the creation of the HistManager.
        The user can redefine this function to manipulate the HistManager
        histogram configuration to add customizations directly to the histogram
        objects before the filling.

        This function should also be redefined to fill the `self.custom_histogram_fields`
        that are passed to the histogram filling.
        '''
        # Here we define a new HistManager for PLJ sample
        # We use PocketCoffea HistManager with new subsample
        self.PLJ_category = StandardSelection({"SR": [PLJ_cut]})
        self.PLJ_subsample = {"PLJ": StandardSelection({"JetFakePhoton": [passthrough]})}
        
        if not self._isMC:
            
            self.PLJsample_hists_manager = CustomHistManager(
                self.cfg.variables,
                self._year,
                "PLJ",
                self.PLJ_subsample["PLJ"].keys(),
                self.PLJ_category,
                variations_config=self.cfg.variations_config[self._sample] if self._isMC else None,
                processor_params=self.params,
                weights_manager=self.weights_manager,
                custom_axes=self.custom_axes,
                isMC=self._isMC,
            )
        

    # def fill_histograms(self, variation):
    #     '''Function which fill the histograms for each category and variation,
    #     throught the HistManager.
    #     '''
    #     if not len(self.SR_events):
    #         # Filling the autofill=True histogram automatically
    #         # Calling hist manager with the subsample masks
    #         self.hists_manager.fill_histograms(
    #             self.SR_events,
    #             self.SR_category,
    #             subsamples=self._subsamples[self._sample],
    #             shape_variation=variation,
    #             custom_fields=self.custom_histogram_fields,
    #         )
    #         # Saving the output for each sample/subsample
    #         for subs in self._subsamples[self._sample].keys():
    #             # When we loop on all the subsample we need to format correctly the output if
    #             # there are no subsamples
    #             if self._hasSubsamples:
    #                 name = f"{self._sample}__{subs}"
    #             else:
    #                 name = self._sample
    #             for var, H in self.hists_manager.get_histograms(subs).items():
    #                 self.output["variables"][var][name][self._dataset] = H


    def fill_histograms_extra(self, variation):
        '''
        The function get called after the filling of the default histograms.
        Redefine it to fill custom histograms
        '''
        if not self._isMC:
            self.PLJsample_hists_manager.fill_histograms(
                self.events,
                self.PLJ_category,
                subsamples=self.PLJ_subsample["PLJ"],
                shape_variation=variation,
                custom_fields=self.custom_histogram_fields,
                smpl="PLJ"
            )
            # Saving the output for each sample/subsample
            subsample_name = "JetFakePhoton"
            dataset_name = self._dataset + "_PLJ"
            for var, H in self.PLJsample_hists_manager.get_histograms("JetFakePhoton").items():
                self.output["variables"][var][subsample_name][dataset_name] = H
        
    def process(self, events: ak.Array):
        '''
        This function get called by Coffea on each chunk of NanoAOD file.
        The processing steps of PocketCoffea are defined in this function.

        Customization points for user-defined processor are provided as
        `_extra` functions. By redefining those functions the user can change the behaviour
        of the processor in some predefined points.

        The processing goes through the following steps:

          - load metadata
          - Skim events (first masking of events):
              HLT triggers should be applied here, but their use is left to the configuration,
              and not hardcoded in the processor.
          - apply object preselections
          - count objects
          - apply event preselection (secondo masking of events)
          - define categories
          - define weights
          - define histograms
          - count events in each category
        '''
        self.start_time = time.time()
        self.events = events
        # Define the accumulator instance for this chunk
        self.output = copy.deepcopy(self.output_format)

        ###################
        # At the beginning of the processing the initial number of events
        # and the sum of the genweights is stored for later use
        #################
        self.load_metadata()
        self.load_metadata_extra()

        self.nEvents_initial = self.nevents
        self.output['cutflow']['initial'][self._dataset] = self._nEvents_primary
        if self._isMC:
            # This is computed before any preselection
            if not self._isSkim:
                self.output['sum_genweights'][self._dataset] = self._nEvents_primary
            else:
                # If the dataset is a skim, the sumgenweights are rescaled
                self.output['sum_genweights'][self._dataset] = ak.sum(self.events.skimRescaleGenWeight * self.events.genWeight)
            #FIXME: handle correctly the skim for the sum_signOf_genweights
            self.output['sum_signOf_genweights'][self._dataset] = ak.sum(np.sign(self.events.genWeight))
                
        ########################
        # Then the first skimming happens.
        # Events that are for sure useless are removed.
        # The user can specify in the configuration a function to apply for the skiming.
        # BE CAREFUL: objects are not corrected and cleaned at this stage, the skimming
        # selections MUST be loose and inclusive w.r.t the final selections.
        #########################
        # Customization point for derived workflows before skimming
        self.process_extra_before_skim()
        # MET filter, lumimask, + custom skimming function
        self.skim_events()
        if not self.has_events:
            return self.output

        if self.cfg.save_skimmed_files:
            self.export_skimmed_chunk()
            return self.output

        #########################
        # After the skimming we apply the object corrections and preselection
        # Doing so we avoid to compute them on the full NanoAOD dataset
        #########################

        self.process_extra_after_skim()
        # Define and load the weights manager
        self.define_weights()
        # Create the HistManager and ColumnManager before systematic variations
        self.define_custom_axes_extra()
        self.define_histograms()
        self.define_histograms_extra()
        self.define_column_accumulators()
        self.define_column_accumulators_extra()

        for variation in self.get_shape_variations():
            # Apply preselections
            self.apply_object_preselection(variation)
            self.count_objects(variation)
            # Compute variables after object preselection
            self.define_common_variables_before_presel(variation)
            # Customization point for derived workflows after preselection cuts
            self.process_extra_before_presel(variation)

            # This will remove all the events not passing preselection
            # from further processing
            self.apply_preselections(variation)

            # If not events remains after the preselection we skip the chunk
            if not self.has_events:
                continue

            ##########################
            # After the preselection cuts has been applied more processing is performend
            ##########################
            # Customization point for derived workflows after preselection cuts
            self.define_common_variables_after_presel(variation)
            self.process_extra_after_presel(variation)

            # This function applies all the cut functions in the cfg file
            # Each category is an AND of some cuts.
            self.define_categories(variation)

            # Weights
            self.compute_weights(variation)
            self.compute_weights_extra(variation)

            # Fill histograms
            self.fill_histograms(variation)
            self.fill_histograms_extra(variation)
            self.fill_column_accumulators(variation)
            self.fill_column_accumulators_extra(variation)

            # Count events
            if variation == "nominal":
                self.count_events(variation)

        self.stop_time = time.time()
        self.save_processing_metadata()
        return self.output

    def postprocess(self, accumulator):
        '''
        The function is called by coffea at the end of the processing.
        The default function calls the `rescale_sumgenweights` function to rescale the histograms
        and `sumw` metadata using the sum of the genweights computed without preselections
        for each dataset.

        Moreover the function saves in the output a dictionary of metadata
        with the full description of the datasets taken from the configuration.

        To add additional customatizaion redefine the `postprocessing` function,
        but remember to include a super().postprocess() call.
        '''
        
        if not self.cfg.do_postprocessing:
            return accumulator

                
        # Saving dataset metadata directly in the output file reading from the config
        dmeta = accumulator["datasets_metadata"] = {
            "by_datataking_period": {},
            "by_dataset": defaultdict(dict)
        }

        for dataset in accumulator["cutflow"]["initial"].keys():
            df = self.cfg.filesets[dataset]
            #copying the full metadata of the used samples in the output per direct reference
            dmeta["by_dataset"][dataset] = df["metadata"]
            # now adding by datataking period
            sample = df["metadata"]["sample"]
            year = df["metadata"]["year"]
            if year not in dmeta["by_datataking_period"]:
                dmeta["by_datataking_period"][year] = defaultdict(set)

            if self.cfg.has_subsamples[sample]:
                for subsam in self._subsamples[sample].keys():
                    dmeta["by_datataking_period"][year][f"{sample}__{subsam}"].add(dataset)
            else:
                dmeta["by_datataking_period"][year][sample].add(dataset)

        # Rescale the histograms and sumw using the sum of the genweights
        if not self.workflow_options.get("donotscale_sumgenweights", False):
            self.rescale_sumgenweights(accumulator)

        for dataset in accumulator["cutflow"]["initial"].keys():
            df = self.cfg.filesets[dataset]
            #copying the full metadata of the used samples in the output per direct reference
            dmeta["by_dataset"][dataset] = df["metadata"]
            year = df["metadata"]["year"]
            
            # These lines add information and description for JetFaksePhoton Sample
            is_mc = df["metadata"]["isMC"]
            if is_mc == "False":
                subsample_name = "JetFakePhoton"
                dataset_name = dataset + "_PLJ"
                PLJ_metadata = {
                    "sample": "JetFaksePhoton",
                    "year": year,
                    "isMC": "True"
                }
                dmeta["by_dataset"][dataset_name] = PLJ_metadata
                if year not in dmeta["by_datataking_period"]:
                    dmeta["by_datataking_period"][year] = defaultdict(set)
                dmeta["by_datataking_period"][year]["JetFakePhoton"].add(dataset_name)

        diction = {pt:{} for _, pt_int in accumulator["nevents"].items() for pt in pt_int}
        for region, dic in accumulator["nevents"].items():
            for pt_interval, _ in dic.items():
                diction[pt_interval][region] = dic[pt_interval]
        
        EF = accumulator["EF"] = {}
        for pt_interval, dicti in diction.items():
            EF[pt_interval] = ((dicti["CRB"]*dicti["CRC"])/dicti["CRD"])/dicti["PLJ"]

        return accumulator