import awkward as ak
import copy
# import os
# import logging
import time


from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.hist_manager import Axis
from object_selector import *
from pocket_coffea.lib.objects import (
    jet_correction,
    jet_selection,
    btagging,
    get_dilepton
)


class TopPartnerBaseProcessor(BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.final_state = "Had"

    def load_metadata_extra(self):
        self._nEvents_primary = self.events.metadata["nevents"]

    def apply_object_preselection(self, variation):

        #####################################################                      #############################################
                                                               #Hadronic Channel#
        #####################################################                      #############################################
        if self.final_state == "Had":
            self.events["MuonLoose"] = lepton_selection(
                self.events, "Muon", self.params, id="loose"
            )
            self.events["ElectronVeto"] = lepton_selection(
                self.events, "Electron", self.params, id="veto"
            )
            self.events["PhotonSR"] = photon_selection(
                self.events, "Photon", self.params, region = "SR"
            )
            self.events["PhotonCRA"] = photon_selection(
                self.events, "Photon", self.params, region = "CRB"
            )
            self.events["PhotonCRB"] = photon_selection(
                self.events, "Photon", self.params, region = "CRB"
            )
            self.events["PhotonCRC"] = photon_selection(
                self.events, "Photon", self.params, region = "CRC"
            )
            self.events["PhotonCRD"] = photon_selection(
                self.events, "Photon", self.params, region = "CRD"
            )
            self.events["PhotonPLJ"] = photon_selection(
                self.events, "Photon", self.params, "PLJ"
            )

            
            self.events["JetGood"], self.jetGoodMask = jet_selection(
                self.events, "Jet", self.params, self._year, "PhotonSR"
            )
            # index = ak.Array([[j for j in range(num)] for num in ak.num(self.events.JetGood)]) # ak.local_index does the same
            # self.events["JetGood"] = ak.with_field(self.events["JetGood"], index, "index")
            self.events["BJetGood"] = btagging(
                self.events["JetGood"], self.params.btagging.working_point[self._year], wp=self.params.object_preselection.Jet.btag.wp
            )
            self.events["JetGood_NotB"] =btagging(
                self.events["JetGood"], self.params.btagging.working_point[self._year], wp=self.params.object_preselection.Jet.btag.wp, veto=True
            )
                
        #####################################################                      #############################################
                                                                 #Muon Channel#
        #####################################################                      #############################################
        if self.final_state == "Mu":
            self.events["MuonGood"] = lepton_selection(
                self.events, "Muon", self.params, id="tight"
            )
            self.events["MuonLoose"] = lepton_selection(
                self.events, "Muon", self.params, id="loose"
            )
            self.events["ElectronVeto"] = lepton_selection(
                self.events, "Electron", self.params, id="veto"
            )
            self.events["PhotonSR"] = photon_selection(
                self.events, "Photon", self.params, "SR", "MuonGood"
            )
            self.events["PhotonCRB"] = photon_selection(
                self.events, "Photon", self.params, "CRB", "MuonGood"
            )
            self.events["PhotonCRC"] = photon_selection(
                self.events, "Photon", self.params, "CRC", "MuonGood"
            )
            self.events["PhotonCRD"] = photon_selection(
                self.events, "Photon", self.params, "CRD", "MuonGood"
            )
            self.events["PhotonPLJ"] = photon_selection(
                self.events, "Photon", self.params, "PLJ", "MuonGood"
            )
            MuPho = ak.with_name(
                ak.concatenate((self.events.MuonGood, self.events.PhotonSR), axis=1),
                name='PtEtaPhiMCandidate',
            )
            self.events["MuPho"] = MuPho[ak.argsort(MuPho.pt, ascending=False)]
            self.events["JetGood"], self.jetGoodMask = jet_selection(
                self.events, "Jet", self.params, self._year, "MuPho"
            )
            self.events["BJetGood"] = btagging(
                self.events["JetGood"], self.params.btagging.working_point[self._year], wp=self.params.object_preselection.Jet.btag.wp
            )
        #####################################################                      #############################################
                                                               #Electron Channel#
        #####################################################                      #############################################
        if self.final_state == "Ele":
            self.events["MuonLoose"] = lepton_selection(
                self.events, "Muon", self.params, id="loose"
            )
            self.events["ElectronGood"] = lepton_selection(
                self.events, "Electron", self.params, id="tight"
            )
            self.events["ElectronVeto"] = lepton_selection(
                self.events, "Electron", self.params, id="veto"
            )
            self.events["PhotonSR"] = photon_selection(
                self.events, "Photon", self.params, "SR", "ElectronGood"
            )
            self.events["PhotonCRB"] = photon_selection(
                self.events, "Photon", self.params, "CRB", "ElectronGood"
            )
            self.events["PhotonCRC"] = photon_selection(
                self.events, "Photon", self.params, "CRC", "ElectronGood"
            )
            self.events["PhotonCRD"] = photon_selection(
                self.events, "Photon", self.params, "CRD", "ElectronGood"
            )
            self.events["PhotonPLJ"] = photon_selection(
                self.events, "Photon", self.params, "PLJ", "ElectronGood"
            )
            ElePho = ak.with_name(
                ak.concatenate((self.events.ElectronGood, self.events.PhotonSR), axis=1),
                name='PtEtaPhiMCandidate',
            )
            self.events["ElePho"] = ElePho[ak.argsort(ElePho.pt, ascending=False)]
            self.events["JetGood"], self.jetGoodMask = jet_selection(
                self.events, "Jet", self.params, self._year, "ElePho"
            )
            self.events["BJetGood"] = btagging(
                self.events["JetGood"], self.params.btagging.working_point[self._year], self.params
            )


    def count_objects(self, variation):
        self.events["nPhotonSR"] = ak.num(self.events.PhotonSR)
        self.events["nPhotonCRB"] = ak.num(self.events.PhotonCRB)
        self.events["nPhotonCRC"] = ak.num(self.events.PhotonCRC)
        self.events["nPhotonCRD"] = ak.num(self.events.PhotonCRD)
        self.events["nPhotonPLJ"] = ak.num(self.events.PhotonPLJ)
        # self.events["nMuonGood"] = ak.num(self.events.MuonGood)
        self.events["nMuonLoose"] = ak.num(self.events.MuonLoose)
        self.events["nElectronVeto"] = ak.num(self.events.ElectronVeto)
        self.events["nJetGood"] = ak.num(self.events.JetGood)
        self.events["nBJetGood"] = ak.num(self.events.BJetGood)
        self.events["nJetGood_NotB"] = ak.num(self.events.JetGood_NotB)

    # Function that defines common variables employed in analyses and save them as attributes of `events`
    def define_common_variables_before_presel(self, variation):
        # non_btag_mask = []
        # btagged_indices = self.events.BJetGood.index
        # for outer_idx, inner_array in enumerate(self.events.JetGood.index):
        #     inner_btagged_indices = btagged_indices[outer_idx]
        #     mask = ak.Array([i not in inner_btagged_indices for i in range(len(inner_array))])
        #     non_btag_mask.append(mask)
        # non_btag_mask = ak.Array(non_btag_mask)
        
        # # events["JetGood"] = ak.with_field(events["JetGood"], btagg_mask, "b_tagged")
        # # combined_jets = ak.combinations(self.events.JetGood, 3, fields=['j0', 'j1', 'j2'])
        
        # self.events["JetGood_NotB"] = self.events.JetGood[non_btag_mask]
        # self.events["nJetGood_NotB"] = ak.num(self.events.JetGood_NotB)
        pass

    def define_common_variables_after_presel(self, variation):        
        self.events["W"] = self.events.JetGood_NotB[:, 0] + self.events.JetGood_NotB[:, 1]
        self.events["top"] = self.events.JetGood_NotB[:, 0] + self.events.JetGood_NotB[:, 1] + self.events.BJetGood[:, 0]
        
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