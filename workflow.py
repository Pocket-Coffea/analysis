import awkward as ak
import copy
# import os
# import logging
import time
from collections import defaultdict
import numpy as np


from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.lib.categorization import StandardSelection
from pocket_coffea.parameters.cuts import passthrough
from custom_hist_manager import CustomHistManager
from object_selector import *
from custom_cut_functions import *
from helper_functions import *
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.objects import (
    jet_correction,
    jet_selection,
    get_dilepton
)


class TopPartnerBaseProcessor(BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.lepton = "Electron"

    def apply_object_preselection(self, variation):
        self.events = ak.with_field(self.events, self.lepton, 'flavor')
        self.events["Photon"] = ak.with_field(self.events.Photon, self.events.Photon.pt * self.events.Photon.pfRelIso03_chg, "chIso")
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

    def define_common_variables_after_presel(self, variation):        

        self.events["neutrino"] = calculateNu4vec(self.events.LeptonGood, self.events.MET)
        self.events["W_transMass"] = np.sqrt(2*self.events.LeptonGood.pt*self.events.MET.pt*(1-np.cos(self.events.LeptonGood.delta_phi(self.events.MET))))
        
        # top reconstruction: W(mu+nu)+b
        # self.events["top"] = self.events.LeptonGood
        # self.events["VLT"] = self.events.LeptonGood

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
            self._hasSubsamples,
            self._subsamples[self._sample].keys(),
            self._categories,
            variations_config=self.cfg.variations_config[self._sample] if self._isMC else None,
            processor_params=self.params,
            weights_manager=self.weights_manager,
            calibrators_manager=self.calibrators_manager,
            custom_axes=self.custom_axes,
            isMC=self._isMC,
        )
