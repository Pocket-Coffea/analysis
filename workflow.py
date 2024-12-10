import awkward as ak

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
            self.events["PhotonCRB"] = photon_selection(
                self.events, "Photon", self.params, region = "CRB"
            )
            self.events["PhotonCRC"] = photon_selection(
                self.events, "Photon", self.params, region = "CRC"
            )
            self.events["PhotonCRD"] = photon_selection(
                self.events, "Photon", self.params, region = "CRD"
            )            
            self.events["JetGood"], self.jetGoodMask = jet_selection(
                self.events, "Jet", self.params, self._year, "PhotonSR"
            )
            self.events["BJetGood"] = btagging(
                self.events["JetGood"], self.params.btagging.working_point[self._year], wp=self.params.object_preselection.Jet.btag.wp
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


    # Function that defines common variables employed in analyses and save them as attributes of `events`
    def define_common_variables_before_presel(self, variation):
        #self.events["JetGood_Ht"] = ak.sum(abs(self.events.JetGood.pt), axis=1)
        pass
