import awkward as ak

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.hist_manager import Axis
from pocket_coffea.lib.objects import (
    jet_correction,
    jet_selection,
    btagging,
    get_dilepton,
)


class TopPartnerBaseProcessor(BaseProcessorABC):
    def __init__(self, cfg: Configurator, final_state="Mu"):
        super().__init__(cfg)
        self.final_state = final_state

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
                self.events, "Jet", self.params, year=self._year, "PhotonSR"
            )
            self.events["BJetGood"] = btagging(
                self.events["JetGood"], self.params.btagging.working_point[self._year], wp=self.params.object_preselection.Jet.btag.wp
            )
    
        #####################################################                      #############################################
                                                                 #Muon Channel#
        #####################################################                      #############################################
        if self.final_state == "Mu":
            self.events["MuonTight"] = lepton_selection(
                self.events, "Muon", self.params, id="tight"
            )
            self.events["MuonLoose"] = lepton_selection(
                self.events, "Muon", self.params, id="loose"
            )
            self.events["ElectronVeto"] = lepton_selection(
                self.events, "Electron", self.params, id="veto"
            )
            self.events["PhotonSR"] = photon_selection(
                self.events, "Photon", self.params, region = "SR", "MuonTight"
            )
            self.events["PhotonCRB"] = photon_selection(
                self.events, "Photon", self.params, region = "CRB", "MuonTight"
            )
            self.events["PhotonCRC"] = photon_selection(
                self.events, "Photon", self.params, region = "CRC", "MuonTight"
            )
            self.events["PhotonCRD"] = photon_selection(
                self.events, "Photon", self.params, region = "CRD", "MuonTight"
            )
            MuPho = ak.with_name(
                ak.concatenate((self.events.MuonTight, self.events.PhotonSR), axis=1),
                name='PtEtaPhiMCandidate',
            )
            self.events["MuPho"] = MuPho[ak.argsort(MuPho.pt, ascending=False)]
            self.events["JetGood"], self.jetGoodMask = jet_selection(
                self.events, "Jet", self.params, year=self._year, "MuPho"
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
            self.events["ElectronTight"] = lepton_selection(
                self.events, "Electron", self.params, id="tight"
            )
            self.events["ElectronVeto"] = lepton_selection(
                self.events, "Electron", self.params, id="veto"
            )
            self.events["PhotonSR"] = photon_selection(
                self.events, "Photon", self.params, region = "SR", "ElectronTight"
            )
            self.events["PhotonCRB"] = photon_selection(
                self.events, "Photon", self.params, region = "CRB", "ElectronTight"
            )
            self.events["PhotonCRC"] = photon_selection(
                self.events, "Photon", self.params, region = "CRC", "ElectronTight"
            )
            self.events["PhotonCRD"] = photon_selection(
                self.events, "Photon", self.params, region = "CRD", "ElectronTight"
            )
            ElePho = ak.with_name(
                ak.concatenate((self.events.ElectronTight, self.events.PhotonSR), axis=1),
                name='PtEtaPhiMCandidate',
            )
            self.events["ElePho"] = ElePho[ak.argsort(ElePho.pt, ascending=False)]
            self.events["JetGood"], self.jetGoodMask = jet_selection(
                self.events, "Jet", self.params, year=self._year, "ElePho"
            )
            self.events["BJetGood"] = btagging(
                self.events["JetGood"], self.params.btagging.working_point[self._year], self.params
            )


    def count_objects(self, variation):
        self.events["nPhotonGood"] = ak.num(self.events.PhotonGood)
        self.events["nMuonGood"] = ak.num(self.events.MuonGood)
        self.events["nLooseMuon"] = ak.num(self.events.LooseMuon)
        self.events["nVetoElectrons"] = ak.num(self.events.VetoElectron)
        self.events["nJetGood"] = ak.num(self.events.JetGood)
        self.events["nBJetGood"] = ak.num(self.events.BJetGood)


    # Function that defines common variables employed in analyses and save them as attributes of `events`
    def define_common_variables_before_presel(self, variation):
        self.events["JetGood_Ht"] = ak.sum(abs(self.events.JetGood.pt), axis=1)
