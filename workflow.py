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
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)


    def photon_selection(events, photon, params):

        photons = events["Photons"]
        cuts = params.object_preselection[photon]
        # Requirements on pT and eta
        passes_eta = abs(photons.eta) < cuts["eta"]
        passes_transition = np.invert(( abs(photons.eta) >= 1.4442) & (abs(photons.eta) <= 1.5660))
        passes_pt = photons.pt > cuts["pt"]
        passes_id = photons.cutBased >= 2
        passes_pixelseed = ~photons.pixelSeed

        good_photons = passes_eta & passes_pt & passes_pixelseed & passes_id

        return photons[good_photons]


    def lepton_selection(events, lepton_flavour, params, Tight = False, Loose = False, Veto = False):

        leptons = events[lepton_flavour]
        cuts = params.object_preselection[lepton_flavour]
        # Requirements on pT and eta
        passes_eta = abs(leptons.eta) < cuts["eta"]
        
        if lepton_flavour == "Electron":
            # Requirements on SuperCluster eta, isolation and id
            etaSC = abs(leptons.deltaEtaSC + leptons.eta)
            passes_SC = np.invert((etaSC >= 1.4442) & (etaSC <= 1.5660))

            passes_iso = True
            
            if Tight:
                passes_pt = leptons.pt > 35
                passes_iso = leptons.pfRelIso03_all < 0.15
                passes_id = leptons.cutBased == 4
            elif Veto:
                passes_pt = leptons.pt > 20
                passes_iso = leptons.pfRelIso03_all < 0.25
                passes_id = leptons.cutBased >= 1

            good_leptons = passes_eta & passes_pt & passes_SC & passes_iso & passes_id

        elif lepton_flavour == "Muon":
            # Requirements on isolation and id
            passes_iso = leptons.pfRelIso04_all < 0.1
            if Tight:
                passes_id = leptons[cuts['tightId']] == True
                passes_pt = leptons.pt > 30
            elif Loose:
                passes_id = leptons[cuts['looseId']] == True
                passes_pt = leptons.pt > 15

            good_leptons = passes_eta & passes_pt & passes_iso & passes_id

        return leptons[good_leptons]

    def apply_object_preselection(self, variation):
        '''
        The ttHbb processor cleans
          - Electrons
          - Muons
          - Jets -> JetGood
          - BJet -> BJetGood

        '''
        # Include the supercluster pseudorapidity variable
        electron_etaSC = self.events.Electron.eta + self.events.Electron.deltaEtaSC
        self.events["Electron"] = ak.with_field(
            self.events.Electron, electron_etaSC, "etaSC"
        )
        # Build masks for selection of photons, muons, electrons, jets
        
        self.events["PhotonGood"] = photon_selection(
            self.events, "Photon", self.params
        )
        self.events["MuonGood"] = lepton_selection(
            self.events, "Muon", self.params, Tight=True
        )
        self.events["LooseMuon"] = lepton_selection(
            self.events, "Muon", self.params, loose=True
        )
        self.events["vetoelectron"] = lepton_selection(
            self.events, "electron", self.params, Veto=True
        )
        self.events["JetGood"], self.jetGoodMask = jet_selection(
            self.events, "Jet", self.params, "MuonGood"
        )
        self.events["BJetGood"] = btagging(
            self.events["JetGood"], self.params.btagging.working_point[self._year], wp=self.params.object_preselection.Jet.btag.wp)


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
