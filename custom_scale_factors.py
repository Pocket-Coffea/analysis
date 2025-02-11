# from pocket_coffea.lib.scale_factors import (
#     get_ele_sf,
#     get_mu_sf
# )
import awkward as ak
import numpy as np
import correctionlib

def get_ele_sf(
        params, year, pt, eta, phi, counts=None, key='', pt_region=None, variations=["nominal"]
):
    '''
    This function computes the per-electron reco or id SF.
    If 'reco', the appropriate corrections are chosen by using the argument `pt_region`.
    '''
    electronSF = params["lepton_scale_factors"]["electron_sf"]
    # translate the `year` key into the corresponding key in the correction file provided by the EGM-POG
    year_pog = electronSF["era_mapping"][year]

    if key in ['reco', 'id']:
        electron_correctionset = correctionlib.CorrectionSet.from_file(
            electronSF.JSONfiles[year]["file"]
        )
        map_name = electronSF.JSONfiles[year]["name"]

        if key == 'reco':
            sfname = electronSF.JSONfiles[year]["reco"][pt_region]
        elif key == 'id':
            sfname = params.object_preselection["Electron"]["id"]
        
        if year in ["2023_preBPix", "2023_postBPix"]:
            # Starting from 2023 SFs require the phi:
            if key == 'reco' and pt_region == 'pt_lt_20':
                # It also appears that for RecoBelow20 SFs the eta must be positive (absolute value).
                eta_np = np.abs(eta.to_numpy())
            else:
                eta_np = eta.to_numpy()
                
            sf = electron_correctionset[map_name].evaluate(
                year_pog, "sf", sfname, eta_np, pt.to_numpy(),  phi.to_numpy()
            )
            sfup = electron_correctionset[map_name].evaluate(
                year_pog, "sfup", sfname, eta_np, pt.to_numpy(), phi.to_numpy()
            )
            sfdown = electron_correctionset[map_name].evaluate(
                year_pog, "sfdown", sfname, eta_np, pt.to_numpy(), phi.to_numpy()
            )
        else:
            # All other eras do not need phi:    
            sf = electron_correctionset[map_name].evaluate(
                year_pog, "sf", sfname, eta.to_numpy(), pt.to_numpy()
            )
            sfup = electron_correctionset[map_name].evaluate(
                year_pog, "sfup", sfname, eta.to_numpy(), pt.to_numpy()
            )
            sfdown = electron_correctionset[map_name].evaluate(
                year_pog, "sfdown", sfname, eta.to_numpy(), pt.to_numpy()
            )
        # The unflattened arrays are returned in order to have one row per event.
        return (
            ak.unflatten(sf, counts),
            ak.unflatten(sfup, counts),
            ak.unflatten(sfdown, counts),
        )
    else:
        raise Exception(f"Invalid key `{key}` for get_ele_sf. Available keys are 'reco', 'id'.")



def sf_ele_reco(params, events, year):
    '''
    This function computes the per-electron reco SF and returns the corresponding per-event SF, obtained by multiplying the per-electron SF in each event.
    Additionally, also the up and down variations of the SF are returned.
    Electrons are split into two categories based on a pt cut depending on the Run preiod, so that the proper SF is applied.
    '''
    coll = "LeptonGood"
    ele_pt = events[coll].pt
    ele_eta = events[coll].etaSC # This is added on top of NanoAOD
    ele_phi = events[coll].phi

    pt_ranges = []
    if year in ['2016_PreVFP', '2016_PostVFP','2017','2018']:
        pt_ranges += [("pt_lt_20", (ele_pt < 20)), 
                      ("pt_gt_20", (ele_pt >= 20))]
    elif year in ["2022_preEE", "2022_postEE", "2023_preBPix", "2023_postBPix"]:
        pt_ranges += [("pt_lt_20", (ele_pt < 20)), 
                      ("pt_gt_20_lt_75", (ele_pt >= 20) & (ele_pt < 75)), 
                      ("pt_gt_75", (ele_pt >= 75))]
    else:
        raise Exception("For chosen year "+year+" sf_ele_reco are not implemented yet")
    
    sf_reco, sfup_reco, sfdown_reco = [], [], []

    for pt_range_key, pt_range in pt_ranges:
        ele_pt_inPtRange = ak.flatten(ele_pt[pt_range])
        ele_eta_inPtRange = ak.flatten(ele_eta[pt_range])
        ele_phi_inPtRange = ak.flatten(ele_phi[pt_range])
        ele_counts_inPtRange = ak.num(ele_pt[pt_range])

        sf_reco_inPtRange, sfup_reco_inPtRange, sfdown_reco_inPtRange = get_ele_sf(
            params,
            year,
            ele_pt_inPtRange,
            ele_eta_inPtRange,
            ele_phi_inPtRange,
            ele_counts_inPtRange,
            'reco',
            pt_range_key,
        )
        
        sf_reco.append(sf_reco_inPtRange)
        sfup_reco.append(sfup_reco_inPtRange)
        sfdown_reco.append(sfdown_reco_inPtRange)

    sf_reco = ak.prod(
        ak.concatenate(sf_reco, axis=1), axis=1
    )
    sfup_reco = ak.prod(
        ak.concatenate(sfup_reco, axis=1), axis=1
    )
    sfdown_reco = ak.prod(
        ak.concatenate(sfdown_reco, axis=1), axis=1
    )

    return sf_reco, sfup_reco, sfdown_reco


def sf_ele_id(params, events, year):
    '''
    This function computes the per-electron id SF and returns the corresponding per-event SF, obtained by multiplying the per-electron SF in each event.
    Additionally, also the up and down variations of the SF are returned.
    '''
    coll = "LeptonGood"
    ele_pt = events[coll].pt
    ele_eta = events[coll].etaSC
    ele_phi = events[coll].phi

    ele_pt_flat, ele_eta_flat, ele_phi_flat, ele_counts = (
        ak.flatten(ele_pt),
        ak.flatten(ele_eta),
        ak.flatten(ele_phi),
        ak.num(ele_pt),
    )

    sf_id, sfup_id, sfdown_id = get_ele_sf(
        params, year, ele_pt_flat, ele_eta_flat, ele_phi_flat, ele_counts, 'id'
    )
        

    # The SF arrays corresponding to the electrons are multiplied along the electron axis in order to obtain a per-event scale factor.
    return ak.prod(sf_id, axis=1), ak.prod(sfup_id, axis=1), ak.prod(sfdown_id, axis=1)


def sf_mu(params, events, year, key=''):
    '''
    This function computes the per-muon id SF and returns the corresponding per-event SF, obtained by multiplying the per-muon SF in each event.
    Additionally, also the up and down variations of the SF are returned.
    '''
    coll = "LeptonGood"
    mu_pt = events[coll].pt
    mu_eta = events[coll].eta

    # Since `correctionlib` does not support jagged arrays as an input, the pt and eta arrays are flattened.
    mu_pt_flat, mu_eta_flat, mu_counts = (
        ak.flatten(mu_pt),
        ak.flatten(mu_eta),
        ak.num(mu_pt),
    )
    sf, sfup, sfdown = get_mu_sf(params, year, mu_pt_flat, mu_eta_flat, mu_counts, key)

    # The SF arrays corresponding to all the muons are multiplied along the
    # muon axis in order to obtain a per-event scale factor.
    return ak.prod(sf, axis=1), ak.prod(sfup, axis=1), ak.prod(sfdown, axis=1)

def get_pho_sf(params, year, pt= None, eta=None, r9=None, counts=None, key='', pxbin= None):

    photonSF = params["photon_scale_factors"]
    # translate the `year` key into the corresponding key in the correction file provided by the EGM-POG
    year_pog = photonSF["era_mapping"][year]

    if key in ['pxseed', 'id']:
        photon_correctionset = correctionlib.CorrectionSet.from_file(
            photonSF.JSONfiles[year]["file"]
        )
        map_name = photonSF.JSONfiles[year][key]
        sfname = params.object_preselection["Photon"]["id"]

        if key == 'pxseed':        
            sf = photon_correctionset[map_name].evaluate(
                year_pog, "sf", sfname, pxbin
            )
            sfup = photon_correctionset[map_name].evaluate(
                year_pog, "sfup", sfname, pxbin
            )
            sfdown = photon_correctionset[map_name].evaluate(
                year_pog, "sfdown", sfname, pxbin
            )
            # The unflattened arrays are returned in order to have one row per event.
            return (
                sf,
                sfup,
                sfdown
            )
        elif key == 'id':
            # All other eras do not need phi:    
            sf = photon_correctionset[map_name].evaluate(
                year_pog, "sf", sfname, eta.to_numpy(), pt.to_numpy()
            )
            sfup = photon_correctionset[map_name].evaluate(
                year_pog, "sfup", sfname, eta.to_numpy(), pt.to_numpy()
            )
            sfdown = photon_correctionset[map_name].evaluate(
                year_pog, "sfdown", sfname, eta.to_numpy(), pt.to_numpy()
            )
            # The unflattened arrays are returned in order to have one row per event.
            return (
                ak.unflatten(sf, counts),
                ak.unflatten(sfup, counts),
                ak.unflatten(sfdown, counts),
            )
    else:
        raise Exception(f"Invalid key `{key}` for get_pho_sf. Available keys are 'pxseed', 'id'.")



def sf_pho_pxseed(params, events, year):
    
    coll = "PhotonGood"
    pho_pt = events[coll].pt
    pho_eta = events[coll].eta
    pho_r9 = events[coll].r9

    pxbins = {
        "EBHighR9": [(pho_eta < 1.4442), (pho_r9 >= 0.96)],
        "EBLowR9": [(pho_eta < 1.4442), (pho_r9 < 0.96)],
        "EEHighR9": [(pho_eta > 1.566), (pho_r9 >= 0.96)],
        "EELowR9": [(pho_eta > 1.566), (pho_r9 < 0.96)],
    }

    sf_pxseed, sfup_pxseed, sfdown_pxseed = [], [], []

    for pxbin, mask in pxbins.items():
        identity = ak.ones_like(pho_pt[mask[0] & mask[1]])
    
        sf_pxseed_pxbin, sfup_pxseed_pxbin, sfdown_pxseed_pxbin = get_pho_sf(
            params, year, key='pxseed', pxbin=pxbin
        )

        sf_pxseed.append(sf_pxseed_pxbin * identity)
        sfup_pxseed.append(sfup_pxseed_pxbin * identity)
        sfdown_pxseed.append(sfdown_pxseed_pxbin * identity)

    sf_pxseed = ak.prod(
        ak.concatenate(sf_pxseed, axis=1), axis=1
    )
    sfup_pxseed = ak.prod(
        ak.concatenate(sfup_pxseed, axis=1), axis=1
    )
    sfdown_pxseed = ak.prod(
        ak.concatenate(sfdown_pxseed, axis=1), axis=1
    )

    return sf_pxseed, sfup_pxseed, sfdown_pxseed


def sf_pho_id(params, events, year):
    
    coll = "PhotonGood"
    pho_pt = events[coll].pt
    pho_eta = events[coll].eta
    pho_r9 = events[coll].r9

    pho_pt_flat, pho_eta_flat, pho_r9_flat, pho_counts = (
        ak.flatten(pho_pt),
        ak.flatten(pho_eta),
        ak.flatten(pho_r9),
        ak.num(pho_pt),
    )

    sf_id, sfup_id, sfdown_id = get_pho_sf(
        params, year, pho_pt_flat, pho_eta_flat, pho_r9_flat, pho_counts, 'id'
    )
        

    # The SF arrays corresponding to the photons are multiplied along the photon axis in order to obtain a per-event scale factor.
    return ak.prod(sf_id, axis=1), ak.prod(sfup_id, axis=1), ak.prod(sfdown_id, axis=1)
