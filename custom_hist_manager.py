from pocket_coffea.lib.hist_manager import Axis, HistManager, get_hist_axis_from_config
from custom_weight_manager import ExtrapolationFactor
import hist
import awkward as ak
from collections import defaultdict
from coffea.analysis_tools import PackedSelection
from typing import List, Tuple
from dataclasses import dataclass, field
from copy import deepcopy
import logging


class CustomHistManager(HistManager):

    def __prefetch_weights(self, category, shape_variation):
        return self._HistManager__prefetch_weights(category, shape_variation)
        
    def fill_histograms(
        self,
        events,
        categories,
        shape_variation="nominal",
        subsamples=None,  # This is a dictionary with name:ak.Array(bool)
        custom_fields=None,
        custom_weight=None,  # it should be a dictionary {variable:weight}
    ):
        '''
        We loop on the configured histograms only
        Doing so the catergory, sample, variation selections are handled correctly (by the constructor).

        Custom_fields is a dict of additional array. The expected lenght of the first dimension is the number of
        events. The categories mask will be applied.
        '''

        # Preloading weights BOTH FOR data and MC
        weights = {}
        for category in self.available_categories:
            weights[category] = self.__prefetch_weights(category, shape_variation)
            
        # Cleaning the weights cache decorator between calls.
        self._weights_cache.clear()
        # Looping on the histograms to read the values only once
        # Then categories, subsamples and weights are applied and masked correctly

        # ASSUNTION, the histograms are the same for each subsample
        # we can take the configuration of the first subsample
        for name, histo in self.histograms[self.subsamples[0]].items():
            # logging.info(f"\thisto: {name}")
            if not histo.autofill:
                continue
            if histo.metadata_hist:
                continue  # TODO dedicated function for metadata histograms

            # Check if a shape variation is under processing and
            # if the current histogram does not require that
            if (
                shape_variation != "nominal"
                and shape_variation not in histo.hist_obj.axes["variation"]
            ):
                continue

            # Get the filling axes --> without any masking.
            # The flattening has to be applied as the last step since the categories and subsamples
            # work at event level

            fill_categorical = {}
            fill_numeric = {}
            data_ndim = None

            
            # Custom part
            ###########################################################
            ###########################################################
            # We need now to iterate on categories and subsamples
            # Mask the events, the weights and then flatten and remove the None correctly
            for category, cat_mask in categories.get_masks():
                # loop directly on subsamples
                for subsample, subs_mask in subsamples.get_masks():
                    # logging.info(f"\t\tcategory {category}, subsample {subsample}")
                    mask = cat_mask & subs_mask
                    # Skip empty categories and subsamples
                    if ak.sum(mask) == 0:
                        continue

                    new_events = events[mask]
                    if category == "PLJ":
                        new_events["PhotonGood"] = new_events["PhotonPLJ"]
                        new_events["BJetGood"] = new_events["BJetGoodPLJ"]
                        new_events["JetGood"] = new_events["JetGoodPLJ"]
                    if category == "CRB":
                        new_events["PhotonGood"] = new_events["PhotonCRB"]
                        new_events["BJetGood"] = new_events["BJetGoodCRB"]
                        new_events["JetGood"] = new_events["JetGoodCRB"]
                    if category == "CRC":
                        new_events["PhotonGood"] = new_events["PhotonCRC"]
                        new_events["BJetGood"] = new_events["BJetGoodCRC"]
                        new_events["JetGood"] = new_events["JetGoodCRC"]
                    if category == "CRD":
                        new_events["PhotonGood"] = new_events["PhotonCRD"]
                        new_events["BJetGood"] = new_events["BJetGoodCRD"]
                        new_events["JetGood"] = new_events["JetGoodCRD"]                        
                    new_events["top"] = new_events.LeptonGood + new_events.BJetGood + new_events.neutrino
                    new_events["VLT"] = new_events.LeptonGood + new_events.BJetGood + new_events.neutrino + new_events.PhotonGood

                    for ax in histo.axes:
                        # Checkout the collection type
                        if ax.type in ["regular", "variable", "int"]:
                            if ax.coll == "events":
                                # These are event level information
                                data = new_events[ax.field]
                            elif ax.coll == "metadata":
                                data = new_events.metadata[ax.field]
                            elif ax.coll == "custom":
                                # taking the data from the custom_fields argument
                                # IT MUST be a per-event number, so we expect an array to mask
                                data = custom_fields[ax.field]
                            else:
                                if ax.coll not in new_events.fields:
                                    raise ValueError(
                                        f"Collection {ax.coll} not found in events!"
                                    )
                                # General collections
                                if ax.pos == None:
                                    data = getattr(new_events[ax.coll], ax.field)
                                elif ax.pos >= 0:
                                    data = ak.pad_none(
                                        new_events[ax.coll][ax.field], ax.pos + 1, axis=1
                                    )[:, ax.pos]
                                else:
                                    raise Exception(
                                        f"Invalid position {ax.pos} requested for collection {ax.coll}"
                                    )
        
                            # Flattening
                            if data_ndim == None:
                                data_ndim = data.ndim
                            elif data_ndim != data.ndim:
                                raise Exception(
                                    f"Incompatible shapes for Axis {ax} of hist {histo}"
                                )
                            # If we have multidim data we need to flatten it
                            # but we will do it after the event masking of each category
        
                            # Filling the numerical axes
                            fill_numeric[ax.name] = data
        
                        #### --> end of numeric axes
                        # Categorical axes (not appling the mask)
                        else:
                            if ax.coll == "metadata":
                                data = new_events.metadata[ax.field]
                                fill_categorical[ax.name] = data
                            elif ax.coll == "custom":
                                # taking the data from the custom_fields argument
                                data = custom_fields[ax.field]
                                fill_categorical[ax.name] = data
                            else:
                                raise NotImplementedError()
            
            ############################################################
            ############################################################
                    
                    # Check if the required data is dim=1, per event,
                    # and the mask is by collection.
                    # In this case the mask is reduced to per-event mask
                    # doing a logical OR only if explicitely allowed by the user
                    # WARNING!! POTENTIAL PROBLEMATIC BEHAVIOUR
                    # The user must be aware of the behavior.

                    if data_ndim == 1 and mask.ndim > 1:
                        if histo.collapse_2D_masks:
                            if histo.collapse_2D_masks_mode == "OR":
                                mask = ak.any(mask, axis=1)
                            elif histo.collapse_2D_masks_mode == "AND":
                                mask = ak.all(mask, axis=1)
                            else:
                                raise Exception(
                                    "You want to collapse the 2D masks on 1D data but the `collapse_2D_masks_mode` is not 'AND/OR'"
                                )

                        else:
                            raise Exception(
                                "+++++ BE AWARE! This is a possible mis-behavior! +++++\n"
                                + f"You are trying to fill the histogram {name} with data of dimention 1 (variable by event)"
                                + "and masking it with a mask with more than 1 dimension (e.g. mask on Jets)\n"
                                + "This means that you are either performing a cut on a collections (e.g Jets),"
                                + " or you are using subsamples with cuts on collections.\n"
                                + "\n As an example of a bad behaviour would be saving the pos=1 of a collection e.g. `JetGood.pt[1]`\n"
                                + "while also having a 2D cut on the `JetGood` collection --> this is not giving you the second jet passing the cut!\n"
                                + "In that case the 2nd JetGood.pt will be always plotted even if masked by the 2D cut: in fact "
                                + "the 2D masks would be collapsed to the event dimension. \n\n"
                                + "If you really wish to save the histogram with a single value for event (data dim=1)"
                                + "you can do so by configuring the histogram with `collapse_2D_masks=True\n"
                                + "The 2D masks will be collapsed on the event dimension (axis=1) doing an OR (default) or an AND\n"
                                + "You can configure this behaviour with `collapse_2D_masks_mode='OR'/'AND'` in the histo configuration."
                            )

                    # Mask the variables and flatten them
                    # save the isnotnone and datastructure
                    # to be able to broadcast the weight
                    has_none_mask = False
                    all_axes_isnotnone = None
                    has_data_structure = False
                    data_structure = None
                    fill_numeric_masked = {}
                    # loop on the cached numerical filling
                    for field, data in fill_numeric.items():
                        # Custom part
                        ###########################################################
                        ###########################################################
                        # Now mask and data dimension doesn't match. data is mask by default in new fill() structure.
                        masked_data = data
                        ###########################################################
                        ###########################################################
                        # For each field we need to mask and flatten
                        if data_ndim > 1:
                            # We need to flatten and
                            # save the data structure for weights propagation
                            if not has_data_structure:
                                data_structure = ak.ones_like(masked_data)
                                has_data_structure = True
                            # flatten the data in one dimension
                            masked_data = ak.flatten(masked_data)

                        # check isnotnone AFTER the flattening
                        if not has_none_mask:  # this is the first axis analyzed
                            all_axes_isnotnone = ~ak.is_none(masked_data)
                            has_none_mask = True
                        else:
                            all_axes_isnotnone = all_axes_isnotnone & (
                                ~ak.is_none(masked_data)
                            )
                        # Save the data for the filling
                        fill_numeric_masked[field] = masked_data

                    # Now apply the isnone mask to all the numeric fields already masked
                    for key, value in fill_numeric_masked.items():
                        # we also convert it to numpy to speedup the hist filling
                        fill_numeric_masked[key] = ak.to_numpy(
                            value[all_axes_isnotnone], allow_missing=False
                        )

                    # Ok, now we have all the numerical axes with
                    # data that has been masked, flattened
                    # removed the none value --> now we need weights for each variation
                    if not histo.no_weights and self.isMC:
                        if shape_variation == "nominal":
                            # if we are working on nominal we fill all the weights variations
                            for variation in histo.hist_obj.axes["variation"]:
                                if variation in self.available_shape_variations:
                                    # We ignore other shape variations when
                                    # we are already working on a shape variation
                                    continue
                                # Only weights variations, since we are working on nominal sample
                                # Check if this variation exists for this category
                                if variation not in weights[category]:
                                    # it means that the variation is in the axes only
                                    # because it is requested for another category
                                    # In this case we fill with the nominal variation
                                    # We get the weights for the current category
                                    weight_varied = weights[category]["nominal"]
                                else:
                                    # We get the weights for the current category
                                    weight_varied = weights[category][variation]

                                # Broadcast and mask the weight (using the cached value if possible)
                                weight_varied = self.mask_and_broadcast_weight(
                                    category,
                                    subsample,
                                    variation,
                                    weight_varied,
                                    mask,
                                    data_structure,
                                )
                                if custom_weight != None and name in custom_weight:
                                    weight_varied = weight_varied * self.mask_and_broadcast_weight(
                                        category + "customW",
                                        subsample,
                                        variation,
                                        custom_weight[
                                            name
                                        ],  # passing the custom weight to be masked and broadcasted
                                        mask,
                                        data_structure,
                                    )

                                # Then we apply the notnone mask
                                weight_varied = weight_varied[ak.to_numpy(all_axes_isnotnone)]
                                # Fill the histogram
                                try:
                                    self.histograms[subsample][name].hist_obj.fill(
                                        cat=category,
                                        variation=variation,
                                        weight=weight_varied,
                                        **{**fill_categorical, **fill_numeric_masked},
                                    )
                                except Exception as e:
                                    raise Exception(
                                        f"Cannot fill histogram: {name}, {histo} {e}"
                                    )
                        else:
                            # Check if this shape variation is requested for this category
                            if shape_variation not in self.available_shape_variations_bycat[category]:
                                # it means that the variation is in the axes only
                                # because it is requested for another category.
                                # We cannot fill just with the nominal, because we are running the shape
                                # variation and the observable hist will be different, also if with nominal weights.
                                continue
                                
                            # Working on shape variation! only nominal weights
                            # (also using the cache which is cleaned for each shape variation
                            # at the beginning of the function)
                            weights_nom = self.mask_and_broadcast_weight(
                                category,
                                subsample,
                                "nominal",
                                weights[category]["nominal"],
                                mask,
                                data_structure,
                            )
                            if custom_weight != None and name in custom_weight:
                                weight_nom = weight_nom * self.mask_and_broadcast_weight(
                                    category + "customW",
                                    subsample,
                                    "nominal",
                                    custom_weight[
                                        name
                                    ],  # passing the custom weight to be masked and broadcasted
                                    mask,
                                    data_structure,
                                )
                            # Then we apply the notnone mask
                            weights_nom = weights_nom[all_axes_isnotnone]
                            # Fill the histogram
                            try:
                                self.histograms[subsample][name].hist_obj.fill(
                                    cat=category,
                                    variation=shape_variation,
                                    weight=weights_nom,
                                    **{**fill_categorical, **fill_numeric_masked},
                                )
                            except Exception as e:
                                raise Exception(
                                    f"Cannot fill histogram: {name}, {histo} {e}"
                                )
                    ##################################################################################
                    elif not histo.no_weights and not self.isMC:   #DATA
                        # Broadcast and mask the weight (using the cached value if possible)
                        weight_data = weights[category]["nominal"]
                        weight_data = self.mask_and_broadcast_weight(
                            category,
                            subsample,
                            "nominal",
                            weight_data,
                            mask,
                            data_structure,
                        )
                        if custom_weight != None and name in custom_weight:
                            weight_data = weight_data * self.mask_and_broadcast_weight(
                                category + "customW",
                                subsample,
                                "nominal",
                                custom_weight[
                                    name
                                ],  # passing the custom weight to be masked and broadcasted
                                mask,
                                data_structure,
                            )

                        # Custom part
                        ###########################################################
                        ###########################################################
                        if category == "PLJ":
                            EF = ExtrapolationFactor(events)
                            EF_weight = EF.compute_EF(self.year)
                            weight_data = EF_weight * weight_data
                        ###########################################################
                        ###########################################################

                        # Then we apply the notnone mask
                        weight_data = weight_data[ak.to_numpy(all_axes_isnotnone)]
                        # Fill the histogram
                        try:
                            # Data histograms don't have variations but now can be weighted
                            self.histograms[subsample][name].hist_obj.fill(
                                cat=category,
                                weight=weight_data,
                                **{**fill_categorical, **fill_numeric_masked},
                            )
                        except Exception as e:
                            raise Exception(
                                f"Cannot fill histogram for Data: {name}, {histo} {e}"
                            )

                    ######################################################
                    elif (
                        histo.no_weights and self.isMC
                    ):  # NO Weights modifier for the histogram
                        try:
                            self.histograms[subsample][name].hist_obj.fill(
                                cat=category,
                                variation="nominal",
                                **{**fill_categorical, **fill_numeric_masked},
                            )
                        except Exception as e:
                            raise Exception(
                                f"Cannot fill histogram: {name}, {histo} {e}"
                            )

                    elif histo.no_weights and not self.isMC:
                        # Fill histograms for Data
                        try:
                            self.histograms[subsample][name].hist_obj.fill(
                                cat=category,
                                **{**fill_categorical, **fill_numeric_masked},
                            )
                        except Exception as e:
                            raise Exception(
                                f"Cannot fill histogram: {name}, {histo} {e}"
                            )
                    else:
                        raise Exception(
                            f"Cannot fill histogram: {name}, {histo}, not implemented combination of options"
                        )
