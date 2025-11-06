import numpy as np
import math

from dataclasses import dataclass, field
from typing import Dict, List, Optional

@dataclass
class HistInfo:
    var_name: str
    ax_label: str
    cat: str
    year: int
    bins: List = field(default_factory=list)
    backgrounds: Dict[str, List] = field(default_factory=dict)
    data: List = field(default_factory=list)
    signals: Dict[str, List] = field(default_factory=dict)
    # syst: List = field(default_factory=list)
    has_plot: bool = False 

    @property
    def stat(self) -> List[float]:
        """Compute statistical errors dynamically."""
        return [math.sqrt(abs(x)) for x in self.data]

    # Temporary syst
    @property
    def syst(self) -> List[float]:
        """Compute syst errors dynamically."""
        return [math.sqrt(abs(x)) for x in self.total_backgrounds]

    @property
    def total_backgrounds(self) -> List[float]:
        return sum(list(self.backgrounds.values()))

    @property
    def bin_centers(self) -> List:
        centers = []
        for i in range(len(self.bins)-1):
            center = (self.bins[i] + self.bins[i+1]) / 2
            centers.append(center)
        return centers

    @property
    def bin_widths(self) -> List:
        return np.diff(self.bins)


class ExtractHistData:
    """
    This class handles coffea output reading and creating HistInfo for HistogramPlotter
    """
    def __init__(self, coffea_output, year):
        self.output = coffea_output
        self.year = year
        self.hist_info = {
            "year": year,
            "backgrounds": {},
            "signals": {},
            "has_plot": False
        }

    def _add_bins_info(self, var:str):
        smpl = list(self.output["variables"][var].keys())[0]
        hist = list(self.output["variables"][var][smpl].values())[0]
        if "DATA_" in smpl:
            self.hist_info["bins"] = hist.axes[1].edges
            self.hist_info["ax_label"] = hist.axes[1].label
        else:
            self.hist_info["bins"] = hist.axes[2].edges
            self.hist_info["ax_label"] = hist.axes[2].label

    def _add_data_info(self, var, cat, data):
        data_hist = None
        for smpl in data[self.year]:
            for hist in self.output["variables"][var][smpl].values():
                if data_hist is not None:
                    data_hist += hist
                else:
                    data_hist = hist
        self.hist_info["data"] = (data_hist[{"cat": cat[-3:]}] + data_hist[{"cat": cat[:-4]}]).values()
        self.hist_info["has_plot"] = True

    def _add_MC_BCs_info(self, var, cat, MC_BCs):
        for bc, bc_lst in MC_BCs.items():
            mc_hist = None
            for smpl in bc_lst:
                for hist in self.output["variables"][var][smpl].values():
                    if mc_hist is not None:
                        mc_hist += hist
                    else:
                        mc_hist = hist
            self.hist_info["backgrounds"][bc] = (mc_hist[{"cat": cat[-3:], "variation": "nominal"}] +
                                                     mc_hist[{"cat": cat[:-4], "variation": "nominal"}]).values()
            self.hist_info["has_plot"] = True

    def _add_signals_info(self, var, cat, signals):
        for signal in signals:
            hist = self.output["variables"][var][signal][f"{signal}_{self.year}"]
            if cat in hist.axes[0]:
                self.hist_info["signals"][signal] = hist[{"cat": cat, "variation": "nominal"}].values()*0.15
                self.hist_info["has_plot"] = True

    def _add_jet_fake_photon_info(self, var, cat, dd_bcs):
        dd_hist = None
        for smpl in dd_bcs[self.year]:
            for hist in self.output["variables"][var][smpl].values():
                if dd_hist is not None:
                    dd_hist += hist
                else:
                    dd_hist = hist
        if "b0_" in cat:
            fake_smpl = "b0_PLJ"
        else:
            fake_smpl = "PLJ"
        if cat in dd_hist.axes[0]:
            self.hist_info["backgrounds"][r"$Jets\rightarrow\gamma$"] = dd_hist[{"cat": fake_smpl}].values()
            self.hist_info["has_plot"] = True

    def extract_hist_info(self, var, cat, data=None, MC_BCs=None, signals=None, **DD_BCs):
        self.hist_info["var_name"] = var
        self.hist_info["cat"] = cat
        self._add_bins_info(var)
        if data is not None:
            self._add_data_info(var, cat, data)
        if MC_BCs is not None:
            self._add_MC_BCs_info(var, cat, MC_BCs)
        if signals is not None:
            self._add_signals_info(var, cat, signals)
        if (r"$Jets\rightarrow\gamma$" in DD_BCs) and ("SR" in cat):
            self._add_jet_fake_photon_info(var, cat, DD_BCs[r"$Jets\rightarrow\gamma$"])

        return HistInfo(**self.hist_info)
    