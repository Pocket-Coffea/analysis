import os
from sample import *
import uproot 
import json
import subprocess

class FilesetCreator:
    def __init__(self, crab_output, dataset_file, xsections_file):
        if not crab_output.endswith("/"):
            self.base_dir = crab_output + "/"
        else:
            self.base_dir = crab_output
        self.sample_list = SampleList(dataset_file)
        self.xsections_file = xsections_file
        self.fileset = {}
    
    def extract_fileset(self, on_GRID=True):
        for root, _, _ in os.walk(self.base_dir):
            if root.endswith("0000"):
                unique_name = self.extract_unique_name(root)
                print(unique_name)
                dataset_info = self.sample_list.get_dataset_info(unique_name)
                if dataset_info["isData"]:
                    sample_name = "{}_{}_{}".format(dataset_info["sample"], dataset_info["year"], dataset_info["era"])
                elif dataset_info["part"] is not None:
                    sample_name = "{}_{}_{}".format(dataset_info["sample"], dataset_info["part"], dataset_info["year"])
                else:
                    sample_name = "{}_{}".format(dataset_info["sample"], dataset_info["year"])
                combined_file_path = self.hadd_files(root, target_file="{}.root".format(sample_name))
                nevents = self.n_events(combined_file_path)
                if sample_name not in self.fileset:
                    self.fileset[sample_name]={}
                    self.fileset[sample_name]["files"] = []
                    self.fileset[sample_name]["metadata"] = {}
                    self.fileset[sample_name]["metadata"]["das_name"] = dataset_info["dataset"]
                    self.fileset[sample_name]["metadata"]["sample"] = "DATA_" + dataset_info["sample"] if dataset_info["isData"] else dataset_info["sample"]
                    self.fileset[sample_name]["metadata"]["year"] = dataset_info["year"]
                    self.fileset[sample_name]["metadata"]["isMC"] = "False" if dataset_info["isData"] else "True"
                    self.fileset[sample_name]["metadata"]["nevents"] = nevents
                    if dataset_info["isData"]:
                        self.fileset[sample_name]["metadata"]["primaryDataset"] = dataset_info["sample"]
                        self.fileset[sample_name]["metadata"]["era"] = dataset_info["era"]
                    else:
                        xsec = self.xsection(dataset_info["dataset"])
                        self.fileset[sample_name]["metadata"]["xsec"] = xsec
                        self.fileset[sample_name]["metadata"]["part"] = dataset_info["part"] if dataset_info["part"] is not None else "nopart"
                        
                    if on_GRID:
                        combined_file_path = "root://eosuser.cern.ch/" + combined_file_path
                    self.fileset[sample_name]["files"].append(combined_file_path)
        return self.fileset
            
    def hadd_files(self, dir, target_file):
        if target_file not in os.listdir(dir):
            files = [os.path.join(dir,f) for f in os.listdir(dir) if ".root" in f]
            cmd = ["hadd", os.path.join(dir, target_file)] + files
            subprocess.run(cmd, check=True)
        return os.path.join(dir, target_file)
    
    def extract_unique_name(self, dir):
        search_string = "crab_"
        index = dir.find(search_string)
        name_index = index + len(search_string)
        unique_name = dir[name_index:].split("/")[0]
        return unique_name

    def file_has_entry(self, file_name):
        with uproot.open(file_name) as file:
            tree = file['Events']
            n_entries = tree.num_entries        
            if n_entries == 0:
                print("{} Has No Entry".format(file_name))
                return False
            else:  
                return True

    def n_events(self, file):
        with uproot.open(file) as f:
            cutflow_hist = f["cutflow"]
            nevents = cutflow_hist.values()[0]
        return nevents.item()

    def xsection(self, ds):
        with open(self.xsections_file) as file:
            xsec_dict = json.load(file)
        return xsec_dict[ds]

    def create_fileset_json(self, fileset_name="fileset.json", on_GRID=True):
        fileset = self.extract_fileset(on_GRID)
        try:  
            with open(fileset_name, 'w') as json_file:  
                json.dump(fileset, json_file, indent=4)
        except Exception as e:  
            print(f"An error occurred while saving to JSON: {e}")


if __name__ == "__main__":
    path_to_crab_output = "/eos/user/m/msahraei/PocketCoffea_input/"
    dataset_file = "Had_NanoAODv9.lst"
    xsec_file = "nano_xsec.json"
    fileset_creator = FilesetCreator(path_to_crab_output, dataset_file, xsec_file)
    fileset_creator.create_fileset_json("local_fileset.json", False)