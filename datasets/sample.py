import subprocess
import re
from xml.etree.ElementTree import ElementTree, Element, SubElement, Comment, tostring

class SampleNameParser:
    def __init__(self):
        self.regexps = []
        self.regexps.append( re.compile( r'^/(?P<sample>[^/_]*)_.*'))
        self.regexps.append( re.compile( r'.*_Tune(?P<tune>[^_/]*)[_/].*' ) )
        self.regexps.append( re.compile( r'.*/RunII(?P<campaign>[^-]*(?P<year>[0-9][0-9])).*' ) )
        self.regexps.append( re.compile( r'.*/(?P<isData>Run20)(?P<year>[0-9][0-9])(?P<era>.)[-_].*' ) )
        self.regexps.append( re.compile( r'.*-Nano(?P<nanotag>[^-_]*)[-_].*' ) )
        self.regexps.append( re.compile( r'.*NanoAODv(?P<nanoversion>[0-9]*).*' ))
        self.regexps.append( re.compile( r'.*NanoAODAPVv(?P<nanoversion>[0-9]*).*' ))
        self.regexps.append( re.compile( r'.*(?P<prevfp>NanoAODAPV).*' )  )
        self.regexps.append( re.compile( r'.*_ext(?P<ext>[0-9]*).*' )  )
        self.regexps.append( re.compile( r'.*-v(?P<version>[0-9]*).*' )  )
        self.regexps.append( re.compile( r'.*-ver(?P<ver>[0-9]*).*' )  )
        self.regexps.append( re.compile( r'.*(?P<madgraph>madgraph).*' )  )
        self.regexps.append( re.compile( r'.*(?P<herwig>herwig).*' )  )
        self.regexps.append( re.compile( r'.*(?P<amcatnlo>amcatnlo).*' )  )
        self.regexps.append( re.compile( r'.*(?P<pythia>pythia).*' )  )
        self.regexps.append( re.compile( r'.*(?P<sherpa>[sS]herpa).*' )  )
        self.regexps.append( re.compile( r'.*(?P<pmx>new_pmx).*' )  )
        self.regexps.append( re.compile( r'.*(?P<backup>backup).*' )  )
        self.regexps.append( re.compile( r'.*(?P<ps>[^_]*PS[^_/]*).*' )  )
        self.regexps.append( re.compile( r'.*(?P<hipm>HIPM).*' )  )
        self.regexps.append( re.compile( r'/ST_(?P<ST>(t-channel|tW|s-channel)_(top|antitop|4f)).*' )  )
        self.regexps.append( re.compile( r'/TTGamma_(?P<TTG>(Dilept|Hadronic|SingleLept)).*' )  )
        self.regexps.append( re.compile( r'/TTTo(?P<TT>(2L2Nu|SemiLeptonic|Hadronic)).*' )  )
        self.regexps.append( re.compile( r'/GJets_HT-(?P<GJets>40To100|100To200|200To400|400To600|600ToInf).*' )  )
        self.regexps.append( re.compile( r'/WGJets_MonoPhoton_PtG-(?P<WGJets>40to130|130)_TuneCP5.*' )  )
        self.regexps.append( re.compile( r'/DYJetsToLL_M-(?P<DYJets>50|10to50)_.*TuneCP5.*' ) )
        # self.regexps.append( re.compile( r'/W(?P<WJets>[1-4])JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.*' ) )
        self.regexps.append( re.compile( r'/VLTquarkToPhotonTop_M-(?P<Signal>600|700|800|900|1000|1100|1200|1300|1400|1500|1600|1700|1800|1900|2000)_PtG-10_TuneCP5_13TeV-madgraph-pythia8/.*' ) )

    def parse(self , sample):
        ret = [ re.match( sample ) for re in self.regexps ]
        info = {}
        for regexp in [ r for r in ret if r ]:
            for group, value in regexp.groupdict().items():
                info[group] = value
        print("info", info)
        return ret , info
    
class Sample:
    def __init__(self , ds , parser = None):
        self.ds = ds
        self.parser = parser if parser else SampleNameParser()
        self.regexp_results , self.info = self.parser.parse( ds )
        self.sample_dict ={
            "ST": "ST",
            "TTTo2L2Nu": "TT",
            "TTToHadronic": "TT",
            "TTToSemiLeptonic": "TT",
            "TTGamma": "TTG",
            "GJets": "GJets",
            "TGJets": "TGJets",
            "WJetsToLNu" : "WJets",
            "DYJetsToLL" : "DYJets",
            "WGToLNuG" : "WG",
            "ZGToLLG": "ZG",
            "WWG": "WWG",
            "WZG" : "WZG",
            "ZZGTo4L": "ZZG",
            "WW": "WW",
            "ZZ": "ZZ",
            "WZ": "WZ",
            # "WGammaToJJGamma": "WG",
            # "WGJets": "WGJets",
            # "ZGammaToJJGamma": "ZG",
            # "DiBoson" : "diboson",
            # "WZTo3LNu" : "WZ"
        }

    def sample_name(self):
        sample = self.info.get("sample", "data")
        if sample == "data":
            smpl_name = self.ds.split("/")[1]
        elif sample == "VLTquarkToPhotonTop":
            smpl_name = "Signal_" + self.info.get("Signal")
            print("smpl_name", smpl_name)
        else:
            smpl_name = self.sample_dict.get(sample)
            print("smpl_name", smpl_name)
        return smpl_name

    def prevfp(self):
        if "hipm" in self.info.keys() or "prevfp" in self.info.keys():
            return True
        else:
            return False

    def isData(self):
        return self.info.get('isData' , False) != False

    def year(self):
        return 2000 + int( self.info['year'] )
    
    def tune(self):
        return self.info.get('tune' , 'noTune' )

    def makeUniqueName(self):
        uName = self.ds.split('/')[1]
        if self.isData():
            uName = "{0}_{1}_{2}".format( uName , self.year() , self.info.get('era') )
            if "ver" in self.info:
                uName += "_ver" + self.info["ver"]
            if "hipm" in self.info:
                uName += "_" + self.info["hipm"]
        else:
            uName = "{0}_{1}".format( uName , self.year() )
            if 'pmx' in self.info :
                uName += "_pmx"
            if 'ext' in self.info :
                uName += '_ext{0}'.format( self.info['ext'] )
            if 'backup' in self.info :
                uName += '_backup'
            if 'ps' in self.info:
                uName += '_' + self.info['ps']
            if 'prevfp' in self.info:
                uName += '_preVFP'
        if 'version' in self.info:
            uName += '_v' + self.info['version']

        return uName

    def get_dataset_info(self):
        dataset_info = {}
        dataset_info["dataset"] = self.ds
        dataset_info["sample"] = self.sample_name()
        dataset_info["isData"] = self.isData()
        dataset_info["era"] = self.info.get("era", "isMC")
        dataset_info["part"] = self.info.get(self.sample_name(), None)
        if self.year() == 2016:
            if "hipm" in self.info.keys() or "prevfp" in self.info.keys():
                dataset_info["year"] = str(self.year()) + "_PreVFP"
            else:
                dataset_info["year"] = str(self.year()) + "_PostVFP"
        else:
            dataset_info["year"] = str(self.year())
        return dataset_info

    def GetParent(self):
        if not hasattr( self , 'parents' ):
            process = subprocess.Popen( [ '/cvmfs/cms.cern.ch/common/dasgoclient', "-query=parent dataset={0}".format(self.ds) ], stdout=subprocess.PIPE)
            self.parents = []
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    s = output.strip()
                    self.parents.append( s )
        if len(self.parents) == 1:
            return Sample( self.parents[0] )
        elif len(self.parents) > 1:
            print( 'here is the list of found parents for {0} : {1}'.format( self.ds , self.parents ) )
            # Filter parents that contain 'MINIAOD'                   
            return Sample(self.parents[0])
        else:
            print("No parents found for {0}".format (self.ds))
            return None

class SampleList:
    def __init__(self, dataset_file):
        self.dataset_file = dataset_file
        self.dataset_list = self.get_all_datasets()
        self.unique_names_dict = self.match_unique_name_dataset()
        
    def get_all_datasets(self):
        with open(self.dataset_file) as f:
            return [a[:-1] for a in f.readlines() if a.strip() and not a.startswith('#')]
    
    def match_unique_name_dataset(self):
        unique_name_dict = {}
        for ds in self.dataset_list:
            smpl = Sample(ds)
            unique_name = smpl.makeUniqueName()
            unique_name_dict[unique_name] = ds
        return unique_name_dict

    def match_sample_dataset(self):
        pass

    def get_dataset_info(self, unique_name):
        dataset = self.unique_names_dict[unique_name]
        smpl = Sample(dataset)
        dataset_info = smpl.get_dataset_info()
        return dataset_info

    def get_sample_info(self, sample_name):
        pass
