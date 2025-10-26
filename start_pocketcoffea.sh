apptainer shell  -B /afs -B /cvmfs/ -B /tmp  -B /eos/cms/ -B /eos/user/m/msahraei/ -B /etc/sysconfig/ngbauth-submit -B ${XDG_RUNTIME_DIR}  --env KRB5CCNAME="FILE:${XDG_RUNTIME_DIR}/krb5cc" /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-el9-stable

# Activate the virtual environment
cd ../PocketCoffea
source myenv/bin/activate
export PYTHONPATH=`pwd`
cd -
