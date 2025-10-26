ulimit -n 4098
pocket-coffea run --cfg ele_config.py --executor futures -o output --ignore-grid-certificate

pocket-coffea run --cfg ele_config.py -o output_config -e condor@lxplus --scaleout 100 --chunksize 200000 --job-dir jobs --job-name skim --queue workday --local-virtualenv

pocket-coffea run --cfg example_config.py --test  -o output_test
