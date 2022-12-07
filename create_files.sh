#!/bin/bash

# 0.1 %
bash simulate_contamination.sh -p 0.1 -r 80 | tee log.cont01_r80M
mkdir /analysis/contamination/simulation/run_01
cd /analysis/contamination/simulation/run_01
bash ../run_ngs.sh cont0.1_reads80000000

# 1 %
#bash simulate_contamination.sh -p 1 -r 80 | tee log.cont1_r80M
#mkdir -p /analysis/contamination/simulation/run_1
#cd /analysis/contamination/simulation/run_1
#bash ../run_ngs.sh cont1_reads80000000

# 5 %
#bash simulate_contamination.sh -p 5 -r 80
#mkdir /analysis/contamination/simulation/run_5
#cd /analysis/contamination/simulation/run_5
#bash ../run_ngs.sh cont5_reads80000000

# 0.2 %
#bash simulate_contamination.sh -p 0.2 -r 80 | tee log.cont02_r80M
#mkdir /analysis/contamination/simulation/run_02
#cd /analysis/contamination/simulation/run_02
#bash ../run_ngs.sh cont0.2_reads80000000

# 10 %
#bash simulate_contamination.sh -p 10 -r 80
#mkdir /analysis/contamination/simulation/run_10
#cd /analysis/contamination/simulation/run_10
#bash ../run_ngs.sh cont10_reads80000000
