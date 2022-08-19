#!bin/sh

source ~/.bashrc
source /opt/intel/bin/compilervars.sh intel64


# Taylor-Green
./test.ex -k 2 -cs 0.5 -mf 1 -bn 20 -f ./data/0.05/Compute_Model_For_TriAngels.txt 2>&1 |tee -a run-Taylor-Green-h005.log
# ./test.ex -k 2 -cs 0.5 -mf $mf -bn $bn -f ./data/0.025/Compute_Model_For_TriAngels.txt 2>&1 |tee -a run-Taylor-Green-h0025.log
# ./test.ex -k 2 -cs 0.5 -mf $mf -bn $bn -f ./data/0.0125/Compute_Model_For_TriAngels.txt 2>&1 |tee -a run-Taylor-Green-h00125.log


sh move.sh
