#!bin/sh

source ~/.bashrc
source /opt/intel/bin/compilervars.sh intel64


### Taylor-Green
./test.ex -tg 1 -k 1 -cs 0.5 -mf 1 -bn 20 -f ./data/Taylor-Green/0.05/Compute_Model_For_TriAngels.txt
# ./test.ex -tg 1 -k 1 -cs 0.5 -mf 1 -bn 20 -f ./data/Taylor-Green/0.025/Compute_Model_For_TriAngels.txt
# ./test.ex -tg 1 -k 1 -cs 0.5 -mf 1 -bn 20 -f ./data/Taylor-Green/0.0125/Compute_Model_For_TriAngels.txt 


### Gresho vortex
# ./test.ex -go 1 -cs 0.5 -k 1 -td 0 -mf 1 -bn 20 -f ./data/Gresho/0.025/Compute_Model_For_TriAngels.txt 
# ./test.ex -go 1 -cs 0.1 -k 0.2 -td 0 -mf 1 -bn 20 -f ./data/Gresho/0.025/Compute_Model_For_TriAngels.txt 


### Noh
# with heat flux and with matter-flow
# ./test.ex -cs 0.5 -k 1 -td 0.1 -mf 1 -bn 20 -f ./data/Noh/32x32/Compute_Model_For_TriAngels.txt

# with heat flux and without matter-flow
# ./test.ex -cs 0.5 -k 1 -td 0.1 -mf 0 -bn 20 -f ./data/Noh/32x32/Compute_Model_For_TriAngels.txt

# without heat flux and without matter-flow
# ./test.ex -cs 0.5 -k 1 -td 0.0 -mf 0 -bn 20 -f ./data/Noh/32x32/Compute_Model_For_TriAngels.txt


### Sedvo
# ./test.ex -cs 0.1 -k 1 -td 0.05 -mf 1 -bn 20 -f ./data/Sedvo/32x32/Compute_Model_For_TriAngels.txt

sh move.sh
