#!bin/sh

source ~/.bashrc
source /opt/intel/bin/compilervars.sh intel64


### Taylor-Green
# ./test.ex -tg 1 -k 0.5 -mf 1 -bn 20 -f ./data/Taylor-Green/0.05/Compute_Model_For_TriAngels.txt


### Gresho vortex
# without heat flux and without matter-flow
# ./test.ex -go 1 -k 0.1 -mf 0 -f ./data/Gresho/0.025/Compute_Model_For_TriAngels.txt 
# without heat flux and with matter-flow
./test.ex -go 1 -k 0.1 -mf 1 -f ./data/Gresho/0.025/Compute_Model_For_TriAngels.txt


### Noh
# with heat flux and with matter-flow
# ./test.ex -k 0.5 -td 0.1 -mf 1 -bn 20 -f ./data/Noh/32x32/Compute_Model_For_TriAngels.txt

# with heat flux and without matter-flow
# ./test.ex -k 0.5 -td 0.1 -mf 0 -bn 20 -f ./data/Noh/32x32/Compute_Model_For_TriAngels.txt

# without heat flux and without matter-flow
# ./test.ex -k 0.5 -td 0.0 -mf 0 -bn 20 -f ./data/Noh/32x32/Compute_Model_For_TriAngels.txt


### Sedvo
# ./test.ex -k 0.5 -td 0.05 -mf 1 -bn 20 -f ./data/Sedvo/32x32/Compute_Model_For_TriAngels.txt

sh move.sh
