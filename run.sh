#!bin/sh

source ~/.bashrc
source /opt/intel/bin/compilervars.sh intel64

############ argparser  ###########
until [ $# -eq 0 ]
do
  name=${1:1}; shift;
  if [[ -z "$1" || $1 == -* ]] ; then eval "export $name=true"; else eval "export $name=$1"; shift; fi
done
########## set default  ###########
if [ ${#nt} == 0 ]
	then
	    nt=1;  # default value to give    
else
	    nt=$nt;
fi
########## set default  ###########
if [ ${#mf} == 0 ]
	then
	    mf=1;  # default value to give    
else
	    mf=$mf;
fi
########## set default  ###########
if [ ${#bn} == 0 ]
	then
	    bn=20;  # default value to give    
else
	    bn=$bn;
fi
########## set default  ###########
if [ ${#k} == 0 ]
	then
	    k=2;  # default value to give    
else
	    k=$k;
fi



# Taylor-Green
./test.ex -k 2 -cs 0.5 -mf $mf -mfm 0 -bn $bn -f ./data/0.05/Compute_Model_For_TriAngels.txt 2>&1 |tee -a run-Taylor-Green-h005.log
# ./test.ex -k 2 -cs 0.5 -mf $mf -mfm 0 -bn $bn -f ./data/0.025/Compute_Model_For_TriAngels.txt 2>&1 |tee -a run-Taylor-Green-h0025.log
# ./test.ex -k 2 -cs 0.5 -mf $mf -mfm 0 -bn $bn -f ./data/0.0125/Compute_Model_For_TriAngels.txt 2>&1 |tee -a run-Taylor-Green-h00125.log


sh move.sh
