#!bin/sh

if [ ! -d "results" ];
then
    mkdir results
else
    rm -rf results
    mkdir results
fi


mv Physical_Quantity_Statistics.txt results
mv -f TecPlot results

