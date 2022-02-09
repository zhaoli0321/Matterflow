#!bin/sh

if [ ! -d "results" ];
then
    mkdir results
else
    rm -rf results
    mkdir results
fi


mv Physical_Quantity_Statistics.txt results
mv Velocity_x_BottomEdge.txt results
mv Pressure_BottomEdge.txt results
mv *.log results
mv -f TecPlot results

