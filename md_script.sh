#!/bin/bash

####   It is a script for for XO-SINR thermostat  ####
####   Written by : Ritama Kar ####
MTS="XO"
del_t=0.1

#Remove the previous output folder
rm -rf out_new

#Create a new folder output
mkdir out_new

#enter to that directory
cd out_new

#create another directory to store probability distribution files
mkdir outnew_prob
cd ..

#compile the program
make
echo "MTS method is=$MTS"
echo "time        Kintetic energy" >> energy_${del_t}.dat

#Run the program for different time step
#for del_t in `seq 0.2 0.05 1.2` ;do
    
    echo "Time step size=$del_t"

#copy the template input file to input file
    cp sinr_template.in sinr.in

#Search XR and replace MTS method in the input file
    sed -i "s/XR/$MTS/g" sinr.in

#search XXX and replace with delta t
    sed -i "s/XXX/$del_t/g" sinr.in

#execute the program
    ./sinr.exe >> energy_${del_t}.dat

#move the trajectory output to output folder
    mv Trajectory_XO_atom_1.out out_new/Trajectory_XO_${del_t}.out

#move the conserved output file to output folder
    mv Conserved_XO.out out_new/Conserved_XO_${del_t}.out
    mv XYZ_XO.out out_new/XYZ_XO_${del_t}.out

#enter to the output directory
    cd out_new

#copy the output trajectory file to new input file
    cp Trajectory_XO_${del_t}.out Trajectory.in

#This section is about to plot histogram	

    ../hist.exe

#Copy the Probability output file to Probaility input file to plot against every time step
    cp Probability_x.out Probability_x.in
    ../plot.exe >> gamma_${del_t}.out
   # lg=`grep "gamma" gamma_${del_t}.out | awk '{print $3}'`
    
   # echo $del_t $lg >> plot.dat

#Move the output probability file to outnew_probability folder
    mv Probability_x.out outnew_prob/Probability_x_${del_t}_${MTS}.out


#back to main directory
    cd ..

echo "Complete for $del_t"

#done

echo "Job done"
