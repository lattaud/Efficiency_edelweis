#!/bin/bash

echo "---------------------------------"
date
echo "---------------------------------"

#
#  Do calibration on CCA
#




STUDYDIRNAME="FID848-event-bank" # repertoire simu



CCDIR="/sps/edelweis/rootDataRun317/streams/Efficiency-FID848"
export CCDATADIR="$CCDIR/$STUDYDIRNAME"
export CALIBFILEDIR="/sps/edelweis/rootDataRun317/streams/prodk" # repertoire avec les .calib
export bigRunNb="317" # Run number

export CC_EVEREST_PROG_DIR="/pbs/home/h/hlattaud/Everest_from_git" # changer pour mon everest  
CC_EVEREST_SCRIPT_DIR="/sps/edelweis/rootDataRun317/streams/Selection10keV_Migdal" # endroit où tu mets everest simu Hugues

SimuType="Line" # pour line mono energetique

############### DETECTOR LIST :  detlist=("*") for All
detlist=("FID848") # list des detecteur
#detlist=("RED30" "NbSi209")
############### RUN LIST :  runlist=("*") for All
runlist=("tk*")
#runlist=("te03a000" "te04a001" "te05a000")

for (( d=0; d<${#detlist[@]}; d++ )); do 
    for (( r=0; r<${#runlist[@]}; r++ )); do 

        for toto in ` ls -1d $CCDATADIR/Simu/${runlist[$r]}/${detlist[$d]}/$SimuType* ` ; do
                export pathwheretocreatedir=$toto"/CalibratedData"
                if [ ! -d $pathwheretocreatedir ]; then
                echo "mkdir -p "$pathwheretocreatedir
                mkdir -p $pathwheretocreatedir
                fi                
          
        done

        for LineDirectory in ` ls -1d $CCDATADIR/Simu/${runlist[$r]}/${detlist[$d]}/$SimuType* | awk -F'_' 'BEGIN{OFS=FS}{NF--}1' | uniq` ; do

        export LineDirectory=$LineDirectory
        export runNb=`echo $LineDirectory | awk -F"/" '{print $(NF-2)}' `
        export det=`echo $LineDirectory | awk -F"/" '{print $(NF-1)}' `        

        echo "det=$det , run=$runNb"
        echo "LineDirectory : $LineDirectory"

        #sh $CC_EVEREST_SCRIPT_DIR/testsimuCC.sh
        LineNB=`echo $LineDirectory | awk -F"$det/" '{print $2}' `
        export jobNb="${STUDYDIRNAME}_${LineNB}_${runNb}_${det}"
        echo "do it $jobNb"
        
        qsub -v CCDATADIR,CALIBFILEDIR,bigRunNb,CC_EVEREST_PROG_DIR,LineDirectory,runNb,det,jobNb -N $jobNb -P P_edelweis -l vmem=4G -l sps=1 -j y -o /pbs/home/h/hlattaud/Everest_from_git/logfiles/job_${jobNb}.log $CC_EVEREST_SCRIPT_DIR/EVEREST_SIMU_HUGUES.sh
          
         
        done
       
     done
done 

exit


