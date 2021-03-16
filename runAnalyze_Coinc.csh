#!/bin/bash
#clear

#export workDir=/sps/edelweis/rootDataRun317/streams/Selection10keV_Migdal/DEBBUG_DIR_3
#mkdir  $workDir
#cd $workDir
#pwd

 mkdir  subData
 mkdir  subData/$run
 mkdir  Processing
 cp $SOFT_DIR/*.cxx  Processing/.
 cp $SOFT_DIR/*.csh  Processing/.
 cp $SOFT_DIR/Makefile  Processing/.
 cd Processing
 mkdir  ProcessedData
 mkdir  ProcessedData/$run
 mkdir  ProcessedData/$run/$detector
 mkdir  ProcessedData/$run/$detector/CalibratedData
 cp -r $ProcessedData_Dir/$run/$detector/*.root ProcessedData/$run/$detector/.  ## attention il y a le directory RootFiles 
 cp -r $ProcessedData_Dir/$run/$detector/CalibratedData/*.root ProcessedData/$run/$detector/CalibratedData/.  ## attention il y a le directory RootFiles
 mkdir  Simulated_data
 mkdir  Simulated_data/Simu
 mkdir  Simulated_data/SimuCoinc
 #cp -r $SimuData_Dir/Simu/$run/$detector/Line_* Simulated_data/Simu/.#!!!!!!!!!!!!!
 cp -r $SimuData_Dir/SimuCoinc/$run/$detector/Line_* Simulated_data/SimuCoinc/.
 mkdir Figures
 mkdir Bank_Lines
 pwd
 ls -lart 
 cd ..
 pwd
 ls -lart 
 ####
 echo "End copying"
 
echo "rootsys avant ccenv: "$ROOTSYS
ccenv root 6.08.02
#export ROOTSYS /pbs/software/centos-7-x86_64/root/6.08.02/
echo "rootsys apres ccenv: "$ROOTSYS
pwd
cd Processing
pwd
chmod 775 *
echo " //ooooooooo Avant le make "
ls -lart
make clean
make
echo " //ooooooooo Apres le make "
ls -lart
PROCESSED_DATA="../Processing/ProcessedData"
	
echo "-----------> done !"
ls -lrt
echo " "
echo " "
pwd
echo "./Analyze_Coinc_cleaned_BATCH.exe $run "
./Analyze_Coinc_cleaned_BATCH.exe $run $Tension 1 $Ion_cut

echo "... Saving output files ... #####################################"
echo "Copying files ..."
echo "On commence les cp"


mkdir -p $SOFT_DIR/Efficiencies_batch_return
   #mkdir -p $OUT_DIR/$run/$detector/StreamOF
echo "cp -r *.root $SOFT_DIR/."
cp -r Efficiencies_*.root $SOFT_DIR/Efficiencies_batch_return/.
cp Figures/* $SOFT_DIR/Figures/.
chmod -R g+rw $SOFT_DIR/Efficiencies_batch_return
chmod -R g+rw $SOFT_DIR/Efficiencies_batch_return
chmod -R g+rw $SOFT_DIR/Figures

exit


