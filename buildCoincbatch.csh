#!/bin/bash
#clear

#export workDir=/sps/edelweis/rootDataRun317/streams/Selection10keV_Migdal/DEBBUG_DIR_2
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
 cp -r $SimuData_Dir/Simu/$run/$detector/Line_* Simulated_data/Simu/.
 pwd
 ls -lart 
 cd ..
 pwd
 ls -lart 
 ####
 echo "End copying"
#fi
#################<-----   Fin Non fabrication fichier tar

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
make
echo " //ooooooooo Apres le make "
ls -lart
PROCESSED_DATA="../Processing/ProcessedData"
	
echo "-----------> done !"
ls -lrt
echo " "
echo " "


pwd
echo "./BuildSimuCoincArgs_Batch.exe $run"
./BuildSimuCoincArgs_Batch.exe $run

echo "... Saving output files ... #####################################"
echo "Copying files ..."
echo "On commence les cp"


mkdir -p $OUT_DIR/$run
mkdir -p $OUT_DIR/$run/$detector
   #mkdir -p $OUT_DIR/$run/$detector/StreamOF
echo "cp -r Simulated_data/SimuCoinc/* $OUT_DIR/."
cp -r Simulated_data/SimuCoinc/Line_* $OUT_DIR/$run/$detector/.

chmod -R g+rw $OUT_DIR/$run
chmod -R g+rw $OUT_DIR/$run

exit
 
