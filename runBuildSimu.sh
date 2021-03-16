#!/bin/bash
j=0
 j=$[$j +1]
 export jobTyp=${!j}    
 j=$[$j +1]
 export DATA_DIR=${!j}    
 j=$[$j +1]
 export run=${!j}    
 j=$[$j +1]
 export detector=${!j}    
 j=$[$j +1]
 export trigchan=${!j}    
 j=$[$j +1]
 export SOFT_DIR=${!j}     
 j=$[$j +1]
 export OUT_DIR=${!j}    
 j=$[$j +1]
 export ProcessedData_Dir=${!j}    
 j=$[$j +1]
 export SimuData_Dir=${!j}  
    
 
#################---> Check OUTPUT Directory
 if [[ -d "$OUT_DIR" ]] ; then
    echo  "Be carefull output $OUT_DIR directory already exists  "    
  else  
    echo  "Creating  $OUT_DIR directory   "
    mkdir -p $OUT_DIR
    chmod -R g+rw $OUT_DIR
 fi
 
#################---> List of Variables
  echo  "	-------------    Inputs   -----------       "  
  echo  "	jobTyp=$jobTyp" 
  echo  "	DATA_DIR=$DATA_DIR"
  echo  "	run=$run"
  echo  "	Detector=$detector"
  echo  "	trigchan=$trigchan"
  echo  "	SOFT_DIR=$SOFT_DIR"
  echo  "	OUT_DIR=$OUT_DIR"
  echo  "	ProcessedData_Dir=$ProcessedData_Dir"
  echo  "   SimuData_Dir=$SimuData_Dir"

#################---> List of detectors
if [  $Detector == "All" ]; then
   Detectors=( "RED11" "RED30" "FID803" "FID824" "FID842" "FID847" "FID848" "NbSi206" "NbSi207" "NbSi208" "NbSi209" "CWO" )
elif [ $Detector == "All1" ]
then
# List of detectors with ChalB trigger
   Detectors=( "RED11" "FID803" "FID824" "FID842" "FID847" "FID848" "NbSi207" "NbSi208" "NbSi209" "CWO" )
else
   Detectors=("$4") 
fi
#################<--- End list of detectors


###################---> Fabrique le fichier tar des data
#################<-----   Fin fabrication fichier tar
 
cd $SOFT_DIR
pwd
#echo "Le fichier tar est $tarName.tar"

#################----> Boucle sur tout les detecteurs pour lancer le batch
if [ "$jobTyp" == "Batch" ]; then 
   for Det in "${Detectors[@]}"; do
         export detector=$Det
         export Time_SEED=`date +%s`
         export EXTNAME=$Time_SEED
         export jobName=$run"_"$detector"_Chan"$trigchan"_"$Simu_Type"_TimeSec_"$EXTNAME
     echo "EXTNAME"
     echo $EXTNAME 
     echo "jobName"
     echo $jobName     
	 mkdir -p  $OUT_DIR/logFiles
	 chmod -R g+rw $OUT_DIR/logFiles
	 rm -f $OUT_DIR/logFiles/$jobName.log
         qsub -v run,detector,trigchan  \
          -v SOFT_DIR,DATA_DIR,OUT_DIR,jobName  \
          -v ProcessedData_Dir,SimuData_Dir \
	      -N $jobName -P P_edelweis -l vmem=4G -l sps=1 -j y  \
	      -o $OUT_DIR/logFiles/$jobName.log buildCoincbatch.csh
   done
else
    source buildCoincbatch.csh run detector trigchan SOFT_DIR DATA_DIR OUT_DIR jobName ProcessedData_Dir SimuData_Dir
fi 
#################<---- Boucle Batch qsub
