#!/bin/bash


echo "---------------------------------"
date
echo "---------BEGIN SCRIPT------------"

export JOBDIR="$PWD"
cd $JOBDIR
echo $JOBDIR
ls -lrt
echo "Penser a ENLEVER JOBDIR remplacer par PWD"

echo "rootsys="$ROOTSYS
#export ROOTSYS="/usr/local/root/6.08.02"
export PATH=/pbs/throng/edelweis/bin:$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH



echo "det=$det , run=$runNb"
echo "LineDirectory : $LineDirectory"
echo "copy calibration files"
echo "cp $CALIBFILEDIR/*.NEPALcalib $JOBDIR"
cp $CALIBFILEDIR/*.NEPALcalib $JOBDIR
echo "copy executable"
echo "cp $CC_EVEREST_PROG_DIR/everest.exe $JOBDIR"
cp $CC_EVEREST_PROG_DIR/everest.exe $JOBDIR
echo "copy files to calibrate"
echo "cp -r $LineDirectory*/ $JOBDIR"
cp -r $LineDirectory*/ $JOBDIR

export LineDirectory_withoutpath=` echo $LineDirectory | awk -F"/$det/" '{print $2}'`
echo "LineDirectory_withoutpath : $LineDirectory_withoutpath"

export LineDirectory_JOBDIR=$JOBDIR/$LineDirectory_withoutpath
echo "LineDirectory_JOBDIR=$LineDirectory_JOBDIR"



for i in {1..2} ; do

    for LineDirectory_AND_Iteration_JOBDIR in ` ls -1d $LineDirectory_JOBDIR*` ; do
  
        for LineDirectory_AND_Iteration_AND_File_JOBDIR in ` ls -1d $LineDirectory_AND_Iteration_JOBDIR/Processed*.root` ; do
   

        export inFile=$LineDirectory_AND_Iteration_AND_File_JOBDIR
        export outFile=`echo "$inFile" | sed -e 's/ProcessedData/CalibratedData/g'`
 
        
        if (($i==1))
        then
            echo "################## CALIBRATION ####################################"
            echo "inputfile : $inFile"
            echo "outputfile : $outFile"
    
             echo "$JOBDIR/everest.exe -d $JOBDIR/ -r $bigRunNb -i $inFile -o $outFile -s $runNb"
             $JOBDIR/everest.exe -d $JOBDIR/ -r $bigRunNb -i $inFile -o $outFile -s $runNb
        
        else
            echo "################# COPY BACK TO SPS ##################################"
            echo "inputfile JORBIR : $inFile"
            export outFilenopath=`echo "$outFile" | awk -F"$JOBDIR/" '{print $2}' | sed -e "s/CalibratedData/CalibratedData\/CalibratedData/g"` 
            export outFile_SPS="$CCDATADIR/Simu/$runNb/$det/$outFilenopath"
            echo "outFile_SPS : $outFile_SPS"
            echo "mv $outFile $outFile_SPS"
            mv $outFile $outFile_SPS
        fi

        done

    done

done


echo "end execution script"

exit


