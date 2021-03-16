################## Efficiency_edelweis ##################
#          Soft framework to compute effyciency         #
#          based on NEPAL Pulse simulation output       #
#                     13/03/2021                        #
#          author : Quentin Arnaud and Hugues Lattaud   #
#########################################################

This framework aims to build effyciency curve in the context of
the Edelweiss collaboration data analysis. These code will allow you 
to launch simulation and analyse the output.

The pulse simulation is a software developped by Julien Billard, which 
inject rescaled physical event into data stream. It allow to scan detectors
response in a large energy range. 

The first step consist in building the reference event bank you will later on
inject.

***********************Step 1 : Reference event bank Building*****************

This step is done by  NEPAL_Template_Select, it takes as inputs, the path to the Nepal processed 
data, the path to the triggered Trace,   the path to the templates, the detector name, the channel you trigger on, a list file of the run to scan for event selection, a boolean to set at false (obsolete option).
$make clean
$make
$./NEPAL_Template_Select.exe /sps/edelweis/rootDataRun317/streams/prodk /sps/edelweis/rootDataRun317/streams/prodk /sps/edelweis/rootDataRun317/streams/prodk NbSi209 1 NbSi_acti_66 0

You can adjust the selection criterion by hand in the code (ugly old school way).
