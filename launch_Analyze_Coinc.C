void launch_Analyze_Coinc()
{
   
    std::vector<string> vecdet;
	vecdet.push_back("NbSi209");

	std::vector<string> vecrun;

    vecrun.push_back("tk19a000");
    vecrun.push_back("tk20a002");
    vecrun.push_back("tk22a000");
    vecrun.push_back("tk24a001");
    vecrun.push_back("tk26a000");
    vecrun.push_back("tc17a000");
    vecrun.push_back("tc19a000");
    vecrun.push_back("tc19a000");
    vecrun.push_back("td02a000");
    vecrun.push_back("td04a000");
    vecrun.push_back("td03a000");
    vecrun.push_back("tj27a001");
    vecrun.push_back("tj28a000");
    vecrun.push_back("tj29a000");
    vecrun.push_back("ud21a000");
    vecrun.push_back("ud22a000");
    vecrun.push_back("ud23a000");
    vecrun.push_back("ud24a000");
    vecrun.push_back("ud27a001");
    vecrun.push_back("ud29a000");
    vecrun.push_back("ue25a001");
    vecrun.push_back("ue27a000");
    vecrun.push_back("ue28a000");
    vecrun.push_back("ue29a000");
    vecrun.push_back("ue30a000");
    vecrun.push_back("uf01a000");
    vecrun.push_back("uf02a000");
    vecrun.push_back("uf05a001");
    vecrun.push_back("uf07a000");
    vecrun.push_back("uf08a000");
    vecrun.push_back("uf09a000");
    /*vecrun.push_back("td16a001");
    vecrun.push_back("td17a000");
    vecrun.push_back("td18a000");
    vecrun.push_back("td19a000");
    vecrun.push_back("td20a000");
    vecrun.push_back("td22a000");*/
    
    string coinc_directory="/sps/edelweis/rootDataRun317/streams/Efficiency_rought_Migdal/NbSi_chihalf_15V_Th"; // plusieur repertoire de simu Efficiency_rought_Migdal/NbSi_2
    coinc_directory=coinc_directory+"/SimuCoinc";


    string Tension = "66.";
    string Ion_cut = "0.50" ;//in keV if negative cut -> no cut 
	
	/////////////////////// DO NOT CHANGE THE FOLLOWING //////////////////////

	string script="/sps/edelweis/rootDataRun317/streams/Selection10keV_Migdal/runAnalyze_Coinc.sh";
	string cmdtype="Batch";//"Batch";//
	string rawdatadirectory="/sps/edelweis/RawData/Run317/streams";
	string softdirectory="/sps/edelweis/rootDataRun317/streams/Selection10keV_Migdal";
	
    ////////////// WHAT IS BELOW ONLY MATTERS FOR FAKE DATA //////////////
    string ProcessedData_Dir="/sps/edelweis/rootDataRun317/streams/prodj";
    string SimuData_Dir="/sps/edelweis/rootDataRun317/streams/Efficiency_rought_Migdal/NbSi_chihalf_15V_Th";
    double Template_Energy=10.37; // in keV
    double Simu_Emin=1.;
    double Simu_Emax=1500./1000.;
    double Simu_Rate=0.02 ; //Hz
    string Simu_Type="Line"; // Line or Flat
    int Nloop=1; // ne pas toucher Nloop

    int Num_Iterations=5;
    int Partition_parity = 0 ; // 0 stand for All , 1 for odd , 2 for pair.
    ///////////////////////////////////////
    // THIS MATTERS ONLY IF Type==Line   //
    std::vector<double> vecLines;        // 
    for(int iterator = 1 ; iterator < 1500. ; iterator+=7) {       
       vecLines.push_back((iterator)/1000.);
    }
   
    
    int NLINES=(int)vecLines.size();
	int NDET=(int)vecdet.size();
	int NRUNS=(int)vecrun.size();
	cout<<"JOB TYPE "<<cmdtype<<" \n"<<endl;		
	cout<<"The following values will be used for the COINC BUILDING \n"<<endl;
	//cout<<"Frequency = "<<Frequency<<" Hz"<<endl;
	//cout<<"Time Window = "<<Time_Window<<" s"<<endl;
	cout<<"COINC_DIR : "<<coinc_directory<<endl;
    cout<<"\n"<<endl;
    cout<<"Processed data dir ="<< ProcessedData_Dir<<endl;
    cout<<" "<<endl;
    cout<<"Bias ="<< Tension<<endl;
    cout<<" "<<endl;
    cout<<"Ionization cut ="<< Ion_cut <<endl;
    cout<<" "<<endl;
    cout<<"these are FAKE data"<<endl;
    cout<<"Niterations : "<<Num_Iterations<<endl;
    cout<<"Rate : "<<Simu_Rate<<" Hz"<<endl;
    cout<<"Simu Type : "<<Simu_Type<<endl;
    cout<<"Line Energies :"<<endl;
     for(int l=0;l<NLINES;l++){cout<<vecLines[l]<<" ";}
     cout<<"\n"<<endl;
      cout<<"total of "<<NLINES<<" lines x "<<NRUNS<<" runs = "<<NLINES*NRUNS<<" jobs"<<endl;
      cout<<"estimated time to launch the jobs = "<<NLINES*NRUNS*1.2/60.<<" hours"<<endl;
    
	cout<<"\n Detectors :	";
	for(int d=0;d<NDET;d++){cout<<vecdet[d]<<" ";}
	cout<<"\n Runs :	";
	for(int r=0;r<NRUNS;r++){cout<<vecrun[r]<<" ";}
	cout<<"\n"<<endl;
   
	cout<<"Do you want to continue ? (y or n)"<<endl;
	string answer;
	cin>>answer;
	if(answer!="y"){cout<<"Processing Canceled"<<endl; exit(1);}
    
	
int triggerchannel = 1 ;
Nloop=NLINES;
	
for(int r=0;r<NRUNS;r++){
	for(int d=0;d<NDET;d++){
              
		    string run=vecrun[r];
		    string det=vecdet[d];
		    string cmd1=script+" "+cmdtype+" "+rawdatadirectory+" "+run+" "+det+" "+std::to_string(triggerchannel)+" ";	            
		    string cmd2=softdirectory+" "+coinc_directory;          
            string cmd3=Form(" %s ",ProcessedData_Dir.c_str());
            string cmd4=Form(" %s %s %s ",SimuData_Dir.c_str(), Tension.c_str(), Ion_cut.c_str());          
		    string command=cmd1+cmd2+cmd3+cmd4;                      
   		    cout<<command<<endl;
		    sleep(1); // wait for 2 seconds		    
            // did you check runs, folder, lines and iterations ?
            system(command.c_str());		

        }
	}	
	



//////////////////////////NEW 5 COL/////////////
/*
arguments n0  :  Batch
arguments n1  :  DATADir      (/sps/edelweis/RawData/Run317/streams)
arguments n3  :  RUN          (td16a001)
arguments n4  :  Detector     (RED30)
arguments n5  : SOFT_Dir     (/pbs/throng/edelweis/StreamNEPALModane_beta/streamnepalmodane)
arguments n6  : OUT_Dir      (/sps/edelweis/rootDataRun317/streams/test)
arguments n7  : Processed_Data_Dir   (/sps/edelweis/rootDataRun317/streams/prodg)
arguments n8 : Simu_Dir
*/

}
