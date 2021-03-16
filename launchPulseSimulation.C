void launchPulseSimulationNEGATIVE_prodHugues()
{
   

	std::vector<string> vecdet;
	vecdet.push_back("NbSi209");

	std::vector<string> vecrun;
   
 //------------------------------------RUN 66V NbSi 10keV  ------------------      
   /*vecrun.push_back("td16a001");
   vecrun.push_back("td17a000");
   vecrun.push_back("td18a000");
   vecrun.push_back("td19a000");*/
   vecrun.push_back("td20a000");
   vecrun.push_back("td22a000");
   

//------------------------------------RUN 66V NbSi pour migdal ------------------   
    //vecrun.push_back("tk19a000");
    //vecrun.push_back("tk20a002");
    //vecrun.push_back("tk22a000");
    //vecrun.push_back("tk24a001");
    //vecrun.push_back("tk26a000");
    //vecrun.push_back("tl15a001");
    //vecrun.push_back("tl15a002");
    //vecrun.push_back("tl18a000");
    //vecrun.push_back("tl22a000");
    //vecrun.push_back("tl24a000");
   /* vecrun.push_back("tc17a000");
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
    vecrun.push_back("uf09a000");*/
    



    string templateNormal ="te30a003"; 
    string templateFast   ="te30a004";
    string templateSlow   ="th03a000";
    string templateSpike  ="te30a004";
    string templateNTD    ="te30a004";



    bool force=false; // si force==false

    int boolFakeData=1;
    string N_Part_Num="All";  // "All" or "1" or "2", etc...
    
    string outputdirectory="/sps/edelweis/rootDataRun317/streams/Efficiency_rought_Migdal/NbSi_chihalf_15V_Th"; // plusieur repertoire de simu Efficiency_rought_Migdal/NbSi_2
    if(boolFakeData==0){
          //outputdirectory=outputdirectory+"/Data";
    }else{outputdirectory=outputdirectory+"/Simu";}



	/////////////////////////////////////////////
	
	
	int NSamplesHeat=1024;  // Number of samples in the TimeWindow
	int trigchannel=1;		// 0=ChalA	1=ChalB
	int tunningRMS=0;		// 
    double CutOff_frequency=1;
	int Filter_Order=2;		// 2 usually does a good job
	double Noise_Rate=5;	// [Hz]	usually 0.5 should be of the order of the event rate roughly
	double borne_inf=0.003;	// 
	double borne_sup=0.6;	// 
    
    double conf_level_inf=0.1; // 0.1
    double conf_level_sup=0.8; // 0.8

	double Trigger_NSIG=0;  // 
	
	/////////////////////// DO NOT CHANGE THE FOLLOWING //////////////////////
    int NumOfTemplates=2;
    string TempType="NormalSlow";
    string FirstType="Normal";
    string SecondType="Fast";
    string ThirdType="Slow";
    string FourthType="NTD";
    string FifthType="Spike";
	string script="/pbs/throng/edelweis/StreamNEPALModane_Migdal_analysis/streamnepalmodane/Processing/runStreamProcess3col_Quentin";
	string cmdtype="Batch";
	string rawdatadirectory="/sps/edelweis/RawData/Run317/streams";
	string softdirectory="/pbs/throng/edelweis/StreamNEPALModane_Migdal_analysis/streamnepalmodane";
	string templatedirectory="/sps/edelweis/rootDataRun317/streams/prodj";//_Midgal_eff_estimation";
	
    ////////////// WHAT IS BELOW ONLY MATTERS FOR FAKE DATA //////////////
    string ProcessedData_Dir="/sps/edelweis/rootDataRun317/streams/prodj";//_Midgal_eff_estimation";
	string Library_Dir="/sps/edelweis/rootDataRun317/streams/prodj";
    string Library_Name="NbSi_acti_66_half";
    double Template_Energy=10.37; // in keV
    double Simu_Emin=10.;
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
   // vecLines.push_back((8.)/1000.);
   // vecLines.push_back((1226.)/1000.);
    
    int NLINES=(int)vecLines.size();
	int NDET=(int)vecdet.size();
	int NRUNS=(int)vecrun.size();		
	cout<<"The following values will be used for the processing \n"<<endl;
	//cout<<"Frequency = "<<Frequency<<" Hz"<<endl;
	//cout<<"Time Window = "<<Time_Window<<" s"<<endl;
	cout<<"Trig Channel = "<<trigchannel<<" "<<endl;
	cout<<"NSamplesHeat = "<<NSamplesHeat<<" "<<endl;
	cout<<"Tunning RMS = "<<tunningRMS<<" "<<endl;
	cout<<"CutOff_frequency = "<<Form("%1.0f",CutOff_frequency)<<" Hz"<<endl;
	cout<<"Filter_Order = "<<Filter_Order<<endl;
	cout<<"Noise_Rate = "<<Form("%2.1f",Noise_Rate)<<" Hz"<<endl;
	cout<<"Born_inf = "<<Form("%2.4f",borne_inf)<<" "<<endl;
	cout<<"Born_sup = "<<Form("%2.4f",borne_sup)<<" "<<endl;
    cout<<"conf_level_inf = "<<Form("%2.2f",conf_level_inf)<<" "<<endl;
	cout<<"conf_level_sup = "<<Form("%2.2f",conf_level_sup)<<" "<<endl;
	cout<<"N_SIG_Trigger = "<<Form("%2.1f",Trigger_NSIG)<<endl;
	cout<<"Template Normal used : "<<templateNormal<<endl;
    cout<<"Template Fast used : "<<templateFast<<endl;
    cout<<"Template Slow used : "<<templateSlow<<endl;
    cout<<"Template Spike used : "<<templateSpike<<endl;
    cout<<"Template NTD used : "<<templateNTD<<endl;
	cout<<"OUTPUT_DIR : "<<outputdirectory<<endl;
    cout<<"\n"<<endl;
    if(boolFakeData){
    cout<<"Processed data dir ="<< ProcessedData_Dir<<endl;
    cout<<" "<<endl;
    cout<<"these are FAKE data"<<endl;
    cout<<"Template : "<<Library_Name<<endl;    
    cout<<"Niterations : "<<Num_Iterations<<endl;
    cout<<"Rate : "<<Simu_Rate<<" Hz"<<endl;
    cout<<"Simu Type : "<<Simu_Type<<endl;
    if(Simu_Type=="Line"){
    cout<<"Line Energies :"<<endl;
     for(int l=0;l<NLINES;l++){cout<<vecLines[l]<<" ";}
     cout<<"\n"<<endl;
      cout<<"total of "<<NLINES<<" lines x "<<NRUNS<<" runs = "<<NLINES*NRUNS<<" jobs"<<endl;
      cout<<"estimated time to launch the jobs = "<<NLINES*NRUNS*1.2/60.<<" hours"<<endl;
    }else{
    cout<<Simu_Emin<<" < Energy (keV) < "<<Simu_Emax<<endl;
    }
    
     
    }else{
    cout<<"these are REAL data"<<endl;
    }
	cout<<"\n Detectors :	";
	for(int d=0;d<NDET;d++){cout<<vecdet[d]<<" ";}
	cout<<"\n Runs :	";
	for(int r=0;r<NRUNS;r++){cout<<vecrun[r]<<" ";}
	cout<<"\n"<<endl;
   
	cout<<"Do you want to continue ? (y or n)"<<endl;
	string answer;

    if(force==false){
	    cin>>answer;
	    if(answer!="y"){cout<<"Processing Canceled"<<endl; exit(1);}
    }
	


if(Simu_Type=="Line" && boolFakeData==1){Nloop=NLINES;}
	
for(int r=0;r<NRUNS;r++){
	for(int d=0;d<NDET;d++){
        
        for(int l=0;l<Nloop;l++){
            if(Simu_Type=="Line"){
                Simu_Emin=vecLines[l];
                Simu_Emax=vecLines[l];
            }        
        
		    string run=vecrun[r];
		    string det=vecdet[d];
		    string cmd1=script+" "+cmdtype+" "+rawdatadirectory+" "+run+" "+det+" ";
		    string cmd2=Form("%i %i %i ",trigchannel,NSamplesHeat,tunningRMS);
            string cmdX=Form("%i %s %s %s %s %s ",NumOfTemplates,TempType.c_str(),templateNormal.c_str(),FirstType.c_str(),templateSlow.c_str(),ThirdType.c_str());
            string cmdY=Form("%1.0f %i %2.1f %2.4f %2.4f %2.1f %2.2f %2.2f ",CutOff_frequency,Filter_Order,Noise_Rate,borne_inf,borne_sup,conf_level_inf,conf_level_sup,Trigger_NSIG);//
		    string cmd3=softdirectory+" "+templatedirectory+" "+outputdirectory+" "+N_Part_Num;          
            string cmd4=Form(" %i %s",boolFakeData,ProcessedData_Dir.c_str());
            string cmd5=Form(" %s %s %.4f",Library_Dir.c_str(),Library_Name.c_str(),Template_Energy);
            string cmd6=Form(" %.4f %.4f %.4f %s %i %i",Simu_Emin,Simu_Emax,Simu_Rate,Simu_Type.c_str(),Num_Iterations,Partition_parity);
          
		    string command=cmd1+cmd2+cmdX+cmdY+cmd3+cmd4;
            if(boolFakeData==1){
                command=command+cmd5+cmd6;
            }
            
   		    cout<<command<<endl;
		    sleep(1); // wait for 2 seconds
		    
            // did you check runs, folder, lines and iterations ?
            system(command.c_str());		

        }
	}	
	
}


//////////////////////////NEW 5 COL/////////////
/*
arguments n0 :  Batch
arguments n1 :  DATADir      (/sps/edelweis/RawData/Run317/streams)
arguments n3 :  RUN          (td16a001)
arguments n4 :  Detector     (RED30)
arguments n5 :  TrigChan     (0)
arguments n6 :  NSampHeat    (1024)
arguments n7 :  TuningRMS    (0)
arguments n8 : NumOfTemplates (2)
arguments n9 : TempType (NormalFast)
arguments n10 :  TemplateNormal (te30a003)
arguments n11 :  FirstType (Normal)
arguments n12 :   TemplateFast (te30a004)
arguments n13 :   SecondType (Fast)
arguments n14 :   TemplateSlow (te30a004)
arguments n15 :   ThirdType (Slow)
arguments n16 :   TemplateNTD (te30a004)
arguments n17 :   FourthType (NTD)
arguments n18 :   TemplateSpike (te30a004)
arguments n19 :   FifthType (Spike)
arguments n20:  CutOff_Freq  (1)  en regle generale 1./(0.5 TimeWindow)
arguments n21:  filter_order (2) :
arguments n22 : Noise_rate   (5)
arguments n23 : b_inf        (0.003)
arguments n24 : b_sup        (0.6)
arguments n25 : Trig_level   (0)
arguments n26 : SOFT_Dir     (/pbs/throng/edelweis/StreamNEPALModane_beta/streamnepalmodane)
arguments n27 : Template_Dir (/sps/edelweis/rootDataRun317/streams/prodg)
arguments n28 : OUT_Dir      (/sps/edelweis/rootDataRun317/streams/test)
arguments n29 : N_part_Num   (All)
arguments n30 : FakeData     (1)
arguments n31 : Processed_Data_Dir   (/sps/edelweis/rootDataRun317/streams/prodg)
arguments n32 : Library_Dir   /sps/edelweis/rootDataRun317/streams/prodg )
arguments n33 : Library_Name (RED30-fond-78v-acti)
arguments n34 : Template_Energy (10.37)
arguments n35 : E_min       (0.)
arguments n36 : E_max       (1.)
arguments n37 : Simu_Rate   (0.1)
arguments n38 : Simu_Type   (Flat)
arguments n39 : Num_Iteration   (20) 
*/

}
