#include <algorithm>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <stdlib.h> 
#include <istream>
#include <time.h>
#include <glob.h>
#include <vector>
#include "TROOT.h"
#include "TSystem.h"
#include "THStack.h"
#include "TPaveStats.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixTSym.h"
#include "TDecompLU.h"
#include "TApplication.h"
#include "TLatex.h"
#include "TRint.h"
#include "TLegend.h"
#include "TEventList.h"
#include "TMatrixD.h"
#include "TGraph2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TArrow.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH2F.h"
#include "TF3.h"
#include "Math/Functor.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TMath.h"
#include "TChain.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Config.h"
#include <TPad.h>
#include "TEfficiency.h"
#include "TVirtualFFT.h"
#include "TGraphErrors.h"
#include "TComplex.h"
#include "TVirtualFFT.h"
#include "Math/Functor.h"
#include "TColor.h"
#include <X11/Xlib.h>
using namespace std;

const int N_Max_Channel = 6;

template <typename VecType>
void freeFromMemory(VecType& myVec) {
     myVec.clear();
     VecType vert;
     myVec.swap(vert);
}



double CalculateDeltaT(double const& value,std::vector<double> const& vec);
void LoadBinary(string const& binarypath,string const& run,int const& RunPart,string const& Channel,std::vector <string> const& Detector_List,double *amplitudes,int const& detID,int const& NB_Samples_to_read_per_channel,int const& stampstart);


int main(int argc, char *argv[])
{
	
	// Delta_chi2_Normalized with 1024,  should implement N_Samples
	//string study_index="FOND_LINES_8";
	string study_index="Efficiency_rought_Migdal/NbSi_chihalf_15V_Th";
    cout<<study_index<<endl;
  //  exit(1);
	string DATA_Path="/sps/edelweis/rootDataRun317/streams/prodj"; // repertoire données


	string SIMU_Path=Form("/sps/edelweis/rootDataRun317/streams/%s/Simu",study_index.c_str()); // repertoire SIMU
	string SIMUCOINC_Path=Form("/sps/edelweis/rootDataRun317/streams/%s/SimuCoinc",study_index.c_str());
	
	string det="NbSi209"; // only 1 detector
	int Trigger_channel = 0 ;
	string template_2 = "Fast";
	if(det == "NbSi209" || det == "FID848"){
	   Trigger_channel= 1;
	   template_2 = "Slow";
  }
	bool usechihalf = false;
	double Ndof_chi2 = 1024;
	if(usechihalf) Ndof_chi2 = 512 ;
	
	// not used anymore ///
	string binarypath=Form("/data/PulseSimulationFiles/study1/Stream");
	
	
	// AJOUT POUR BINAIRES ///
	std::vector<string> vecDet;
	vecDet.push_back(det);
	int detID=0;
	//////////////////////////
	
	string run_to_process = argv[1];

	std::vector<string> vecruns={run_to_process}; // runs to processed  //"tl21a000","tl20a000","tg13a000","tj26a000","ub13a000","th23a001","tg13a001","ub12a000"

	const int Nbruns=(int)vecruns.size();
	std::vector<string> vec_Lines;

	// ici
	double SIMU_Emin=0; // should only matter for Flat
	double SIMU_Emax=1; // should only matter for Flat
	string SIMU_Type="Line"; // or line
	int SIMU_NTOT_iterations=10000; // will take the min of chosen and available NTOT_iterations;
	int DATA_NFiles=10000; // will take the min of chosen and available
	
	int fracmult_N_samples=1;
	
	bool useEverest=true;
	
    bool LoadTracesfromBinaries=false;  // won't work if channel=1 , need to work on that

	
		
	//////////////////////////////////////////////////////////

	
		
	if(SIMU_Type=="Line"){
	
	///////////////// LOOK FOR LINES ////////
	for(int r=0;r<Nbruns;r++){
		string commandline=Form("ls -d %s/%s/%s/Line_*_* | awk -F\"/Line_\" '{print $2}' | awk -F\"_\" '{print $1}' | sort -n | uniq > templist_%s_Lines.txt",SIMU_Path.c_str(),vecruns[r].c_str(),det.c_str(),vecruns[r].c_str());
		system(commandline.c_str());
	}
	
	bool difference_N_Lines=false; // always set false
	for(int r=0;r<Nbruns;r++){
		for(int rr=0;rr<Nbruns;rr++){
			if(rr>r){
				string commandline=Form("diff templist_%s_Lines.txt templist_%s_Lines.txt",vecruns[r].c_str(),vecruns[rr].c_str());
				cout<<commandline<<endl;
				system(commandline.c_str());
				string Ndiffcommand=Form("diff templist_%s_Lines.txt templist_%s_Lines.txt | wc -l",vecruns[r].c_str(),vecruns[rr].c_str());
				int Ndiff = gSystem->GetFromPipe(Ndiffcommand.c_str()).Atoi(); // th
				if(Ndiff>0){difference_N_Lines=true;}
				//cout<<"Ndiff="<<Ndiff<<endl;
			}
		}
		
	}
	
	string answer="n";
	if(difference_N_Lines==true){
		cout<<"number of lines different \n"<<endl;
		cout<<"you can't recreate all SimuCoinc files at once it won't work\n"<<endl;
		cout<<"are you sure you want to continue ? enter 'y' if so"<<endl;
		cin>>answer;
		if(answer!="y"){cout<<"exiting program then"<<endl; exit(1);}
	}
	
	
	
	string commandline=Form("ls -d %s/%s/%s/Line_*_* | awk -F\"/Line_\" '{print $2}' | awk -F\"_\" '{print $1}' | sort -n | uniq > templist_Lines.txt",SIMU_Path.c_str(),vecruns[0].c_str(),det.c_str());
	cout<<commandline<<endl;
	system(commandline.c_str());
	ifstream fline;
	fline.open("templist_Lines.txt");
	
	string buffer;
	while(1){
		fline>>buffer;
		if(!fline.good()){break;}
		vec_Lines.push_back(buffer);
	}
	
	cout<<vec_Lines.size()<< "Lines found"<<endl;
	for(int i=0;i<(int)vec_Lines.size();i++){
		cout<<i<<"	"<<vec_Lines[i]<<endl;
	}
	

	//SIMU_Emin=2.9;
	/////////////////////////////////////////
	}	
	
	

	int NFILES_BUG=0;
	std::vector<string> filesbugged;

	std::vector<string> filesempty;
	




	for(int r=0;r<Nbruns;r++){
	
		string run=vecruns[r];
		cout<<"run "<<run<<endl;
		
		string allrootfiles=Form("%s/%s/%s/ProcessedData_%s_*_%s_ChanTrig%i.root",DATA_Path.c_str(),run.c_str(),det.c_str(),run.c_str(),det.c_str(),Trigger_channel);
		
		string command=Form("ls %s | wc -l",allrootfiles.c_str());
		int nfiles = gSystem->GetFromPipe(command.c_str()).Atoi(); // the bash command returns a TString that is converted into int with Atoi();
	   	//cout<<"nfiles found :"<<nfiles<<endl;
	   	if(nfiles==0){cout<<"file does not exist"<<endl; exit(1);}
		
		nfiles=min(DATA_NFiles,nfiles);
		//cout<<"nfiles used :"<<nfiles<<endl;
		//nfiles = 13;
		for(int filenumber=0;filenumber<nfiles;filenumber++){
			
	   		//cout<<"File Number "<<filenumber<<endl;
	   	string format;
		   if(filenumber<10){
		      format="S0";
		   }else{
		      format="S";
		   }
	   	/////////////////////////// LOAD REAL DATA /////////////////////////////
		   	
		   string DATA_File=Form("ProcessedData_%s_%s%i_%s_ChanTrig%i.root",run.c_str(),format.c_str(),filenumber,det.c_str(),Trigger_channel);
		   string DATA_PathAndFile=Form("%s/%s/%s/%s",DATA_Path.c_str(),run.c_str(),det.c_str(),DATA_File.c_str());
		   	
			TFile *fdata=new TFile(DATA_PathAndFile.c_str(),"READ");
			TTree *tree_data=(TTree*)fdata->Get("EventTree_trig_Normal_filt");
			int Nentries_Data=tree_data->GetEntries();
			int t0_timestamp_Data;
		   	tree_data->SetBranchAddress("MicroStp",&t0_timestamp_Data);
			
			std::vector<double> vectimes_Data;
			
			for(int entry=0;entry<Nentries_Data;entry++){
		   		tree_data->GetEntry(entry);
		   		vectimes_Data.push_back((double)t0_timestamp_Data);
		   	}
			
			TTree *tree_data_info=(TTree*)fdata->Get("RunTree_Normal");
			int N_samples_Heat;
			double TimeWindow_Heat;
			tree_data_info->SetBranchAddress("N_samples_Heat",&N_samples_Heat);
			tree_data_info->SetBranchAddress("TimeWindow_Heat",&TimeWindow_Heat);
			tree_data_info->GetEntry(0);
			int d3=N_samples_Heat/TimeWindow_Heat;
			//int d2=1e5/d3;
			N_samples_Heat=N_samples_Heat*fracmult_N_samples;
			
			
			double Trace_Data_at_t0[N_samples_Heat];
			
			/////////////////////////// INIT SIMU //////////////////////////////////// 
			// determine the number of iterations
			int NTYPES=1;
			
			if(SIMU_Type=="Line"){
			NTYPES=(int)vec_Lines.size();
			}
			
			for(int t=0;t<NTYPES;t++){
				string pattern;
				if(SIMU_Type=="Flat"){
					pattern=Form("%s/%s/%s/%s_%.4f_%.4f_*",SIMU_Path.c_str(),run.c_str(),det.c_str(),SIMU_Type.c_str(),SIMU_Emin,SIMU_Emax);
				}else if(SIMU_Type=="Line"){
					pattern=Form("%s/%s/%s/%s_%s_*",SIMU_Path.c_str(),run.c_str(),det.c_str(),SIMU_Type.c_str(),vec_Lines[t].c_str());
				}else{cout<<"this simutype doesn't exist"<<endl; exit(1);}
				
				
				string command=Form("ls -d %s | wc -l",pattern.c_str());
		     	int niterations_temp = gSystem->GetFromPipe(command.c_str()).Atoi(); // the bash command returns a TString that is converted into int with Atoi();   		
				SIMU_NTOT_iterations=min(niterations_temp,SIMU_NTOT_iterations);
				cout<<SIMU_NTOT_iterations<<" iterations for "<<SIMU_Type.c_str()<<" "<<vec_Lines[t].c_str()<<endl;
				
				
				for (int SIMU_iteration=1;SIMU_iteration<=SIMU_NTOT_iterations;SIMU_iteration++){
					
					cout<<Form("run %d/%d , file %d/%d , Type %d/%d , iteration %d/%d",r,Nbruns,filenumber,nfiles,t,NTYPES,SIMU_iteration,SIMU_NTOT_iterations)<<endl;
                   
				   	string SIMU_condition;
				   	if(SIMU_Type=="Flat"){
				   	   SIMU_condition=Form("%s_%.4f_%.4f_%i",SIMU_Type.c_str(),SIMU_Emin,SIMU_Emax,SIMU_iteration);
				   	}else if(SIMU_Type=="Line"){
				   	   SIMU_condition=Form("%s_%s_%i",SIMU_Type.c_str(),vec_Lines[t].c_str(),SIMU_iteration);
				   	}else{
				   	   cout<<"this simutype doesn't exist"<<endl; exit(1);
				   	 }
					
					   string SIMU_RELATIVE_Path=Form("%s/%s/%s/%s",SIMU_Path.c_str(),run.c_str(),det.c_str(),SIMU_condition.c_str());
					   string SIMU_File=Form("ProcessedData_%s_%s%i_%s_%s_ChanTrig%i.root",run.c_str(),format.c_str(),filenumber,det.c_str(),SIMU_condition.c_str(),Trigger_channel);
					   string SIMU_PathAndFile=SIMU_RELATIVE_Path+"/"+SIMU_File;
					//cout<<"simu : "<<SIMU_PathAndFile.c_str()<<endl;
					   string SIMU_EVEREST_File=Form("CalibratedData_%s_%s%i_%s_%s_ChanTrig%i.root",run.c_str(),format.c_str(),filenumber,det.c_str(),SIMU_condition.c_str(),Trigger_channel);
					   string SIMU_EVEREST_PathAndFile=SIMU_RELATIVE_Path+"/CalibratedData/"+SIMU_EVEREST_File;
					//cout<<"simu : "<<SIMU_EVEREST_PathAndFile.c_str()<<endl;
				   	string SIMUCOINC_File=Form("SIMUCOINC_%s",SIMU_File.c_str());
				   	string SIMUCOINC_RELATIVE_Path=Form("%s/%s/%s/%s",SIMUCOINC_Path.c_str(),run.c_str(),det.c_str(),SIMU_condition.c_str());
				   	string SIMUCOINC_PathAndFile=SIMUCOINC_RELATIVE_Path+"/"+SIMUCOINC_File;
				   	//cout<<"simu coinc :"<<SIMUCOINC_PathAndFile<<endl;
				   	////////////////// COPYING ROOTFILES IN SIMULATIONSCOINC DIRECTORY /////////////////////////////////////////
				   	string commandcopyfiles=Form("mkdir -p %s ; cp %s %s",SIMUCOINC_RELATIVE_Path.c_str(),SIMU_PathAndFile.c_str(),SIMUCOINC_PathAndFile.c_str());	 
				   	//cout<<"SIMU COINC PATH AND FILE : "<<SIMUCOINC_PathAndFile<<endl;  	
				   	system(commandcopyfiles.c_str());
			 	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
				   	   
				   	   
				   	TFile *fsimu=new TFile(SIMU_PathAndFile.c_str(),"READ");
				   	   	
					
				   	////////////////////////////////// LOAD SIMU ///////////////////////////////////////////
				   	
				   	
				   	TTree *tree_Trig_SimuOnly_t0=(TTree*)fsimu->Get("tree_Trig_SimuOnly_t0");
				   	int t0_timestamp_SimuOnly;
				   	tree_Trig_SimuOnly_t0->SetBranchAddress("t0_timestamp",&t0_timestamp_SimuOnly);
				   	
				   	TTree *tree_Input_SimuOnly_Energy=(TTree*)fsimu->Get("tree_Input_SimuOnly_Energy");
				   	double Energy_In;
				   	tree_Input_SimuOnly_Energy->SetBranchAddress("Energy_In",&Energy_In);

                  TTree *tree_Input_SimuOnly_EnergyOFInput=(TTree*)fsimu->Get("tree_Input_SimuOnly_EnergyOFInput");
	                double Energy_OFInput[N_Max_Channel];
				   	tree_Input_SimuOnly_EnergyOFInput->SetBranchAddress("Energy_OFInput",&Energy_OFInput);
				   	
				   	////////////////////////////////// LOAD PROCESSED ROOT FILE OF SIMU EVENTS ////////////////////////
				   	TTree *EventTree_trig_Normal_filt=(TTree*)fsimu->Get("EventTree_trig_Normal_filt");
				   	double Energy_OF_Normal[N_Max_Channel];
				   	double chi2_OF_Normal[N_Max_Channel];
				   	double chi2_OF_half[N_Max_Channel];
				   	double Time_modulo_reset;
				   	int Pass_Cut_Maintenance;
				   	double Delta_chi2_Normalized;
				   	
				   	EventTree_trig_Normal_filt->SetBranchAddress("Energy_OF",&Energy_OF_Normal);
				   	EventTree_trig_Normal_filt->SetBranchAddress("chi2_OF",&chi2_OF_Normal);
				   	if(det=="NbSi209")EventTree_trig_Normal_filt->SetBranchAddress("chi2_OF_half",&chi2_OF_half);
				   	EventTree_trig_Normal_filt->SetBranchAddress("Time_modulo_reset",&Time_modulo_reset);
    					EventTree_trig_Normal_filt->SetBranchAddress("Pass_Cut_Maintenance",&Pass_Cut_Maintenance);
				   	
				   	TTree *EventTree_trig_Fast_filt=(TTree*)fsimu->Get(Form("EventTree_trig_%s_filt",template_2.c_str()));
				   	double chi2_OF_Fast[N_Max_Channel];
				   	EventTree_trig_Fast_filt->SetBranchAddress("chi2_OF",&chi2_OF_Fast);
				   	// pour le NbSi penser à utiliser le chi2/2;				   	
				  	
				   	///////////////////////////////// LOAD CALIBRATED ROOT FILE OF SIMU EVENTS ////////////////////////
				   	// EVEREST				   	
				   	TFile *fsimuEVEREST;
				   	TTree *Energies_Trig_Filt;
				   	double EhA;
				   	double EhB;
				    double Eh;
				    double EiA;
				    double EiB;
				    double EiC;
				    double EiD;
				    double Ei;// will be (Eia + EiB)/2 pour NbSi 209
                  cout<<"before"<<endl;
                  cout<<"n events in EvenTree_trig_Normal ="<<EventTree_trig_Normal_filt->GetEntries()<<endl;
                  bool empty=false; 
                  if(EventTree_trig_Normal_filt->GetEntries()==0){
                     filesempty.push_back(SIMU_PathAndFile);
                     empty=true;
                  }
                  if((int)filesempty.size()>=1){
                      cout<<"list of files empty:"<<endl;
                      for(int xx=0;xx<(int)filesempty.size();xx++){
                         cout<<filesempty[xx]<<endl;
                   }
                  }
                    
				   	if(useEverest==true && empty==false){
                     cout<<"non empty"<<endl;
				   		fsimuEVEREST=new TFile(SIMU_EVEREST_PathAndFile.c_str(),"READ");
				   		Energies_Trig_Filt=(TTree*)fsimuEVEREST->Get("Energies_Trig_Filt");
				   		Energies_Trig_Filt->SetBranchAddress("EhA",&EhA);
				   		Energies_Trig_Filt->SetBranchAddress("EhB",&EhB);
				   		Energies_Trig_Filt->SetBranchAddress("Eh",&Eh);
				   		Energies_Trig_Filt->SetBranchAddress("EiA",&EiA);
				   		Energies_Trig_Filt->SetBranchAddress("EiB",&EiB);
				   		Energies_Trig_Filt->SetBranchAddress("EiC",&EiC);
				   		Energies_Trig_Filt->SetBranchAddress("EiD",&EiD);
				   		Energies_Trig_Filt->SetBranchAddress("Ei",&Ei);
                         if(Energies_Trig_Filt->GetEntries()!=EventTree_trig_Normal_filt->GetEntries()){
                           cout<<"problem, Number of events different in calibrated files"<<endl; exit(1);
                         }
				   	}                    
				   	////////////////// OPENING COPY OF SIMU CALLED COINC SIMU AND ADD NEW TREE 
				   	TFile *fcoinc=new TFile(SIMUCOINC_PathAndFile.c_str(),"UPDATE");
				   	fcoinc->cd();
				   	TTree *tcoinc = new TTree("tsimu","All info");
                  if(empty==true){
                    	fcoinc->cd();
				     	   tcoinc->Write("",TObject::kOverwrite);
				    	   fcoinc->Close();
                     fsimu->Close();
                     cout<<"empty simu tree"<<endl;
                     cout<<"Writing empty tsimucoinc tree"<<endl;
                     continue;
                  }			    	  
				   	double deltaT_simu_datasimu;
				   	double deltaT_simu_data;
				   	
				   	tcoinc->Branch("deltaT_simu_datasimu",&deltaT_simu_datasimu,"deltaT_simu_datasimu/D");
				   	tcoinc->Branch("deltaT_simu_data",&deltaT_simu_data,"deltaT_simu_data/D");				   	
				   	if(useEverest==true){
				   		tcoinc->Branch("EhA",&EhA,"EhA/D");
				   		tcoinc->Branch("EhB",&EhB,"EhB/D");
				   		tcoinc->Branch("Eh",&Eh,"Eh/D");
				   		tcoinc->Branch("EiA",&EiA);
				   		tcoinc->Branch("EiB",&EiB);
				   		tcoinc->Branch("EiC",&EiC);
				   		tcoinc->Branch("EiD",&EiD);
				   		tcoinc->Branch("Ei",&Ei);
				   	}
				   	tcoinc->Branch("Energy_OF_Normal",&Energy_OF_Normal,Form("Energy_OF_Normal[%i]/D",N_Max_Channel));
                    tcoinc->Branch("Energy_OFInput",&Energy_OFInput,Form("Energy_OFInput[%i]/D",N_Max_Channel));
				   	tcoinc->Branch("Time_modulo_reset",&Time_modulo_reset,"Time_modulo_reset/D");
				   	tcoinc->Branch("Pass_Cut_Maintenance",&Pass_Cut_Maintenance,"Pass_Cut_Maintenance/I");
				   	tcoinc->Branch("chi2_OF_Normal",&chi2_OF_Normal,Form("chi2_OF_Normal[%i]/D",N_Max_Channel));
				   	if(det=="NbSi209")tcoinc->Branch("chi2_OF_half",&chi2_OF_half,Form("chi2_OF_half[%i]/D",N_Max_Channel));
				   	tcoinc->Branch(Form("chi2_OF_%s",template_2.c_str()),&chi2_OF_Fast,Form("chi2_OF_%s[%i]/D",template_2.c_str(),N_Max_Channel));
				   	tcoinc->Branch("Delta_chi2_Normalized",&Delta_chi2_Normalized,"Delta_chi2_Normalized/D");
				   	tcoinc->Branch("Energy_In",&Energy_In,"Energy_In/D");
		  			   tcoinc->Branch("t0_timestamp",&t0_timestamp_SimuOnly,"t0_timestamp/I"); // IT IS THE EQUIVALENT TO MicroStp
		  			
		  			   if(LoadTracesfromBinaries==true){
		  			      tcoinc->Branch("Trace_Data_at_t0",Trace_Data_at_t0,Form("Trace_Data_at_t0[%i]/D",N_samples_Heat));
		  			   }
				   	//////////////////////////////// LOAD REAL DATA AND SIMU /////////////////////////////////////////
				
		  	
				   	int t0_timestamp_RealAndSimu;
				   	TTree *tree_Trig_RealAndSimu_t0=(TTree*)fsimu->Get("tree_Trig_RealAndSimu_t0");
				   	tree_Trig_RealAndSimu_t0->SetBranchAddress("t0_timestamp",&t0_timestamp_RealAndSimu);
				   	
				   	int Nentries_RealAndSimu=tree_Trig_RealAndSimu_t0->GetEntries();
				   	std::vector<double> vectimes_Trig_RealAndSimu_t0;
				
				   	for(int entry=0;entry<Nentries_RealAndSimu;entry++){
				   		tree_Trig_RealAndSimu_t0->GetEntry(entry);
				   		vectimes_Trig_RealAndSimu_t0.push_back((double)t0_timestamp_RealAndSimu);
				   	}				   
					// FILL NEW TREE				
				   	int Nentries_SimuOnly=tree_Trig_SimuOnly_t0->GetEntries();
				   	std::vector<double> vectimes_Trig_SimuOnly_t0;	
				   	 				   
				   	for(int entry=0;entry<Nentries_SimuOnly;entry++){				  
				   	//cout<<"entry ="<<entry<<" over "<<Nentries_SimuOnly<<" entries in simu only"<<endl;				   		
				   	tree_Trig_SimuOnly_t0->GetEntry(entry);			
				   	double t0=(double)t0_timestamp_SimuOnly;	
				   	/// AJOUT TRACE ///				   		
				   	if(LoadTracesfromBinaries==true){
				   		   LoadBinary(binarypath,run,filenumber,"chalA",vecDet,Trace_Data_at_t0,detID,N_samples_Heat,t0_timestamp_SimuOnly); cout<<"end load binary"<<endl;
				      }				   						   				   		
				   	// Calculate Delta T for each event				   		
				   	// impossible de pas trigger en principe sur le stream data + simu donc il y a un bug
				   	if(vectimes_Trig_RealAndSimu_t0.size()==0){
				   			deltaT_simu_datasimu=99999;	
				   			if(entry==0){NFILES_BUG++;filesbugged.push_back(SIMU_PathAndFile);}  
				   			
				   	}else{
				  
				   			deltaT_simu_datasimu=CalculateDeltaT(t0,vectimes_Trig_RealAndSimu_t0);
				   	}
				   	deltaT_simu_data=CalculateDeltaT(t0,vectimes_Data);				 		
				 		// Before Fill() GetEntry() for all trees in order to load the variables I want to save 	
				 		EventTree_trig_Normal_filt->GetEntry(entry);
				 		EventTree_trig_Fast_filt->GetEntry(entry);
				 		Delta_chi2_Normalized=(chi2_OF_Normal[0]-chi2_OF_Fast[0])/Ndof_chi2; 
				 		tree_Input_SimuOnly_Energy->GetEntry(entry);
				 		if(useEverest==true){
				   			Energies_Trig_Filt->GetEntry(entry);
				      }				 		
                  tree_Input_SimuOnly_EnergyOFInput->GetEntry(entry);
                  tcoinc->Fill();
				   	} // end loop entries simu
				  
				   	cout<<"out of loop over "<<Nentries_SimuOnly<<" simu entries"<<endl;
				   	fcoinc->cd();
				   	tcoinc->Write("",TObject::kOverwrite);
				   	fcoinc->Close();
				   	
				   	freeFromMemory(vectimes_Trig_RealAndSimu_t0);
				   	fsimu->Close();
				   	fsimuEVEREST->Close();
                   
			   	} // end loop simu iterations
               cout<<"out of loop over "<<SIMU_NTOT_iterations<<" simu iterations"<<endl;
		   	}// end loop over line types
            cout<<"out of loop over "<<NTYPES<<" types"<<endl;
		   	freeFromMemory(vectimes_Data);
            cout<<"closing datafile"<<endl;
		   	fdata->Close();	
            cout<<"datafile closed"<<endl;
		   //	}
	   	}// end loop runs
	    cout<<"out of loop over "<<Nbruns<<" runs"<<endl;


	
	}

	cout<<"END"<<endl;
	return 0;


}   


double CalculateDeltaT(double const& value,std::vector<double> const& vec)
{
	auto first=lower_bound(vec.begin(),vec.end(),value);
	
	double t1;
	double t2;
	double deltaT;
	if(first==vec.begin()){
		t2=*(first);
		deltaT=t2-value;
	}else if(first==vec.end()){
		t1=*(first-1);
		deltaT=t1-value;
	}else{
		t1=*(first-1);
		t2=*(first);
		if(abs(t1-value)<abs(t2-value)) {
			deltaT=t1-value;
		}else{
			deltaT=t2-value;
		}
	}

	return deltaT;

}






void LoadBinary(string const& binaryfilepath,string const& run,int const& RunPart,string const& Channel,std::vector <string> const& Detector_List,double *amplitudes,int const& detID,int const& NB_Samples_to_read_per_channel,int const& stampstart)
{

	const int NDET=(int)Detector_List.size();
	int dictionnary[NDET];
	

	//std::vector<int> vec_octet_offset;
	long int octet_offset;

	std::vector<string> vec_Channels;
	
	string format;
	if(RunPart<10){format="S0";}else{format="S";}
	
	// example of binary file ./CRYO_IPNL/rawdata/AC/2019/td16l001/td16l001_S02
	string binaryfile=Form("%s/%s/%s_%s%d",binaryfilepath.c_str(),run.c_str(),run.c_str(),format.c_str(),RunPart);
	
	/*
	// example of logfile file ./CRYO_IPNL/rawdata/AC/2019/td16l001/td16l001_log
	string logfile=Form("%s/%s/%s_log",binaryfilepath.c_str(),run.c_str(),run.c_str());
	
	/////// octet offset in log file in line : Les donnees commencent dans le fichier a l'octet 15374 (0x3c0e)
	ifstream fichierlog(logfile);
	if (!fichierlog){cout << " ==> ERROR: log file "<<logfile<<" does not exist !!" << endl;}
	string line;
	while(!fichierlog.eof()){
	
		getline(fichierlog,line);
		if((line.find("donnees", 0)) != string::npos && (line.find("commencent", 0)) != string::npos){
                string valOct = line.substr(line.find("octet"));
                int temp_octet_offset = std::stof(valOct.erase(0,6));
                vec_octet_offset.push_back(temp_octet_offset);
		}
		
	}
	fichierlog.close();
	*/
	
	///// find the different channels ///////
	ifstream fichier;
	fichier.open(binaryfile, ios::binary);
	if (!fichier){cout << " ==> ERROR: binary file "<<binaryfile<<" does not exist !!" << endl;}
 	string line;
 
 	int countdictionnary=-1;
 	while(1){
		
		getline(fichier,line);
		if((line.find("* Voie", 0)) != string::npos){
			countdictionnary+=1;
			vec_Channels.push_back(line);
			for(unsigned int k=0;k<Detector_List.size();k++){
				if(line.find(Channel,0) != string::npos && line.find(Detector_List[k],0) != string::npos){
					dictionnary[k]=countdictionnary;
				}
			}
			
		}
		

		if((line.find("* Donnees", 0)) != string::npos){
			octet_offset=fichier.tellg();
			break;
		} // Stop read before the binary part	
	}
 
 

	fichier.seekg(0,ios::end);  //vec_octet_offset[RunPart]
	//long int end = fichier.tellg();	// fichier.tellg() returns the number of bytes between offset and end of file
	
	int Nchannels=(int)vec_Channels.size();
	int single_length = sizeof(short)*Nchannels;
    short *buffer = new short[single_length/sizeof(short)];
	long int octet_start=octet_offset+(stampstart-NB_Samples_to_read_per_channel/2)*single_length;
   
	fichier.seekg(octet_start, std::ios::beg);

	for(int s=0;s<NB_Samples_to_read_per_channel;s++){
	
		fichier.read((char*)buffer,single_length);
		
		if(!memcmp(buffer,"* Arret ",8) || !memcmp(buffer,"* Arret ",3)){cout<<"End of Binary"<<endl; break;}
	
		amplitudes[s]=buffer[dictionnary[detID]];
	}

 	fichier.close();
}




