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
//#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
//#include "Math/GSLMinimizer.h"
//#include "Math/GSLSimAnMinimizer.h"
#include "TColor.h"
#include <X11/Xlib.h>
using namespace std;
string coincidencemethod="fast";

const int N_Max_Channel = 6;

template <typename VecType>
void freeFromMemory(VecType& myVec) {
     myVec.clear();
     VecType vert;
     myVec.swap(vert);
}

// convert anything into TString
template<typename T>
TString to_TStr(const T & aconvertir)
{
    ostringstream oss;
    oss << aconvertir;
    string sfromstream=oss.str();
    TString TSfroms=sfromstream;
    return TSfroms;
}


bool FileExists(std::string const & filename);
int color_rainbow(int const & i,int  const & nbval);
bool IsItaCoincidence(double const & t,std::vector<double> const & vec,double const & tolerance);
bool IsItaCoincidenceAndAntiCoincidence(double const & t,std::vector<double> const & veccoinc,std::vector<double> const & vecanticoinc,double const & tolerance);
double GetBinomialUncertainty(int const & ntrue,int const & nfalse);
void linspace(double const & min,double const & max,int const & N,std::vector<double> &vec);
void logspace(double const & min,double const & max,int const & N,std::vector<double> &vec);
bool AreTheseFilesDifferent(string  const & file1,string  const & file2);
double CalculateDeltaT(double const & value,std::vector<double> const & vec);
void buildefficiency(TH1D *h1,TH1D *h0,TH1D *he,int const & NBINS);
void LoadBinary(string const & binarypath,string const & run,int const & RunPart,string const & Channel,std::vector <string> const & Detector_List,double *amplitudes,int const &  detID,int const &  NB_Samples_to_read_per_channel,int const & stampstart);
void LoadFilteredBinary(string const &  fileandpath,double *amplitudes,int const &  NB_Samples_to_read_per_channel,int const & stampstart);
void itx(double const & x,double const & y,TString const &  mess,int const & color, double const & size,int  const & opt);
TCanvas* createcanvas(string const & namecanvas,string  const & style);
void BuildEfficiencyIterativeCuts(TGraphErrors *geff_Trigger,double *N_Trigger,double *N_all,double *Eline,int const &  Nlines,string const & scale);
void BuildEfficiencyIterativeCutsfromTH1D(TGraphErrors *geff_Trigger,TH1D *hcut,TH1D *hall);
int main(int argc, char *argv[])
{
	gStyle->SetOptStat(0);	
	TApplication *toto=new TApplication("toto",0,0);
	std::vector<string> vec_Directories;	
	vec_Directories.push_back("NbSi_chihalf_15V_Th");	
	const int Ndirectories=(int)vec_Directories.size();
	string masterdirectory="Efficiency_rought_Migdal";
	string masterpath=Form("/sps/edelweis/rootDataRun317/streams/%s",masterdirectory.c_str());	
	string DATA_Path="/sps/edelweis/rootDataRun317/streams/prodj";
	string study_index="vstudytestHugues";	
	bool onlycutRT=true;
	if(onlycutRT==true){study_index+="_RTonly";}	
	string det="NbSi209";	
	string Runs_to_analyse = argv[1];
	double Voltage = stod(argv[2]) ;
	double cutIon  = stod(argv[4]);
	int verbose = stoi(argv[3]);
	// not used anymore //
//-------------------------------------------------------------------------------------------------------------
	string binarypath="notusedanymore";
	string Filtered_binarypath="notusedanymore";
//-------------------------------------------------------------------------------------------------------------	
	
	
	string Figurespath="/sps/edelweis/rootDataRun317/streams/Selection10keV_Migdal/Figures";		
	string LegendcutRT_Maintenance;
	if(onlycutRT){
		LegendcutRT_Maintenance="Cut RT";
	}else{
		LegendcutRT_Maintenance="Cut RT and Maintenance";
	}
	
	// AJOUT POUR BINAIRES ///
//-------------------------------------------------------------------------------------------------------------
	std::vector<string> vecDet;
	vecDet.push_back(det);
	int detID=0;
	//////////////////////////
//-------------------------------------------------------------------------------------------------------------
    std::vector<string> vecruns ;
	std::cout << " opening run to process list : "<<Runs_to_analyse<<std::endl;
	ifstream Runs(Runs_to_analyse);
	std::string ListRun_name[100] = {""};
	int count_line_mc = 0 ;
	if(!Runs.fail()){
	    while(std::getline( Runs,ListRun_name[count_line_mc] )){
	        TString temp_str = ListRun_name[count_line_mc];
		    if(temp_str.Contains("#", TString::kIgnoreCase)) continue ;	
		    std::cout<<" List content "<< ListRun_name[count_line_mc]<<std::endl;
		    vecruns.push_back(ListRun_name[count_line_mc]);
		    count_line_mc++;
	    }
	}
	
	
	//std::vector<string> vecruns={"ud27a001","ud29a000","ud24a000","ud21a000","ud23a000","ud22a000"};//"tk26a000","tk19a000","tl15a001","tl22a000","tl18a000"};//,"tl24a000""tl15a001","tl15a002","tl24a000","tl18a000","tl22a000"};//Run_to_analyse};//"tl20a000","tl21a000"
   string template_2 = "Slow";

	const int Nbruns=(int)vecruns.size();
	std::vector<string> vec_Lines;
	
	/////////////////////////////// COUPURES 10KEV BANK /////////////////////////////////////
	//double EminBank=9.7;
	//double EmaxBank=11.;
	double EminBank=9.6;
	double EmaxBank=12.;
	
	/////////////////////////////// COUPURES CHI2 /////////////////////////////////////
	//channel A
	double par0_Normal=1.15;
	double par1_Normal=100.;
	double par2_Normal=300.;
	if(det == "NbSi209"){	
	   par0_Normal=1.2;
	   par1_Normal=100.;
	   par2_Normal=200.;
	
	}
	
	TF1 *fchi2_ION=new TF1("fchi2_ION","[0]",0,3e3);
	fchi2_ION->SetParameter(0,1.4);
	TF1 *fchi2_Normal=new TF1("fchi2_Normal","[0]+[1]*TMath::Power(x*((1+fabs([3])/3.))/[2],3)",0,18e3);
	fchi2_Normal->SetParameters(par0_Normal,par1_Normal,par2_Normal,Voltage);
	
	//channel B
	double par0B_Normal=1.15;
	double par1B_Normal=100.;
	double par2B_Normal=300.;
	if(det == "NbSi209"){	
	   par0B_Normal=1.2;
	   par1B_Normal=100.;
	   par2B_Normal=100.;
	
	}
	TF1 *fchi2B_Normal=new TF1("fchi2B_Normal","[0]+[1]*TMath::Power(x*((1+fabs([3])/3.))/[2],3)",0,18e3);
	fchi2B_Normal->SetParameters(par0B_Normal,par1B_Normal,par2B_Normal,Voltage);
	
	double par0_DELTA=0.001;//1000;//0.000;//0.001;
	double par1_DELTA=0.;
	double par2_DELTA=1.;
	TF1 *fchi2_DELTA=new TF1("fchi2_DELTA","[0]+[1]*TMath::Power(x/((1+fabs([3])/3.))/[2],3)",0.,1.5e3);
	fchi2_DELTA->SetParameters(par0_DELTA,par1_DELTA,par2_DELTA,Voltage);

	// ici
	double SIMU_Emin=0;
	double SIMU_Emax=2.935;
	string SIMU_Type="Line"; // or line
	int SIMU_NTOT_iterations_MAX=2000; // will take the min of chosen and available NTOT_iterations;
	int DATA_NFiles_MAX=2000; // will take the min of chosen and available
	
	int fracmult_N_samples=1;
	
	bool useEverest=true;
	
	double conversionkeV2ADU=780; // ATTENTION !!!! c'est une valeur approximative pour runs de May
	if(useEverest==true){
		conversionkeV2ADU=1; // no conversion when Everest becauze already calibrated
	}
	
	
	string BankLineNameBegin=Form("Everestized_Chi2Cut_%.4f",par0_DELTA);

	////////////////// ATTENTION QUENTIN ENLEVE IF(1) a la place de PAssCurRt//////////////////////////
	bool saveLines=true;	
	bool LoadTracesfromBinaries=false;
	bool showtraces=false;
	bool usechihalf = true;
	int N_samples_Heat=1024; //  N_samples_Heat to show when show traces  otherwise useless
	int ndofchi2 = 1024 ;
	if(usechihalf) ndofchi2=512;
	/////// below, it only matters if showtraces==true ///////
	bool showfiltered=false;
	double plotE_In_min=0.01; // only show events with E_In>plotE_In_min		
	bool show=false;	// if show=true show all
	bool showcoinc_simu_datasimu=false;
	//////////////////////////////////////////////////////////	
	string common_name="";
	if(Nbruns==1){
		common_name=vecruns[0];
	}else{
		common_name="";
		for(int r=0;r<Nbruns;r++){
			if(r==0){
				common_name+=vecruns[r];
			}else{
				common_name+="-"+vecruns[r];
			}
		}
	}	
	string rootfilesaveLines="/sps/edelweis/rootDataRun317/streams/Selection10keV_Migdal/Bank_Lines/"+study_index+"_"+BankLineNameBegin+"_"+common_name+"_Bank_Lines.root";	
	if(SIMU_Type!="Line"){saveLines=false;}
	if(verbose > 0) cout<<rootfilesaveLines<<endl;		
	if(SIMU_Type=="Line"){	
///////////////// LOOK FOR LINES ////////
	   string commandline=Form("ls -d %s/%s/SimuCoinc/*/*/Line* | awk -F\"/Line_\" '{print $2}' | awk -F\"_\" '{print $1}' | sort -n | uniq > templist_Lines.txt",masterpath.c_str(),vec_Directories.at(0).c_str());	
	   if(verbose > 0) cout<<commandline<<endl;
	   //system(commandline.c_str());
	   ifstream fline;
	   fline.open("templist_Lines.txt");	
	   string buffer;
	   while(1){
		   fline>>buffer;
		   if(!fline.good()){
		      break;
		   }
		   vec_Lines.push_back(buffer);
	   }	
	   if(verbose > 0) cout<<vec_Lines.size()<< "Lines found"<<endl;
	   for(int i=0;i<(int)vec_Lines.size();i++){
		   if(verbose > 0) cout<<i<<"	"<<vec_Lines[i]<<endl;
	   }
	}	
	int NFILES_BUG=0;
	std::vector<string> filesbugged;
	bool optionlogmonitoring=false;	
	double ANALYSIS_Emax=1499/1000.;
	double ANALYSIS_Emin=0.0010;	
	if(SIMU_Type!="Line"){
		if(ANALYSIS_Emin<SIMU_Emin || ANALYSIS_Emax>SIMU_Emax){
	      if(verbose > 0) cout<<"Change the Analysis Energy range"<<endl; 
	      exit(1);
	   }
	}
	double ANALYSIS_Ewidth=0.007;
	const int Energy_NBINS=(int)((ANALYSIS_Emax-ANALYSIS_Emin)/ANALYSIS_Ewidth);
	if(floor(((ANALYSIS_Emax-ANALYSIS_Emin)/ANALYSIS_Ewidth))!=((ANALYSIS_Emax-ANALYSIS_Emin)/ANALYSIS_Ewidth)){
	   if(verbose > 0) cout<<"Choose a bin width that gives rise to (int) Nbins here : "<<((ANALYSIS_Emax-ANALYSIS_Emin)/ANALYSIS_Ewidth)<< " condition  "<<floor(((ANALYSIS_Emax-ANALYSIS_Emin)/ANALYSIS_Ewidth))<<endl; 
	   //exit(1);
	}
	if(verbose > 0) cout<<"Energy_NBINS="<<Energy_NBINS<<endl;
	
	std::vector<double> vectolerances;
	double tolmin=2;
	double tolmax=5;
	int Ntol=1;
	bool logtol=false;
	if(logtol==true){
	   logspace(tolmin,tolmax,Ntol,vectolerances);
	}else{
	   linspace(tolmin,tolmax,Ntol,vectolerances);
	}
	const int N_tolerances=(int)vectolerances.size();
	TCanvas* ctest=createcanvas("ctest","implement");	
	if(showfiltered){ctest->Divide(1,2);}
	ctest->cd(1);
	
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////////////////// 		PROGRAM 	//////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
	
 	TFile *foutsaveLines;
	if(saveLines){
		foutsaveLines=new TFile(rootfilesaveLines.c_str(),"RECREATE");
		foutsaveLines->cd();
	}	
	TH1D *h1_simu_coinc_data[Ntol];
	TH1D *h0_simu_coinc_data[Ntol];
	TH1D *he_simu_coinc_data[Ntol];	
		
	TH1D *h1_simu_coinc_datasimu[Ntol];
	TH1D *h0_simu_coinc_datasimu[Ntol];
	TH1D *he_simu_coinc_datasimu[Ntol];
	
	TH1D *CorrectedE_h1_simu_coinc_datasimu_anticoinc_data[Ntol];
	TH1D *CorrectedE_h0_simu_coinc_datasimu_anticoinc_data[Ntol];
	TH1D *CorrectedE_he_simu_coinc_datasimu_anticoinc_data[Ntol];
	
	TH1D *CorrectedE_h20_simu_coinc_datasimu_anticoinc_data_cut1[Ntol];
	TH1D *CorrectedE_h2_simu_coinc_datasimu_anticoinc_data_cut1[Ntol];
	TH1D *CorrectedE_h2eff_simu_coinc_datasimu_anticoinc_data_cut1[Ntol];

	TH1D *h1_simu_coinc_datasimu_anticoinc_data[Ntol];
	TH1D *h0_simu_coinc_datasimu_anticoinc_data[Ntol];
	TH1D *he_simu_coinc_datasimu_anticoinc_data[Ntol];
	
	TH1D *h0_simu_coinc_datasimu_anticoinc_data_cut1[Ntol];
	TH1D *h1_simu_coinc_datasimu_anticoinc_data_cut1[Ntol];
	TH1D *he_simu_coinc_datasimu_anticoinc_data_cut1[Ntol];
		
	// 2D Plots 					   	
	double Emin2D=ANALYSIS_Emin;
	double Emax2D=ANALYSIS_Emax*conversionkeV2ADU;	
	double NORMALChi2min2D=0.3;
	double NORMALChi2max2D=10;
	int Nbins2D=1000;	
	double DELTAChi2min2D=-0.5;
	double DELTAChi2max2D=0.5;
	if(ANALYSIS_Emax<(100./1000.)){
		DELTAChi2min2D=-0.05;
		DELTAChi2max2D=0.05;
		NORMALChi2max2D=2;
	}	
	int DELTANbins2D=2000;	
	Double_t xEdges[Nbins2D + 1];
  	Double_t NORMALyEdges[Nbins2D + 1]; 
	double NORMALlog_min=TMath::Log10(NORMALChi2min2D);
    double NORMALlog_max=TMath::Log10(NORMALChi2max2D);
    double NORMALlog_delta=(NORMALlog_max-NORMALlog_min)/Nbins2D;   	    	
	double E_width=(Emax2D-Emin2D)/Nbins2D;
	for(int i=0;i<=Nbins2D;i++){
   		xEdges[i]=Emin2D+i*E_width;
   		NORMALyEdges[i]=TMath::Power(10,NORMALlog_min+i*NORMALlog_delta);
  	}	
	TH2D *h2D_2_NORMALchi2vsE[Ntol]; // coinc and pass the cut
	TH2D *h2D_1_NORMALchi2vsE[Ntol]; // coinc but doesnt pass the cut 
	TH2D *h2D_0_NORMALchi2vsE[Ntol]; // not in coinc	
	TH2D *h2D_2_NORMALchi2BvsE[Ntol]; // coinc and pass the cut
	TH2D *h2D_1_NORMALchi2BvsE[Ntol]; // coinc but doesnt pass the cut 
	TH2D *h2D_0_NORMALchi2BvsE[Ntol];	
	TH2D *h2D_2_DELTAchi2vsE[Ntol]; // coinc and pass the cut
	TH2D *h2D_1_DELTAchi2vsE[Ntol]; // coinc but doesnt pass the cut 
	TH2D *h2D_0_DELTAchi2vsE[Ntol]; // not in coinc	
	TH2D *h2D_2_EoutvsE[Ntol]; // coinc and pass the cut
	TH2D *h2D_1_EoutvsE[Ntol]; // coinc but doesnt pass the cut 
	TH2D *h2D_0_EoutvsE[Ntol]; // not in coinc	
	int Nsupposedlines=(int)vec_Lines.size();
	int Nlines;
	if(SIMU_Type!="Line"){
		Nlines=1;
	}else{
		Nlines=Nsupposedlines;
	}	
	TH1D *h1EoutLine_2Normalized_EoutvsE[Ntol][Nlines]; // all coinc and pass the cut
	TH1D *h1EoutLine_2_EoutvsE[Ntol][Nlines]; // all coinc and pass the cut
	TH1D *h1EoutLine_1_EoutvsE[Ntol][Nlines]; // all coinc
	TH1D *h1EoutLine_0_EoutvsE[Ntol][Nlines]; // all
	
	TH1D *h1EoutLine_2Normalized_EoutvsEB[Ntol][Nlines]; // all coinc and pass the cut
	TH1D *h1EoutLine_2_EoutvsEB[Ntol][Nlines]; // all coinc and pass the cut
	TH1D *h1EoutLine_1_EoutvsEB[Ntol][Nlines]; // all coinc
	TH1D *h1EoutLine_0_EoutvsEB[Ntol][Nlines];
	
	TH1D *h1Chi2Line_2_EoutvsChi2[Ntol][Nlines][6]; // all coinc and pass the cut
	TH1D *h1Chi2Line_1_EoutvsChi2[Ntol][Nlines][6]; // all coinc
	TH1D *h1Chi2Line_0_EoutvsChi2[Ntol][Nlines][6];
	
	//Eout chal vs EionB -> check reco 0 -> EionB 1 -> EionFid
	TH1D *h1Line_2_EoutvsEion[Ntol][Nlines][2]; // all coinc and pass the cut
	TH1D *h1Line_1_EoutvsEion[Ntol][Nlines][2]; // all coinc
	TH1D *h1Line_0_EoutvsEion[Ntol][Nlines][2];
	
	
	TGraph* chi2_vs_Ein_2 [6];
	TGraph* chi2_vs_Ein_1 [6];
	TGraph* chi2_vs_Ein_0 [6];
	TGraph* Ea_vs_Ein_0 = new TGraph();
	TGraph* Eb_vs_Ein_0 = new TGraph();
	TGraph* Ea_vs_Ein_1 = new TGraph();
	TGraph* Eb_vs_Ein_1 = new TGraph();
	TGraph* Ea_vs_Ein_2 = new TGraph();
	TGraph* Eb_vs_Ein_2 = new TGraph();
	
	TGraph* Et_vs_EionB_2 = new TGraph();
	TGraph* Et_vs_EionB_1 = new TGraph();
	
	TGraph* Et_vs_EionFid_2 = new TGraph();
	TGraph* Et_vs_EionFid_1 = new TGraph();
	
	
	TGraph* sigEa_vs_Ein_1 = new TGraph();
	TGraph* sigEb_vs_Ein_1 = new TGraph();
	TGraph* sigEa_vs_Ein_2 = new TGraph();
	TGraph* sigEb_vs_Ein_2 = new TGraph();
	TGraph* sigchi2_vs_Ein_2 [6];
	TGraph* sigchi2_vs_Ein_1 [6];
	double Eminout=1./1000.; // 0.1 eV
	double Emaxout=1.5; // 13 keV
	double chi2min = 0.1;
	double chi2max = 10;
	int Nbinsout=(Emaxout-Eminout)/(0.1/1000.);		
	for(int tol_index=0;tol_index<N_tolerances;tol_index++){
		for(int l=0;l<Nlines;l++){
			h1EoutLine_2Normalized_EoutvsE[tol_index][l]=new TH1D(Form("h1EoutLine_2Normalized_EoutvsE%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			h1EoutLine_2_EoutvsE[tol_index][l]=new TH1D(Form("h1EoutLine_2_EoutvsE%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			h1EoutLine_1_EoutvsE[tol_index][l]=new TH1D(Form("h1EoutLine_1_EoutvsE%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			h1EoutLine_0_EoutvsE[tol_index][l]=new TH1D(Form("h1EoutLine_0_EoutvsE%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			h1EoutLine_2Normalized_EoutvsEB[tol_index][l]=new TH1D(Form("h1EoutLine_2Normalized_EoutvsE%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			h1EoutLine_2_EoutvsEB[tol_index][l]=new TH1D(Form("h1EoutLine_2_EoutvsEB%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			h1EoutLine_1_EoutvsEB[tol_index][l]=new TH1D(Form("h1EoutLine_1_EoutvsEB%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			h1EoutLine_0_EoutvsEB[tol_index][l]=new TH1D(Form("h1EoutLine_0_EoutvsEB%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			
			h1Line_2_EoutvsEion[tol_index][l][0]=new TH1D(Form("h1Line_2_EoutvsEionB%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			h1Line_1_EoutvsEion[tol_index][l][0]=new TH1D(Form("h1Line_1_EoutvsEionB%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			h1Line_0_EoutvsEion[tol_index][l][0]=new TH1D(Form("h1Line_1_EoutvsEionB%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			h1Line_2_EoutvsEion[tol_index][l][1]=new TH1D(Form("h1Line_2_EoutvsEionFid%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			h1Line_1_EoutvsEion[tol_index][l][1]=new TH1D(Form("h1Line_1_EoutvsEionFid%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			h1Line_0_EoutvsEion[tol_index][l][1]=new TH1D(Form("h1Line_1_EoutvsEionFid%i_%i",tol_index,l),"",Nbinsout,Eminout,Emaxout);
			
			
			
			for(int itchi = 0 ; itchi < 6 ; itchi++){
			   h1Chi2Line_2_EoutvsChi2[tol_index][l][itchi]=new TH1D(Form("h1Chi2Line_2_EoutvsChi2%i_%i_%i",tol_index,l,itchi),"",Nbinsout,chi2min,chi2max);
			   h1Chi2Line_1_EoutvsChi2[tol_index][l][itchi]=new TH1D(Form("h1Chi2Line_1_EoutvsChi2%i_%i_%i",tol_index,l,itchi),"",Nbinsout,chi2min,chi2max);
			   h1Chi2Line_0_EoutvsChi2[tol_index][l][itchi]=new TH1D(Form("h1Chi2Line_0_EoutvsChi2%i_%i_%i",tol_index,l,itchi),"",Nbinsout,chi2min,chi2max);
			   chi2_vs_Ein_2[itchi]= new TGraph();
			   chi2_vs_Ein_1[itchi]= new TGraph();
			   chi2_vs_Ein_0[itchi]= new TGraph();
			   
			   sigchi2_vs_Ein_2[itchi]= new TGraph();
			   sigchi2_vs_Ein_1[itchi]= new TGraph();
			}
			if(h1EoutLine_2Normalized_EoutvsE[tol_index][l]->GetBinWidth(1)!=(0.0001)){if(verbose > 0) cout<<"binning is not 0.1 eV"<<endl; exit(1);}
		}
	}
	double N_all[Nlines]; // all 
	double N_Trigger[Nlines]; // coinc simu anti coinc data
	double N_cut_Bank[Nlines]; // coinc simu anti coinc data + bank selection
	double N_cut_RT_Maintenance[Nlines]; // etc.. 
	double N_cut_Chi2NormIon[Nlines]; // etc.. 
	double N_cut_Chi2Norm[Nlines]; // etc
	double N_cut_Chi2NormB[Nlines];
	double N_cut_DeltaChi2[Nlines];
	double N_cut_Eionization[Nlines];
	double Ncut_FidChal[Nlines];
	double Eline[Nlines];	
	int NtrialsDeltaChi2 = 0;
	double N_cut_TrialDeltaChi2[NtrialsDeltaChi2][Nlines];	
	int ultimatecheckntotentries=0;
	for(int l=0;l<Nlines;l++){
		N_all[l]=0;
		N_Trigger[l]=0;
		N_cut_Bank[l]=0;
		N_cut_RT_Maintenance[l]=0;
		N_cut_Chi2NormIon[l]=0;
		N_cut_Chi2Norm[l] =0;
		N_cut_Eionization[l]=0;
		N_cut_Chi2NormB[l]=0;
		N_cut_DeltaChi2[l]=0;
		Ncut_FidChal[0];
		Eline[l]=(double)atof(vec_Lines[l].c_str());
		if(verbose > 0) cout<<"check conversion into double : "<<Eline[l]<<endl;		
		for(int c=0;c<NtrialsDeltaChi2;c++){
		   N_cut_TrialDeltaChi2[c][l]=0;
		}
	}			
	double eV_Emin=0.001;
	double eV_Emax=1499.;
	int eV_Nbins=215;
	TH1D *eV_hEout_all            =new TH1D("eV_hEout_all","",eV_Nbins,eV_Emin,eV_Emax);
	TH1D *eV_hEout_Trigger        =new TH1D("eV_hEout_Trigger","",eV_Nbins,eV_Emin,eV_Emax);
	TH1D *eV_hEout_Bank           =new TH1D("eV_hEout_Bank","",eV_Nbins,eV_Emin,eV_Emax);
	TH1D *eV_hEout_RT_Maintenance =new TH1D("eV_hEout_RT_Maintenance","",eV_Nbins,eV_Emin,eV_Emax);
	TH1D *eV_hEout_Chi2NormIon    =new TH1D("eV_hEout_Chi2NormIon","",eV_Nbins,eV_Emin,eV_Emax);
	TH1D *eV_hEout_chalfid        =new TH1D("eV_hEout_chalfid","",eV_Nbins,eV_Emin,eV_Emax);
	TH1D *eV_hEout_Eion           =new TH1D("eV_hEout_Eion","",eV_Nbins,eV_Emin,eV_Emax);
	TH1D *eV_hEout_Chi2Norm       =new TH1D("eV_hEout_Chi2Norm","",eV_Nbins,eV_Emin,eV_Emax);
	TH1D *eV_hEout_Chi2NormB      =new TH1D("eV_hEout_Chi2NormB","",eV_Nbins,eV_Emin,eV_Emax);
	TH1D *eV_hEout_DeltaChi2      =new TH1D("eV_hEout_DeltaChi2","",eV_Nbins,eV_Emin,eV_Emax);	
	TH1D *eV_hEout_Trial_DeltaChi2[NtrialsDeltaChi2];
	for(int c=0;c<NtrialsDeltaChi2;c++){
		eV_hEout_Trial_DeltaChi2[c]=new TH1D(Form("eV_hEout_Trial_DeltaChi2_%i",c),"",eV_Nbins,eV_Emin,eV_Emax);
	}
/////////////////////////LOOP TOLERANCE////////////////////////////////////////////////////////////////////////
   for(int tol_index=0;tol_index<N_tolerances;tol_index++){	
	double tolerance=vectolerances[tol_index];
	if(verbose > 0) cout<<"tolerance = "<<tolerance<<endl;		
	h1_simu_coinc_data[tol_index]    =new TH1D(Form("h1_simu_coinc_data_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	h0_simu_coinc_data[tol_index]    =new TH1D(Form("h0_simu_coinc_data_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	he_simu_coinc_data[tol_index]    =new TH1D(Form("he_simu_coinc_data_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);	
	h1_simu_coinc_datasimu[tol_index]=new TH1D(Form("h1_simu_coinc_datasimu_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	h0_simu_coinc_datasimu[tol_index]=new TH1D(Form("h0_simu_coinc_datasimu_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	he_simu_coinc_datasimu[tol_index]=new TH1D(Form("he_simu_coinc_datasimu_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);	
	CorrectedE_h1_simu_coinc_datasimu_anticoinc_data[tol_index]         =new TH1D(Form("CorrectedE_h1_simu_coinc_datasimu_anticoinc_data_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	CorrectedE_h0_simu_coinc_datasimu_anticoinc_data[tol_index]         =new TH1D(Form("CorrectedE_h0_simu_coinc_datasimu_anticoinc_data_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	CorrectedE_he_simu_coinc_datasimu_anticoinc_data[tol_index]         =new TH1D(Form("CorrectedE_he_simu_coinc_datasimu_anticoinc_data_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	CorrectedE_h2_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]    =new TH1D(Form("CorrectedE_h2_simu_coinc_datasimu_anticoinc_data_cut1__tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	CorrectedE_h20_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]   =new TH1D(Form("CorrectedE_h20_simu_coinc_datasimu_anticoinc_data_cut1__tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	CorrectedE_h2eff_simu_coinc_datasimu_anticoinc_data_cut1[tol_index] =new TH1D(Form("CorrectedE_h2eff_simu_coinc_datasimu_anticoinc_data_cut1_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	
	h1_simu_coinc_datasimu_anticoinc_data[tol_index]     =new TH1D(Form("h1_simu_coinc_datasimu_anticoinc_data_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	h0_simu_coinc_datasimu_anticoinc_data[tol_index]     =new TH1D(Form("h0_simu_coinc_datasimu_anticoinc_data_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	he_simu_coinc_datasimu_anticoinc_data[tol_index]     =new TH1D(Form("he_simu_coinc_datasimu_anticoinc_data_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);	
	h0_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]=new TH1D(Form("h0_simu_coinc_datasimu_anticoinc_data_cut1_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	h1_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]=new TH1D(Form("h1_simu_coinc_datasimu_anticoinc_data_cut1_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	he_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]=new TH1D(Form("he_simu_coinc_datasimu_anticoinc_data_cut1_tol%i",tol_index),"",Energy_NBINS,ANALYSIS_Emin,ANALYSIS_Emax);
	h2D_2_NORMALchi2vsE[tol_index] =new TH2D(Form("h2D_2_NORMALchi2vsE_tol%i",tol_index) ,"",Nbins2D,xEdges,Nbins2D,NORMALyEdges);
	h2D_1_NORMALchi2vsE[tol_index] =new TH2D(Form("h2D_1_NORMALchi2vsE_tol%i",tol_index) ,"",Nbins2D,xEdges,Nbins2D,NORMALyEdges);
	h2D_0_NORMALchi2vsE[tol_index] =new TH2D(Form("h2D_0_NORMALchi2vsE_tol%i",tol_index) ,"",Nbins2D,xEdges,Nbins2D,NORMALyEdges);
	h2D_2_NORMALchi2BvsE[tol_index]=new TH2D(Form("h2D_2_NORMALchi2BvsE_tol%i",tol_index),"",Nbins2D,xEdges,Nbins2D,NORMALyEdges);
	h2D_1_NORMALchi2BvsE[tol_index]=new TH2D(Form("h2D_1_NORMALchi2BvsE_tol%i",tol_index),"",Nbins2D,xEdges,Nbins2D,NORMALyEdges);
	h2D_0_NORMALchi2BvsE[tol_index]=new TH2D(Form("h2D_0_NORMALchi2BvsE_tol%i",tol_index),"",Nbins2D,xEdges,Nbins2D,NORMALyEdges);	
	h2D_2_DELTAchi2vsE[tol_index]  =new TH2D(Form("h2D_2_IONchi2vsE_tol%i",tol_index)    ,"",Nbins2D,xEdges,Nbins2D,NORMALyEdges);
	h2D_1_DELTAchi2vsE[tol_index]  =new TH2D(Form("h2D_1_IONchi2vsE_tol%i",tol_index)    ,"",Nbins2D,xEdges,Nbins2D,NORMALyEdges);
	h2D_0_DELTAchi2vsE[tol_index]  =new TH2D(Form("h2D_0_IONchi2vsE_tol%i",tol_index)    ,"",Nbins2D,xEdges,Nbins2D,NORMALyEdges);	
	h2D_2_EoutvsE[tol_index]       =new TH2D(Form("h2D_2_EoutvsE_tol%i",tol_index)       ,"",Nbins2D,ANALYSIS_Emin,ANALYSIS_Emax,Nbins2D,Emin2D,Emax2D);
	h2D_1_EoutvsE[tol_index]       =new TH2D(Form("h2D_1_EoutvsE_tol%i",tol_index)       ,"",Nbins2D,ANALYSIS_Emin,ANALYSIS_Emax,Nbins2D,Emin2D,Emax2D);
	h2D_0_EoutvsE[tol_index]       =new TH2D(Form("h2D_0_EoutvsE_tol%i",tol_index)       ,"",Nbins2D,ANALYSIS_Emin,ANALYSIS_Emax,Nbins2D,Emin2D,Emax2D);
   int triggerchanel = 0 ;
   if(det=="NbSi209") triggerchanel=1;
	for(int dir=0;dir<Ndirectories;dir++){	
	   string directory=vec_Directories[dir];		
	   string SIMUCOINC_Path=Form("%s/%s/SimuCoinc",masterpath.c_str(),directory.c_str());
	   if(verbose > 0) cout<<"before loop over runs"<<endl;
	   if(verbose > 0) cout<<SIMUCOINC_Path<<endl;
	   for(int r=0;r<Nbruns;r++){	
		   string run=vecruns[r];
		   if(verbose > 0) cout<<"run "<<run<<endl;	
//--------------------------CHANGE HERE TO PAIRS PARTITION-----------------------------------------------------------------------------------------------------		   	
		   string allrootfiles=Form("%s/%s/%s/ProcessedData_%s_S*_%s_ChanTrig%i.root",DATA_Path.c_str(),run.c_str(),det.c_str(),run.c_str(),det.c_str(),triggerchanel);		
		   string command=Form("ls %s | wc -l",allrootfiles.c_str());
		int nfiles_TOT = gSystem->GetFromPipe(command.c_str()).Atoi(); // the bash command returns a TString that is converted into int with Atoi();
	   	if(verbose > 0) cout<<"nfiles found :"<<nfiles_TOT<<endl;
	   	if(nfiles_TOT==0){
	   	   if(verbose > 0) cout<<"file does not exist"<<endl;
	   	   exit(1); 
	   	}
		   nfiles_TOT=min(DATA_NFiles_MAX,nfiles_TOT);
		   if(verbose > 0) cout<<"nfiles used :"<<nfiles_TOT<<endl;
		   for(int filenumber=0;filenumber<nfiles_TOT;filenumber+=2){		
	   		string format;
		   	if(filenumber<10){
		   	   format="S0";
		   	}else{
		   	   format="S";
		      }
			   double Trace_Data_at_t0[N_samples_Heat];
			   double Trace_DataAndSimu_Filtered[N_samples_Heat];
			   double Trace_Data_Filtered[N_samples_Heat];
			   double times[N_samples_Heat];
			   for(int i=0;i<N_samples_Heat;i++){ 
			      times[i]=i;
			   }						
/////////////////////////// INIT SIMU //////////////////////////////////// 
// determine the number of iterations			
			   int NTYPES=1;
			   if(SIMU_Type=="Line"){
			      NTYPES=(int)vec_Lines.size();
			   }
			
			   for(int t=0;t<NTYPES;t++){
			
				   string pattern;
				   if(SIMU_Type=="Flat"){
				      pattern=Form("%s/%s/%s/%s_%.4f_%.4f_*",SIMUCOINC_Path.c_str(),run.c_str(),det.c_str(),SIMU_Type.c_str(),SIMU_Emin,SIMU_Emax);
				   }else if(SIMU_Type=="Line"){
				      pattern=Form("%s/%s/%s/%s_%s_*",SIMUCOINC_Path.c_str(),run.c_str(),det.c_str(),SIMU_Type.c_str(),vec_Lines[t].c_str());
				   }else{
				      if(verbose > 0) cout<<"this simutype doesn't exist"<<endl;
				      exit(1);
				   }
				   string command=Form("ls -d %s -R 2>/dev/null | wc -l",pattern.c_str());
		     		int niterations_temp = gSystem->GetFromPipe(command.c_str()).Atoi(); // the bash command returns a TString that is converted into int with Atoi();   		
				   int SIMU_NTOT_iterations=min(niterations_temp,SIMU_NTOT_iterations_MAX);				
				   int YES=0;
				   int NO=0;
				   for (int SIMU_iteration=1;SIMU_iteration<=SIMU_NTOT_iterations;SIMU_iteration++){
					   if(verbose > 0 /*&& t%100 == 0*/) cout<<"Simu iteration "<<SIMU_iteration<<" over "<<SIMU_NTOT_iterations<<endl;					
					   if(verbose > 0 /*&& t%100 == 0*/) cout<<Form("directory %s , run %d/%d , file %d/%d , Type %d/%d , iteration %d/%d",directory.c_str(),r,Nbruns,filenumber,nfiles_TOT,t,NTYPES,SIMU_iteration,SIMU_NTOT_iterations)<<endl;
					
					   string SIMU_condition;
					   if(SIMU_Type=="Flat"){
					      SIMU_condition=Form("%s_%.4f_%.4f_%i",SIMU_Type.c_str(),SIMU_Emin,SIMU_Emax,SIMU_iteration);
					   }else if(SIMU_Type=="Line"){
					      SIMU_condition=Form("%s_%s_%i",SIMU_Type.c_str(),vec_Lines[t].c_str(),SIMU_iteration);
					   }else{
					      if(verbose > 0) cout<<"this simutype doesn't exist"<<endl; 
					      exit(1);
					   }
					string SIMU_File=Form("ProcessedData_%s_%s%i_%s_%s_ChanTrig%i.root",run.c_str(),format.c_str(),filenumber,det.c_str(),SIMU_condition.c_str(),triggerchanel);		//S*[02468]	
					//string SIMU_File=Form("ProcessedData_%s_S*[02468]_%s_%s_ChanTrig%i.root",run.c_str(),format.c_str(),filenumber,det.c_str(),SIMU_condition.c_str(),triggerchanel);	
				   	string SIMUCOINC_File=Form("SIMUCOINC_%s",SIMU_File.c_str());
				   	string SIMUCOINC_RELATIVE_Path=Form("%s/%s/%s/%s",SIMUCOINC_Path.c_str(),run.c_str(),det.c_str(),SIMU_condition.c_str());
				   	string SIMUCOINC_PathAndFile=SIMUCOINC_RELATIVE_Path+"/"+SIMUCOINC_File;						
////////////////////////////////// LOAD SIMU COINC ///////////////////////////////////////////
				   	if(verbose > 0 && t%100 == 0 ) cout<<"load simu coinc"<<endl;				 
				   	if(verbose > 0 /*&& t%100 == 0*/) cout<<SIMUCOINC_PathAndFile<<endl;
				   	TFile *fcoinc=new TFile(SIMUCOINC_PathAndFile.c_str(),"READ");				   	
				   	TTree *tcoinc = (TTree*)fcoinc->Get("tsimu");
				   	TTree *tRuntree = (TTree*)fcoinc->Get("RunTree_Normal");
                    double deltaT_simu_datasimu;
				   	double deltaT_simu_data;
				   	double Energy_OF_Normal[N_Max_Channel];
				   	double Energy_OFInput[N_Max_Channel];
				   	double chi2_OF_Normal[N_Max_Channel];
				   	double chi2_OF_saved_norm[N_Max_Channel];
				   	double chi2_OF_Fast[N_Max_Channel];
				   	double Delta_chi2_Normalized;
				   	double Energy_In;
				   	int t0_timestamp_SimuOnly;
				   	int Pass_Cut_Maintenance;
				   	double Time_modulo_reset;
				   	int MegaStp;
					double EhA;
					double EhB;
					double EiA;
					double EiB;
					double EiC;
					double EiD;
					double Eh;
					Int_t R_t;					
				   	tcoinc->SetBranchAddress("deltaT_simu_datasimu",&deltaT_simu_datasimu);
				   	tcoinc->SetBranchAddress("deltaT_simu_data",&deltaT_simu_data);				  
				   	tcoinc->SetBranchAddress("Energy_OF_Normal",&Energy_OF_Normal);
				   	tRuntree->SetBranchAddress("Rt",&R_t);
				   	if(useEverest==true){
				   		int test=tcoinc->SetBranchAddress("EhA",&EhA);
				   		test=tcoinc->SetBranchAddress("EhB",&EhB);
				   		test=tcoinc->SetBranchAddress("Eh",&Eh);
				   		test=tcoinc->SetBranchAddress("EiA",&EiA);
				   		test=tcoinc->SetBranchAddress("EiB",&EiB);
				   		test=tcoinc->SetBranchAddress("EiC",&EiC);
				   		test=tcoinc->SetBranchAddress("EiD",&EiD);				   		
				   	}
				   	tcoinc->SetBranchAddress("Pass_Cut_Maintenance",&Pass_Cut_Maintenance);
				   	tcoinc->SetBranchAddress("Time_modulo_reset",&Time_modulo_reset);				   	
				   	if(usechihalf){
				   	   tcoinc->SetBranchAddress("chi2_OF_half",&chi2_OF_Normal);
				   	   tcoinc->SetBranchAddress("chi2_OF_Normal",&chi2_OF_saved_norm);
				   	}else{
				   	   tcoinc->SetBranchAddress("chi2_OF_Normal",&chi2_OF_Normal);
				   	}
				   	tcoinc->SetBranchAddress("Energy_OFInput",&Energy_OFInput);
				   	tcoinc->SetBranchAddress(Form("chi2_OF_%s",template_2.c_str()),&chi2_OF_Fast);
				   	tcoinc->SetBranchAddress("Delta_chi2_Normalized",&Delta_chi2_Normalized);
				   	tcoinc->SetBranchAddress("Energy_In",&Energy_In);
		  			   tcoinc->SetBranchAddress("t0_timestamp",&t0_timestamp_SimuOnly);
		  			   if(LoadTracesfromBinaries==true){
		  			      tcoinc->SetBranchAddress("Trace_Data_at_t0",&Trace_Data_at_t0);
		  			   }
//////////////////////////////// LOAD REAL DATA AND SIMU /////////////////////////////////////////
			   	   int Nentries_SimuCoinc=tcoinc->GetEntries();
				   	ultimatecheckntotentries+=Nentries_SimuCoinc;
				      if(verbose > 0 && t%100 == 0) cout<<"loop over coinc "<<endl;
				   	for(int entry=0;entry<Nentries_SimuCoinc;entry++){				   		
				   		tcoinc->GetEntry(entry);			   		
				   		bool coinc_simu_datasimu=(bool)(abs(deltaT_simu_datasimu)<tolerance);
				   		bool coinc_simu_data=(bool)(abs(deltaT_simu_data)<tolerance);
				   		bool coinc_simu_datasimu_anticoinc_data;
				   		if(coinc_simu_datasimu && !coinc_simu_data){
				   		   coinc_simu_datasimu_anticoinc_data=true;
				   		}else{
				   		   coinc_simu_datasimu_anticoinc_data=false;
				   		}
//////////////////////////////////////////////////////////////////////////////////////////////////
				 		   double energy_in_kev;
				 		   double energy_in_kevB;
				 		   double energy_in_kevTot;
				 		   double Lukefactor = (1. +  fabs(Voltage)/3.);
				 		   
				 		   double energy_in_kevphonon;
				 		   double energy_in_kevphononB;
				 		   double energy_in_kevphononTot;
				 		   
				 		   double energy_ion_A = EiA;
				 		   double energy_ion_B = EiB;
				 		   double energy_ion_C = EiC;
				 		   double energy_ion_D = EiD;
				 		   double EiFid ;
				 		   
				 		   if(det== "NbSi209") EiFid   = (EiA + EiB)/2.;
				 		   
				 		   if(useEverest==true){
				 			   energy_in_kev   =EhA;
				 			   energy_in_kevB  =EhB;
				 			   energy_in_kevTot=Eh;
				 			   
				 			   /*energy_in_kevphonon   =EhA*(1+fabs(Voltage)/3.));
				 			   energy_in_kevphononB  =EhB*(1+fabs(Voltage)/3.));
				 			   energy_in_kevphononTot=Eh*(1+fabs(Voltage)/3.));*/
				 		   }else{
				 			   energy_in_kev=Energy_OF_Normal[0]/conversionkeV2ADU;
				 		   }				 		
				 		   double energy_in_ev   =energy_in_kev*1000;
				 		   double energy_in_evB  =energy_in_kevB*1000;
				 		   double energy_in_evTot=energy_in_kevTot*1000;				 		
				 		   double Normal_chi2normalized   =     chi2_OF_saved_norm[0]/1024.;
				 		   double Normal_chi2normalizedB  =     chi2_OF_Normal[1]/ndofchi2;	
				 		   double Normal_chi2normalizedIonA =	chi2_OF_saved_norm[2]/1024;		 		
				 		   double Normal_chi2normalizedIonB =	chi2_OF_saved_norm[3]/1024;
				 		   if(usechihalf){
				 		      Normal_chi2normalizedIonA =	chi2_OF_saved_norm[2]/1024.;		 		
				 		      Normal_chi2normalizedIonB =	chi2_OF_saved_norm[3]/1024.;
				 		   }		 		
				 		   bool passcut10keVBank=false;
				 		   
						   if(Energy_OFInput[0]>=EminBank && Energy_OFInput[0]<=EmaxBank){
						      //std::if(verbose > 0) cout<< " Pass cut bank "<<Energy_OFInput[0]<<std::endl;
						      passcut10keVBank=true;
						   }else{
						      passcut10keVBank=false;
						      //std::if(verbose > 0) cout<< " testing cut bank "<<Energy_OFInput[0]<<std::endl;
						   }				 		
				 		   bool passNormalchi2cut=false;
				 		   bool passIonisationchi2cut = false ;
				 		   if(energy_in_kev>0){
				 		      if(Normal_chi2normalizedIonA <= 1.4 && Normal_chi2normalizedIonB <= 1.4){
				 		         passIonisationchi2cut=true;
				 		         //std::if(verbose > 0) cout<<" passed chi2 ioni A "<< Normal_chi2normalizedIonA<<" B "<< Normal_chi2normalizedIonB << std::endl;
				 		      }else{
				 		         passIonisationchi2cut=false;
				 		         //std::if(verbose > 0) cout<<" failed chi2 ioni A "<< Normal_chi2normalizedIonA<<" B "<< Normal_chi2normalizedIonB << std::endl;
				 		      }
				 			   if(Normal_chi2normalized<fchi2_Normal->Eval(energy_in_kev)){
				 			         passNormalchi2cut=true;
				 			      }else{
				 			         passNormalchi2cut=false;
				 			      }
				 		   }else{
				 			   if(Normal_chi2normalized<fchi2_Normal->Eval(0.)){
				 			         passNormalchi2cut=true;
				 			      }else{
				 			         passNormalchi2cut=false;
				 			      }
				 		   }
				 		   bool passEionCut = false;
				 		   if(EiFid > cutIon || cutIon < 0.){
				 		    //std::cout<< " test ion EB "<<energy_ion_B <<std::endl;
				 		   	passEionCut = true ;
				 		   }else{
				 		   	passEionCut = false;			 		   
				 		   }
				 		   bool passFidchalCut = false;
				 		   if(pow(energy_in_kev*Lukefactor - energy_in_kevB*Lukefactor,2) < pow(0.6,2) + pow(0.045*energy_in_kevTot*Lukefactor,2)   ){
				 		    passFidchalCut=true;
				 		   }else{
				 		    passFidchalCut=false;
				 		   }
				 		   			 		   
				 		   bool passNormalchi2Bcut=false;
				 		   if(energy_in_kevB>0){
				 			   if(Normal_chi2normalizedB<fchi2B_Normal->Eval(energy_in_kevB)){				 			   
				 			         passNormalchi2Bcut=true;
				 			         
				 			      }else{
				 			         passNormalchi2Bcut=false;
				 			      }
				 		   }else{
				 			   if(Normal_chi2normalizedB<fchi2B_Normal->Eval(0.)){
				 			         passNormalchi2Bcut=true;
				 			      }else{
				 			         passNormalchi2Bcut=false;
				 			      }
				 		   }									
						   bool passDELTAchi2cut=false;
						   if(energy_in_kev>0){
							   if(Delta_chi2_Normalized<fchi2_DELTA->Eval(energy_in_kev)){
							      passDELTAchi2cut=true;
							   }else{
							      passDELTAchi2cut=false;
							   }
						   }else{
							   if(Delta_chi2_Normalized<fchi2_DELTA->Eval(0.)){
							      passDELTAchi2cut=true;
							   }else{
							      passDELTAchi2cut=false;
							   }
						   }									
						   bool passMaintenance=false;
						   if(Pass_Cut_Maintenance==1){
						      passMaintenance=true;
						   }else{
						      passMaintenance=false;
						   }			
						   bool passRT=false;
						   if(Time_modulo_reset>2.5 && Time_modulo_reset<(128-0.5)){
						      passRT=true;
						   }else{
						      //std::if(verbose > 0) cout<<" R_t "<<R_t<<std::endl; 
						      passRT=false;
						   }			
						   bool passRTandMaintenance=false;
						   if(passMaintenance && passRT){
						      passRTandMaintenance=true;
						   }else{
						      passRTandMaintenance=false;						  
						   }
						   bool passALLcut=false;					
						   if(onlycutRT==true && passRT==true){
						      passRTandMaintenance=true;
						   } // AJOUT DEBUG			
						   if(det == "RED30")  passNormalchi2Bcut = true; 
						   if(det== "NbSi209") passDELTAchi2cut   = true; 
						   if(passNormalchi2cut && passNormalchi2Bcut && passDELTAchi2cut && passRTandMaintenance && passcut10keVBank && passIonisationchi2cut && passFidchalCut){
						      passALLcut=true; 
						   }else{
						      passALLcut=false;
						   }
            	 		if(tol_index==0){
							   N_all[t]=N_all[t]+1;
							   eV_hEout_all->Fill(energy_in_ev);
							   if(passRTandMaintenance){
								   N_cut_RT_Maintenance[t]++;
								   eV_hEout_RT_Maintenance->Fill(energy_in_ev);								
								   if(passcut10keVBank){								
									   N_cut_Bank[t]++;
									   eV_hEout_Bank->Fill(energy_in_ev);									
									   if(coinc_simu_datasimu_anticoinc_data){									
										   N_Trigger[t]++;
										   eV_hEout_Trigger->Fill(energy_in_ev);	
										   if(passIonisationchi2cut){
										    N_cut_Chi2NormIon[t]++;
										   	eV_hEout_Chi2NormIon->Fill(energy_in_ev);	
										    if(passEionCut){
										        N_cut_Eionization[t]++;
										        eV_hEout_Chi2NormIon->Fill(EiFid);
										        if(passFidchalCut){
										        Ncut_FidChal[t]++;
										        eV_hEout_chalfid->Fill(energy_in_ev);
										        if(passNormalchi2cut){										
											         N_cut_Chi2Norm[t]++;
											        eV_hEout_Chi2Norm->Fill(energy_in_ev);
											         if(passNormalchi2Bcut){									
											             N_cut_Chi2NormB[t]++;
											             eV_hEout_Chi2NormB->Fill(energy_in_evB);
											             
											             if(passDELTAchi2cut){								
												            N_cut_DeltaChi2[t]++;
												            eV_hEout_DeltaChi2->Fill(energy_in_ev);
											            }											
											            for(int c=0;c<NtrialsDeltaChi2;c++){
												             if(/*Delta_chi2_Normalized<vecTrialDeltaChi2[c]*/true){
													             N_cut_TrialDeltaChi2[c][t]++;
													             eV_hEout_Trial_DeltaChi2[c]->Fill(energy_in_ev);
												             }
											         }//loop delta chicut
											      }//chi2B
										      }//chi2A
										      }//chal fid cut
										     }//Eionozation 
										   }//chiIon
									   }//trigger
								   }//bank cut
							   }//RT
						   }//loop over tolerance			 		
				 		   // histo without any cut to normalize Line spectra
				 		   h1EoutLine_0_EoutvsE[tol_index][t]->Fill(energy_in_kev);		
				 		   h1EoutLine_0_EoutvsEB[tol_index][t]->Fill(energy_in_kevB);
				 		   h1Line_0_EoutvsEion[tol_index][t][0]->Fill(energy_ion_B);
						   h1Line_0_EoutvsEion[tol_index][t][1]->Fill((energy_ion_B+energy_ion_A)/2.);
				 		   for(int itchi = 0 ; itchi < 6 ; itchi++){
								   double chi2 = chi2_OF_saved_norm[itchi] / 1024.;
								      if( itchi <= 1){
								         //NDOF = 512 ;
								         chi2 = chi2_OF_Normal[itchi]/ 512 ;
								       }
								      h1Chi2Line_0_EoutvsChi2[tol_index][t][itchi]->Fill(chi2);
								   }						 		
						   if(coinc_simu_datasimu_anticoinc_data){							
				 			   h1EoutLine_1_EoutvsE[tol_index][t]->Fill(energy_in_kev);
				 			   h1EoutLine_1_EoutvsEB[tol_index][t]->Fill(energy_in_kevB);
				 			   h1Line_1_EoutvsEion[tol_index][t][0]->Fill(energy_ion_B);
							   h1Line_1_EoutvsEion[tol_index][t][1]->Fill((energy_ion_B+energy_ion_A)/2.);
				 			   for(int itchi = 0 ; itchi < 6 ; itchi++){
								   double chi2 = chi2_OF_saved_norm[itchi] / 1024.;
								      if( itchi <= 1){
								         //NDOF = 512 ;
								         chi2 = chi2_OF_Normal[itchi]/ 512 ;
								       }
								      h1Chi2Line_1_EoutvsChi2[tol_index][t][itchi]->Fill(chi2);
								   }						 			
							   if(passALLcut){
								   h1EoutLine_2_EoutvsE[tol_index][t]->Fill(energy_in_kev);
								   h1EoutLine_2Normalized_EoutvsE[tol_index][t]->Fill(energy_in_kev);
								   h1EoutLine_2_EoutvsEB[tol_index][t]->Fill(energy_in_kevB);
								   h1EoutLine_2Normalized_EoutvsEB[tol_index][t]->Fill(energy_in_kevB);
								   h1Line_2_EoutvsEion[tol_index][t][0]->Fill(energy_ion_B);
								   h1Line_2_EoutvsEion[tol_index][t][1]->Fill((energy_ion_B+energy_ion_A)/2.);
								   
								   for(int itchi = 0 ; itchi < 6 ; itchi++){
								   double NDOF = 1024.;
								      double chi2 = chi2_OF_saved_norm[itchi] / 1024.;
								      if( itchi <= 1){
								         //NDOF = 512 ;
								         chi2 = chi2_OF_Normal[itchi]/ 512 ;
								       }
								      h1Chi2Line_2_EoutvsChi2[tol_index][t][itchi]->Fill(chi2);								      
								   }
							   }
						   }
////////////////////////////////////////
						   if(coinc_simu_data){
							   h1_simu_coinc_data[tol_index]->Fill(Energy_In);					   	
						   }else{
							   h0_simu_coinc_data[tol_index]->Fill(Energy_In);
						   }
////////////////////////////////////////					
						   if(coinc_simu_datasimu){
							   h1_simu_coinc_datasimu[tol_index]->Fill(Energy_In);					
						   }else{
							   h0_simu_coinc_datasimu[tol_index]->Fill(Energy_In);
						   }
///////////// ADDITIONAL CUTS ONLY FOR SIMU COINC DATA SIMU ANTICOINC DATA//////////////////								
						   double energychaleur ;
						   double energychaleurB;
						   if(useEverest==true){
						   	   energychaleur =EhA;
							   energychaleurB=EhB;
						   }else{
							   energychaleur =Energy_OF_Normal[0];
							   energychaleurB=Energy_OF_Normal[1];
						   }								
						   if(coinc_simu_datasimu_anticoinc_data){
							   h1_simu_coinc_datasimu_anticoinc_data[tol_index]->Fill(Energy_In);					
							   if(passALLcut){
								   h2D_2_NORMALchi2vsE [tol_index]->Fill(energychaleur ,Normal_chi2normalized);
								   h2D_2_NORMALchi2BvsE[tol_index]->Fill(energychaleurB,Normal_chi2normalizedB);
								   h2D_2_DELTAchi2vsE[tol_index]  ->Fill(energy_ion_B,Normal_chi2normalizedIonB);
								   h2D_2_EoutvsE[tol_index]       ->Fill(Energy_In,energychaleur);
							   }else{
								   h2D_1_NORMALchi2vsE [tol_index]->Fill(energychaleur ,Normal_chi2normalized);
								   h2D_1_NORMALchi2BvsE[tol_index]->Fill(energychaleurB,Normal_chi2normalizedB);
								   h2D_1_DELTAchi2vsE[tol_index]  ->Fill(energy_ion_B,Normal_chi2normalizedIonB);
								   h2D_1_EoutvsE[tol_index]       ->Fill(Energy_In,energychaleur);
							   }									
						   }else{
							   h0_simu_coinc_datasimu_anticoinc_data[tol_index]->Fill(Energy_In);
							   h2D_0_NORMALchi2vsE[tol_index] ->Fill(energychaleur ,Normal_chi2normalized);
							   h2D_0_NORMALchi2BvsE[tol_index]->Fill(energychaleurB,Normal_chi2normalizedB);
							   h2D_0_DELTAchi2vsE[tol_index]  ->Fill(energy_ion_B,Normal_chi2normalizedIonB);
							   h2D_0_EoutvsE[tol_index]       ->Fill(Energy_In,energychaleur);
						   }								
						   if(coinc_simu_datasimu_anticoinc_data && passALLcut){
							   h1_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]->Fill(Energy_In);
							   CorrectedE_h2_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]->Fill(energy_in_kev);
						   }else{
							   h0_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]->Fill(Energy_In);
							   CorrectedE_h20_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]->Fill(energy_in_kev);
						   }					
						   if(coinc_simu_datasimu_anticoinc_data){
							   CorrectedE_h1_simu_coinc_datasimu_anticoinc_data[tol_index]->Fill(energy_in_kev);
						   }else{
							   CorrectedE_h0_simu_coinc_datasimu_anticoinc_data[tol_index]->Fill(energy_in_kev);	
						   }															 	
				   	}//endloop over coincidence	   	
                  fcoinc->Close();				
			   	}//endloop simu iteration		   	
		   	}// endloop line types iteration (for flat there is only 1 iteration)	
	   	}//endloop filenumber	
	   }//endloop run	
   }//endloop directories	
}// endloop tolerance
if(verbose > 0) cout<<"calculate"<<endl;
///////////////////////// CALCULATE EFFICIENCIES /////////////////////// 
for(int tol_index=0;tol_index<N_tolerances;tol_index++){	
	buildefficiency(h1_simu_coinc_data[tol_index],h0_simu_coinc_data[tol_index],he_simu_coinc_data[tol_index],Energy_NBINS);
	he_simu_coinc_data[tol_index]->SetMaximum(1.05);
	he_simu_coinc_data[tol_index]->SetMinimum(-0.05);
	he_simu_coinc_data[tol_index]->SetLineWidth(2);
	
	buildefficiency(h1_simu_coinc_datasimu[tol_index],h0_simu_coinc_datasimu[tol_index],he_simu_coinc_datasimu[tol_index],Energy_NBINS);
	he_simu_coinc_datasimu[tol_index]->SetMaximum(1.1);
	he_simu_coinc_datasimu[tol_index]->SetMinimum(-0.050);
	he_simu_coinc_datasimu[tol_index]->SetLineWidth(2);
	
	buildefficiency(h1_simu_coinc_datasimu_anticoinc_data[tol_index],h0_simu_coinc_datasimu_anticoinc_data[tol_index],he_simu_coinc_datasimu_anticoinc_data[tol_index],Energy_NBINS);
	he_simu_coinc_datasimu_anticoinc_data[tol_index]->SetMaximum(1.05);
	he_simu_coinc_datasimu_anticoinc_data[tol_index]->SetMinimum(-0.05);
	he_simu_coinc_datasimu_anticoinc_data[tol_index]->SetLineWidth(2);
	buildefficiency(h1_simu_coinc_datasimu_anticoinc_data_cut1[tol_index],h0_simu_coinc_datasimu_anticoinc_data_cut1[tol_index],he_simu_coinc_datasimu_anticoinc_data_cut1[tol_index],Energy_NBINS);
	he_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]->SetMaximum(1.05);
	he_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]->SetMinimum(-0.050);
	he_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]->SetLineWidth(2);
	buildefficiency(CorrectedE_h1_simu_coinc_datasimu_anticoinc_data[tol_index],CorrectedE_h0_simu_coinc_datasimu_anticoinc_data[tol_index],CorrectedE_he_simu_coinc_datasimu_anticoinc_data[tol_index],Energy_NBINS);
	CorrectedE_he_simu_coinc_datasimu_anticoinc_data[tol_index]->SetMaximum(1.05);
	CorrectedE_he_simu_coinc_datasimu_anticoinc_data[tol_index]->SetMinimum(-0.05);
	CorrectedE_he_simu_coinc_datasimu_anticoinc_data[tol_index]->SetLineWidth(2);	
	buildefficiency(CorrectedE_h2_simu_coinc_datasimu_anticoinc_data_cut1[tol_index],CorrectedE_h20_simu_coinc_datasimu_anticoinc_data_cut1[tol_index],CorrectedE_h2eff_simu_coinc_datasimu_anticoinc_data_cut1[tol_index],Energy_NBINS);
	CorrectedE_h2eff_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]->SetMaximum(1.05);
	CorrectedE_h2eff_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]->SetMinimum(-0.05);
	CorrectedE_h2eff_simu_coinc_datasimu_anticoinc_data_cut1[tol_index]->SetLineWidth(2);
}//endloop tolerance

///// CALCULATE EFFICIENCIES ITERATIVE CUTS ///////////////////
TGraphErrors *eV_geff_hEout_all=new TGraphErrors();
TGraphErrors *eV_geff_hEout_RT_Maintenance=new TGraphErrors();
TGraphErrors *eV_geff_hEout_Bank=new TGraphErrors();
TGraphErrors *eV_geff_hEout_Trigger=new TGraphErrors();
TGraphErrors *eV_geff_hEout_Chi2NormION=new TGraphErrors();
TGraphErrors *eV_geff_hEout_Chi2Norm=new TGraphErrors();
TGraphErrors *eV_geff_hEout_Chi2BNorm=new TGraphErrors();
TGraphErrors *eV_geff_hEout_DeltaChi2=new TGraphErrors();

BuildEfficiencyIterativeCutsfromTH1D(eV_geff_hEout_RT_Maintenance,eV_hEout_RT_Maintenance,eV_hEout_all);
BuildEfficiencyIterativeCutsfromTH1D(eV_geff_hEout_Bank,eV_hEout_Bank,eV_hEout_all);
BuildEfficiencyIterativeCutsfromTH1D(eV_geff_hEout_Trigger,eV_hEout_Trigger,eV_hEout_all);
BuildEfficiencyIterativeCutsfromTH1D(eV_geff_hEout_Chi2NormION ,eV_hEout_Chi2NormIon ,eV_hEout_all);
BuildEfficiencyIterativeCutsfromTH1D(eV_geff_hEout_Chi2Norm ,eV_hEout_Chi2Norm ,eV_hEout_all);
BuildEfficiencyIterativeCutsfromTH1D(eV_geff_hEout_Chi2BNorm,eV_hEout_Chi2NormB,eV_hEout_all);
BuildEfficiencyIterativeCutsfromTH1D(eV_geff_hEout_DeltaChi2,eV_hEout_DeltaChi2,eV_hEout_all);

TGraphErrors *eV_geffrelative_hEout_Bank_over_RT_Maintenance=new TGraphErrors();
TGraphErrors *eV_geffrelative_hEout_Trigger_over_Bank=new TGraphErrors();
TGraphErrors *eV_geffrelative_hEout_Chi2Norm_over_Trigger=new TGraphErrors();
TGraphErrors *eV_geffrelative_hEout_Chi2BNorm_over_Trigger=new TGraphErrors();
TGraphErrors *eV_geffrelative_hEout_DeltaChi2_over_Chi2Norm=new TGraphErrors();
TGraphErrors *eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm[NtrialsDeltaChi2];
for(int c=0;c<NtrialsDeltaChi2;c++){
	eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm[c]=new TGraphErrors();
	//BuildEfficiencyIterativeCutsfromTH1D(eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm[c],eV_hEout_Trial_DeltaChi2[c],eV_hEout_Chi2Norm);
}
//BuildEfficiencyIterativeCutsfromTH1D(eV_geffrelative_hEout_Bank_over_RT_Maintenance,eV_hEout_Bank,eV_hEout_RT_Maintenance);
//BuildEfficiencyIterativeCutsfromTH1D(eV_geffrelative_hEout_Trigger_over_Bank,eV_hEout_Trigger,eV_hEout_Bank);
//BuildEfficiencyIterativeCutsfromTH1D(eV_geffrelative_hEout_Chi2Norm_over_Trigger,eV_hEout_Chi2Norm,eV_hEout_Trigger);
//BuildEfficiencyIterativeCutsfromTH1D(eV_geffrelative_hEout_Chi2BNorm_over_Trigger,eV_hEout_Chi2NormB,eV_hEout_Trigger);
//BuildEfficiencyIterativeCutsfromTH1D(eV_geffrelative_hEout_DeltaChi2_over_Chi2Norm,eV_hEout_DeltaChi2,eV_hEout_Chi2Norm);
TGraphErrors *Trial_eV_geffrelative_cut_DeltaChi2[NtrialsDeltaChi2];
for(int c=0;c<NtrialsDeltaChi2;c++){
	Trial_eV_geffrelative_cut_DeltaChi2[c]=new TGraphErrors();
	BuildEfficiencyIterativeCuts(Trial_eV_geffrelative_cut_DeltaChi2[c],N_cut_TrialDeltaChi2[c],N_cut_Chi2Norm,Eline,Nlines,"eV");	
}
TGraphErrors *eV_geff_cut_RT_Maintenance   =new TGraphErrors();
TGraphErrors *eV_geff_cut_Bank             =new TGraphErrors();
TGraphErrors *eV_geff_Trigger              =new TGraphErrors(); 
TGraphErrors *eV_geff_cut_Chi2NormION      =new TGraphErrors();
TGraphErrors *eV_geff_cut_Eion             =new TGraphErrors();
TGraphErrors *eV_geff_cut_chalfid         =new TGraphErrors();
TGraphErrors *eV_geff_cut_Chi2Norm         =new TGraphErrors();
TGraphErrors *eV_geff_cut_Chi2BNorm        =new TGraphErrors();
TGraphErrors *eV_geff_cut_DeltaChi2        =new TGraphErrors();

TGraphErrors *keV_geff_cut_RT_Maintenance    =new TGraphErrors();
TGraphErrors *keV_geff_cut_Bank              =new TGraphErrors();
TGraphErrors *keV_geff_Trigger               =new TGraphErrors(); 
TGraphErrors *keV_geff_cut_Chi2NormION       =new TGraphErrors();
TGraphErrors *keV_geff_cut_Eion              =new TGraphErrors();
TGraphErrors *keV_geff_cut_chalfid           =new TGraphErrors();
TGraphErrors *keV_geff_cut_Chi2Norm          =new TGraphErrors();
TGraphErrors *keV_geff_cut_Chi2BNorm         =new TGraphErrors();
TGraphErrors *keV_geff_cut_DeltaChi2         =new TGraphErrors();

BuildEfficiencyIterativeCuts(eV_geff_Trigger,N_Trigger,N_all,Eline,Nlines,"eV");
BuildEfficiencyIterativeCuts(eV_geff_cut_Bank,N_cut_Bank,N_all,Eline,Nlines,"eV");
BuildEfficiencyIterativeCuts(eV_geff_cut_RT_Maintenance,N_cut_RT_Maintenance,N_all,Eline,Nlines,"eV");
BuildEfficiencyIterativeCuts(eV_geff_cut_Chi2NormION,N_cut_Chi2NormIon,N_all,Eline,Nlines,"eV");
BuildEfficiencyIterativeCuts(eV_geff_cut_Eion,N_cut_Eionization,N_all,Eline,Nlines,"eV");
BuildEfficiencyIterativeCuts(eV_geff_cut_chalfid,Ncut_FidChal,N_all,Eline,Nlines,"eV");
BuildEfficiencyIterativeCuts(eV_geff_cut_Chi2Norm,N_cut_Chi2Norm,N_all,Eline,Nlines,"eV");
BuildEfficiencyIterativeCuts(eV_geff_cut_Chi2BNorm,N_cut_Chi2NormB,N_all,Eline,Nlines,"eV");
BuildEfficiencyIterativeCuts(eV_geff_cut_DeltaChi2,N_cut_DeltaChi2,N_all,Eline,Nlines,"eV");

BuildEfficiencyIterativeCuts(keV_geff_Trigger,N_Trigger,N_all,Eline,Nlines,"keV");
BuildEfficiencyIterativeCuts(keV_geff_cut_Bank,N_cut_Bank,N_all,Eline,Nlines,"keV");
BuildEfficiencyIterativeCuts(keV_geff_cut_RT_Maintenance,N_cut_RT_Maintenance,N_all,Eline,Nlines,"keV");
BuildEfficiencyIterativeCuts(keV_geff_cut_Chi2NormION ,N_cut_Chi2NormIon,N_all,Eline,Nlines,"keV");
BuildEfficiencyIterativeCuts(keV_geff_cut_Eion,N_cut_Eionization,N_all,Eline,Nlines,"keV");
BuildEfficiencyIterativeCuts(keV_geff_cut_chalfid,Ncut_FidChal,N_all,Eline,Nlines,"keV");
BuildEfficiencyIterativeCuts(keV_geff_cut_Chi2Norm ,N_cut_Chi2Norm,N_all,Eline,Nlines,"keV");
BuildEfficiencyIterativeCuts(keV_geff_cut_Chi2BNorm,N_cut_Chi2NormB,N_all,Eline,Nlines,"keV");
BuildEfficiencyIterativeCuts(keV_geff_cut_DeltaChi2,N_cut_DeltaChi2,N_all,Eline,Nlines,"keV");

TGraphErrors *eV_geffrelative_Bank_over_RT_Maintenance=new TGraphErrors();
TGraphErrors *eV_geffrelative_Trigger_over_Bank=new TGraphErrors();
TGraphErrors *eV_geffrelative_Chi2Norm_over_Trigger=new TGraphErrors();
TGraphErrors *eV_geffrelative_Chi2BNorm_over_Trigger=new TGraphErrors();
TGraphErrors *eV_geffrelative_DeltaChi2_over_Chi2Norm=new TGraphErrors();
/*
BuildEfficiencyIterativeCuts(eV_geffrelative_Bank_over_RT_Maintenance,N_cut_Bank,N_cut_RT_Maintenance,Eline,Nlines,"eV");
BuildEfficiencyIterativeCuts(eV_geffrelative_Trigger_over_Bank,N_Trigger,N_cut_Bank,Eline,Nlines,"eV");
BuildEfficiencyIterativeCuts(eV_geffrelative_Chi2Norm_over_Trigger,N_cut_Chi2Norm,N_Trigger,Eline,Nlines,"eV");
BuildEfficiencyIterativeCuts(eV_geffrelative_Chi2BNorm_over_Trigger,N_cut_Chi2NormB,N_Trigger,Eline,Nlines,"eV");
BuildEfficiencyIterativeCuts(eV_geffrelative_DeltaChi2_over_Chi2Norm,N_cut_DeltaChi2,N_cut_Chi2Norm,Eline,Nlines,"eV");*/

TGraphErrors *keV_geffrelative_Bank_over_RT_Maintenance=new TGraphErrors();
TGraphErrors *keV_geffrelative_Trigger_over_Bank=new TGraphErrors();
TGraphErrors *keV_geffrelative_Chi2Norm_over_Trigger=new TGraphErrors();
TGraphErrors *keV_geffrelative_Chi2BNorm_over_Trigger=new TGraphErrors();
TGraphErrors *keV_geffrelative_DeltaChi2_over_Chi2Norm=new TGraphErrors();

/*BuildEfficiencyIterativeCuts(keV_geffrelative_Bank_over_RT_Maintenance,N_cut_Bank,N_cut_RT_Maintenance,Eline,Nlines,"keV");
BuildEfficiencyIterativeCuts(keV_geffrelative_Trigger_over_Bank,N_Trigger,N_cut_Bank,Eline,Nlines,"keV");
BuildEfficiencyIterativeCuts(keV_geffrelative_Chi2Norm_over_Trigger,N_cut_Chi2Norm,N_Trigger,Eline,Nlines,"keV");
BuildEfficiencyIterativeCuts(keV_geffrelative_Chi2BNorm_over_Trigger,N_cut_Chi2NormB,N_Trigger,Eline,Nlines,"keV");
BuildEfficiencyIterativeCuts(keV_geffrelative_DeltaChi2_over_Chi2Norm,N_cut_DeltaChi2,N_cut_Chi2Norm,Eline,Nlines,"keV");*/

////////////////////////////////////////////////////////////////////////////////////////////////////
/*TCanvas *canvas_eV_geff_hEout=new TCanvas("canvas_eV_geff_hEout","canvas_eV_geff_hEout",67,55,1853,1025);
TH2D *h2_eV_geff_hEout=new TH2D("h2_eV_geff_hEout",Form("%s;Energy Reconstructed [eV];Efficiency",common_name.c_str()),1,10.,1500.,1,-0.01,1.01);
h2_eV_geff_hEout->DrawCopy();

eV_geff_hEout_RT_Maintenance->SetLineColor(color_rainbow(0,6));
eV_geff_hEout_RT_Maintenance->SetMarkerColor(color_rainbow(0,6));
eV_geff_hEout_RT_Maintenance->Draw("LPsame");
eV_geff_hEout_Bank->SetLineColor(color_rainbow(1,6));
eV_geff_hEout_Bank->SetMarkerColor(color_rainbow(1,6));
eV_geff_hEout_Bank->Draw("LPsame");
eV_geff_hEout_Trigger->SetLineColor(color_rainbow(2,6));
eV_geff_hEout_Trigger->SetMarkerColor(color_rainbow(2,6));
eV_geff_hEout_Trigger->Draw("LPsame");
eV_geff_hEout_Chi2NormION->SetLineColor(color_rainbow(5,6));
eV_geff_hEout_Chi2NormION->SetMarkerColor(color_rainbow(5,6));
eV_geff_hEout_Chi2NormION->Draw("LPsame");
eV_geff_hEout_Chi2Norm->SetLineColor(color_rainbow(3,6));
eV_geff_hEout_Chi2Norm->SetMarkerColor(color_rainbow(3,6));
eV_geff_hEout_Chi2Norm->Draw("LPsame");
eV_geff_hEout_Chi2BNorm->SetLineColor(color_rainbow(5,6));
eV_geff_hEout_Chi2BNorm->SetMarkerColor(color_rainbow(5,6));
eV_geff_hEout_Chi2BNorm->Draw("LPsame");
eV_geff_hEout_DeltaChi2->SetLineColor(color_rainbow(4,6));
eV_geff_hEout_DeltaChi2->SetMarkerColor(color_rainbow(4,6));
//eV_geff_hEout_DeltaChi2->Draw("LPsame");
gPad->SetGridx();
gPad->SetGridy();
TLegend *leg_eV_geff_hEout=new TLegend(0.5,0.1,0.9,0.4);
leg_eV_geff_hEout->AddEntry(eV_geff_hEout_RT_Maintenance,LegendcutRT_Maintenance.c_str(),"lep");
leg_eV_geff_hEout->AddEntry(eV_geff_hEout_Bank,"Cut 10 keV Bank","lep");
leg_eV_geff_hEout->AddEntry(eV_geff_hEout_Trigger,"Trigger (coinc Simu anticoinc Data)","lep");
leg_eV_geff_hEout->AddEntry(eV_geff_hEout_Chi2NormION,"Cut Chi2 ION","lep");
leg_eV_geff_hEout->AddEntry(eV_geff_hEout_Chi2Norm,"Cut Chi2 Normal","lep");
leg_eV_geff_hEout->AddEntry(eV_geff_hEout_Chi2BNorm,"Cut Chi2 B Normal","lep");
//leg_eV_geff_hEout->AddEntry(eV_geff_hEout_DeltaChi2,Form("Cut Delta Chi2 < %.3f" ,par0_DELTA),"lep");
leg_eV_geff_hEout->DrawClone();
/*
TCanvas *canvas_eV_geffrelative_hEout=new TCanvas("canvas_eV_geffrelative_hEout","canvas_eV_geffrelative_hEout",67,55,1853,1025);
TH2D *h2_eV_geffrelative_hEout=new TH2D("h2_eV_geffrelative_hEout",Form("%s;Energy Reconstructed [eV];Relative Efficiency",common_name.c_str()),1,10.,1500.,1,-0.01,1.01);
h2_eV_geffrelative_hEout->DrawCopy();

eV_geffrelative_Bank_over_RT_Maintenance->SetLineColor(color_rainbow(1,6));
eV_geffrelative_Bank_over_RT_Maintenance->SetMarkerColor(color_rainbow(1,6));
eV_geffrelative_Bank_over_RT_Maintenance->Draw("LPsame");
eV_geffrelative_hEout_Trigger_over_Bank->SetLineColor(color_rainbow(2,6));
eV_geffrelative_hEout_Trigger_over_Bank->SetMarkerColor(color_rainbow(2,6));
eV_geffrelative_hEout_Trigger_over_Bank->Draw("LPsame");
eV_geffrelative_hEout_Chi2Norm_over_Trigger->SetLineColor(color_rainbow(3,6));
eV_geffrelative_hEout_Chi2Norm_over_Trigger->SetMarkerColor(color_rainbow(3,6));
eV_geffrelative_hEout_Chi2Norm_over_Trigger->Draw("LPsame");
eV_geffrelative_hEout_Chi2BNorm_over_Trigger->SetLineColor(color_rainbow(5,6));
eV_geffrelative_hEout_Chi2BNorm_over_Trigger->SetMarkerColor(color_rainbow(5,6));
eV_geffrelative_hEout_Chi2BNorm_over_Trigger->Draw("LPsame");
eV_geffrelative_hEout_DeltaChi2_over_Chi2Norm->SetLineColor(color_rainbow(4,6));
eV_geffrelative_hEout_DeltaChi2_over_Chi2Norm->SetMarkerColor(color_rainbow(4,6));
//eV_geffrelative_hEout_DeltaChi2_over_Chi2Norm->Draw("LPsame");
gPad->SetGridx();
gPad->SetGridy();
TLegend *leg_eV_geffrelative_hEout=new TLegend(0.5,0.1,0.9,0.4);
leg_eV_geffrelative_hEout->AddEntry(eV_geffrelative_Bank_over_RT_Maintenance,"Cut 10 keV Bank","lep");
leg_eV_geffrelative_hEout->AddEntry(eV_geffrelative_hEout_Trigger_over_Bank,"Trigger (coinc Simu anticoinc Data)","lep");
leg_eV_geffrelative_hEout->AddEntry(eV_geffrelative_hEout_Chi2Norm_over_Trigger,"Cut Chi2 Normal","lep");
leg_eV_geffrelative_hEout->AddEntry(eV_geffrelative_hEout_Chi2BNorm_over_Trigger,"Cut Chi2 B Normal","lep");
//leg_eV_geffrelative_hEout->AddEntry(eV_geffrelative_hEout_DeltaChi2_over_Chi2Norm,Form("Cut Delta Chi2 < %.3f" ,par0_DELTA),"lep");
leg_eV_geffrelative_hEout->DrawClone();
canvas_eV_geffrelative_hEout->Update();
canvas_eV_geffrelative_hEout->SaveAs("effrelative_eV_Eout.png");*/
/////////////////////////////

//TCanvas *canvas_eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm=new TCanvas("canvas_eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm","canvas_eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm",67,55,1853,1025);
TH2D *h2_eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm=new TH2D("h2_eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm",Form("%s;Energy Reconstructed [eV];Relative Efficiency",common_name.c_str()),1,10.,1500.,1,-0.01,1.01);
//h2_eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm->DrawCopy();
TLegend *leg_eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm=new TLegend(0.7,0.1,0.9,0.5);
for(int c=0;c<NtrialsDeltaChi2;c++){
	eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm[c]->SetLineColor(color_rainbow(NtrialsDeltaChi2-1-c,NtrialsDeltaChi2));
	eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm[c]->SetMarkerColor(color_rainbow(NtrialsDeltaChi2-1-c,NtrialsDeltaChi2));
	//eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm[c]->Draw("LPsame");
	//leg_eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm->AddEntry(eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm[c],Form("Delta Chi2 < %.4f",vecTrialDeltaChi2[c]),"lep");
}
//leg_eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm->DrawClone("same");
//canvas_eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm->Update();
//canvas_eV_geffrelative_hEout_Trial_DeltaChi2_over_Chi2Norm->SaveAs("Trial_effrelative_eV_Eout.png");
//gPad->SetGridx();
//gPad->SetGridy();
//////////////////////////////////////////
//TCanvas *Trial_eV_cneweffrelative=new TCanvas("Trial_eV_cneweffrelative","Trial_eV_cneweffrelative",67,55,1853,1025);
TH2D *Trial_eV_h2neweffrelative=new TH2D("Trial_eV_h2neweffrelative",Form("%s;Energy Input [eV];Relative Efficiency",common_name.c_str()),1,10.,1500.,1,-0.01,1.01);
//Trial_eV_h2neweffrelative->DrawCopy();
TLegend *Trial_eV_legneweffrelative=new TLegend(0.7,0.1,0.9,0.5);
for(int c=0;c<NtrialsDeltaChi2;c++){
	Trial_eV_geffrelative_cut_DeltaChi2[c]->SetLineColor(color_rainbow(NtrialsDeltaChi2-1-c,NtrialsDeltaChi2));
	Trial_eV_geffrelative_cut_DeltaChi2[c]->SetMarkerColor(color_rainbow(NtrialsDeltaChi2-1-c,NtrialsDeltaChi2));
	//Trial_eV_geffrelative_cut_DeltaChi2[c]->Draw("LPsame");
	//Trial_eV_legneweffrelative->AddEntry(Trial_eV_geffrelative_cut_DeltaChi2[c],Form("Delta Chi2 < %.4f",vecTrialDeltaChi2[c]),"lep");
}
//gPad->SetGridx();
//gPad->SetGridy();
//Trial_eV_legneweffrelative->DrawClone();
//Trial_eV_cneweffrelative->Update();
//Trial_eV_cneweffrelative->SaveAs("Trial_effrelative_eV_Ein.png");
//TCanvas *eV_cneweffrelative=new TCanvas("eV_cneweffrelative","eV_cneweffrelative",67,55,1853,1025);
TH2D *eV_h2neweffrelative=new TH2D("eV_h2neweffrelative",Form("%s;Energy Input [eV];Relative Efficiency",common_name.c_str()),1,10.,1500.,1,-0.01,1.01);
//eV_h2neweffrelative->DrawCopy();

eV_geffrelative_Bank_over_RT_Maintenance->SetLineColor(color_rainbow(1,6));
eV_geffrelative_Bank_over_RT_Maintenance->SetMarkerColor(color_rainbow(1,6));
//eV_geffrelative_Bank_over_RT_Maintenance->Draw("LPsame");
eV_geffrelative_Trigger_over_Bank->SetLineColor(color_rainbow(2,6));
eV_geffrelative_Trigger_over_Bank->SetMarkerColor(color_rainbow(2,6));
//eV_geffrelative_Trigger_over_Bank->Draw("LPsame");
eV_geffrelative_Chi2Norm_over_Trigger->SetLineColor(color_rainbow(3,6));
eV_geffrelative_Chi2Norm_over_Trigger->SetMarkerColor(color_rainbow(3,6));
//eV_geffrelative_Chi2Norm_over_Trigger->Draw("LPsame");
eV_geffrelative_Chi2BNorm_over_Trigger->SetLineColor(color_rainbow(5,6));
eV_geffrelative_Chi2BNorm_over_Trigger->SetMarkerColor(color_rainbow(5,6));
//eV_geffrelative_Chi2BNorm_over_Trigger->Draw("LPsame");
eV_geffrelative_DeltaChi2_over_Chi2Norm->SetLineColor(color_rainbow(4,6));
eV_geffrelative_DeltaChi2_over_Chi2Norm->SetMarkerColor(color_rainbow(4,6));
//eV_geffrelative_DeltaChi2_over_Chi2Norm->Draw("LPsame");

//gPad->SetGridx();
//gPad->SetGridy();
TLegend *eV_legneweffrelative=new TLegend(0.5,0.1,0.9,0.4);
eV_legneweffrelative->AddEntry(eV_geffrelative_Bank_over_RT_Maintenance,"Cut 10 keV Bank","lep");
eV_legneweffrelative->AddEntry(eV_geffrelative_Trigger_over_Bank,"Trigger (coinc Simu anticoinc Data)","lep");
eV_legneweffrelative->AddEntry(eV_geffrelative_Chi2Norm_over_Trigger,"Cut Chi2 Normal","lep");
eV_legneweffrelative->AddEntry(eV_geffrelative_Chi2BNorm_over_Trigger,"Cut Chi2 B Normal","lep");
//eV_legneweffrelative->AddEntry(eV_geffrelative_DeltaChi2_over_Chi2Norm,Form("Cut Delta Chi2 < %.3f" ,par0_DELTA),"lep");
//eV_legneweffrelative->DrawClone();
//eV_cneweffrelative->Update();
//eV_cneweffrelative->SaveAs("effrelative_eV_Ein.png");

//TCanvas *keV_cneweffrelative=new TCanvas("keV_cneweffrelative","keV_cneweffrelative",67,55,1853,1025);
TH2D *keV_h2neweffrelative=new TH2D("keV_h2neweffrelative",Form("%s;Energy Input [keV];Relative Efficiency",common_name.c_str()),1,1./1000.,11.,1,-0.01,1.01);
//keV_h2neweffrelative->DrawCopy();

/*
keV_geffrelative_Bank_over_RT_Maintenance->SetLineColor(color_rainbow(1,6));
keV_geffrelative_Bank_over_RT_Maintenance->SetMarkerColor(color_rainbow(1,6));
keV_geffrelative_Bank_over_RT_Maintenance->Draw("LPsame");
keV_geffrelative_Trigger_over_Bank->SetLineColor(color_rainbow(2,6));
keV_geffrelative_Trigger_over_Bank->SetMarkerColor(color_rainbow(2,6));
keV_geffrelative_Trigger_over_Bank->Draw("LPsame");
keV_geffrelative_Chi2Norm_over_Trigger->SetLineColor(color_rainbow(3,6));
keV_geffrelative_Chi2Norm_over_Trigger->SetMarkerColor(color_rainbow(3,6));
keV_geffrelative_Chi2Norm_over_Trigger->Draw("LPsame");
keV_geffrelative_Chi2BNorm_over_Trigger->SetLineColor(color_rainbow(5,6));
keV_geffrelative_Chi2BNorm_over_Trigger->SetMarkerColor(color_rainbow(5,6));
keV_geffrelative_Chi2BNorm_over_Trigger->Draw("LPsame");
keV_geffrelative_DeltaChi2_over_Chi2Norm->SetLineColor(color_rainbow(4,6));
keV_geffrelative_DeltaChi2_over_Chi2Norm->SetMarkerColor(color_rainbow(4,6));
//keV_geffrelative_DeltaChi2_over_Chi2Norm->Draw("LPsame");
*/
//gPad->SetLogx();
//gPad->SetGridx();
// gPad->SetGridy();
TLegend *keV_legneweffrelative=new TLegend(0.5,0.1,0.9,0.4);
keV_legneweffrelative->AddEntry(keV_geffrelative_Bank_over_RT_Maintenance,"Cut 10 keV Bank","lep");
keV_legneweffrelative->AddEntry(keV_geffrelative_Trigger_over_Bank,"Trigger (coinc Simu anticoinc Data)","lep");
keV_legneweffrelative->AddEntry(keV_geffrelative_Chi2Norm_over_Trigger,"Cut Chi2 Normal","lep");
keV_legneweffrelative->AddEntry(keV_geffrelative_Chi2BNorm_over_Trigger,"Cut Chi2 B Normal","lep");
//keV_legneweffrelative->AddEntry(keV_geffrelative_DeltaChi2_over_Chi2Norm,Form("Cut Delta Chi2 < %.3f" ,par0_DELTA),"lep");
keV_legneweffrelative->DrawClone();

/////////////////////////////////////////////////////////////
TCanvas *keV_cneweff=new TCanvas("keV_cneweff","keV_cneweff",67,55,1853,1025);
TH2D *keV_h2neweff=new TH2D("keV_h2neweff",Form("%s;Energy Input [keV];Efficiency",common_name.c_str()),1,1./1000.,11.,1,-0.01,1.01);
keV_h2neweff->DrawCopy();

keV_geff_cut_RT_Maintenance->SetLineColor(color_rainbow(0,6));
keV_geff_cut_RT_Maintenance->SetMarkerColor(color_rainbow(0,6));
keV_geff_cut_RT_Maintenance->Draw("LPsame");
keV_geff_cut_Bank->SetLineColor(color_rainbow(1,6));
keV_geff_cut_Bank->SetMarkerColor(color_rainbow(1,6));
keV_geff_cut_Bank->Draw("LPsame");
keV_geff_Trigger->SetLineColor(color_rainbow(2,6));
keV_geff_Trigger->SetMarkerColor(color_rainbow(2,6));
keV_geff_Trigger->Draw("LPsame");
keV_geff_cut_Chi2NormION->SetLineColor(color_rainbow(5,6));
keV_geff_cut_Chi2NormION->SetMarkerColor(color_rainbow(5,6));
keV_geff_cut_Chi2NormION->Draw("LPsame");
keV_geff_cut_Eion->SetLineColor(color_rainbow(1,6));
keV_geff_cut_Eion->SetMarkerColor(color_rainbow(1,6));
keV_geff_cut_Eion->Draw("LPsame");
keV_geff_cut_chalfid->SetLineColor(color_rainbow(0,6));
keV_geff_cut_chalfid->SetMarkerColor(color_rainbow(0,6));
keV_geff_cut_chalfid->Draw("LPsame");
keV_geff_cut_Chi2Norm->SetLineColor(color_rainbow(3,6));
keV_geff_cut_Chi2Norm->SetMarkerColor(color_rainbow(3,6));
keV_geff_cut_Chi2Norm->Draw("LPsame");
keV_geff_cut_Chi2BNorm->SetLineColor(color_rainbow(5,6));
keV_geff_cut_Chi2BNorm->SetMarkerColor(color_rainbow(5,6));
keV_geff_cut_Chi2BNorm->Draw("LPsame");
keV_geff_cut_DeltaChi2->SetLineColor(color_rainbow(4,6));
keV_geff_cut_DeltaChi2->SetMarkerColor(color_rainbow(4,6));
//keV_geff_cut_DeltaChi2->Draw("LPsame");
gPad->SetLogx();
gPad->SetGridx();
gPad->SetGridy();
TLegend *keV_legneweff=new TLegend(0.5,0.1,0.9,0.4);
keV_legneweff->AddEntry(keV_geff_cut_RT_Maintenance,LegendcutRT_Maintenance.c_str(),"lep");
keV_legneweff->AddEntry(keV_geff_cut_Bank,"Cut 10 keV Bank","lep");
keV_legneweff->AddEntry(keV_geff_Trigger,"Trigger (coinc Simu anticoinc Data)","lep");
keV_legneweff->AddEntry(keV_geff_cut_Chi2NormION,"Cut Chi2 ION","lep");
keV_legneweff->AddEntry(keV_geff_cut_Eion,"Cut Eion >0","lep");
keV_legneweff->AddEntry(keV_geff_cut_chalfid,"Fid heat cut");
keV_legneweff->AddEntry(keV_geff_cut_Chi2Norm,"Cut Chi2 Normal","lep");
keV_legneweff->AddEntry(keV_geff_cut_Chi2BNorm,"Cut Chi2 B Normal","lep");
//keV_legneweff->AddEntry(keV_geff_cut_DeltaChi2,Form("Cut Delta Chi2 < %.3f" ,par0_DELTA),"lep");
//keV_legneweff->DrawClone();

TCanvas *eV_cneweff=new TCanvas("eV_cneweff","eV_cneweff",67,55,1853,1025);
TH2D *eV_h2neweff=new TH2D("eV_h2neweff",Form("%s;Energy Input [eV];Efficiency",common_name.c_str()),1,1.,1500.,1,-0.01,1.01);
eV_h2neweff->DrawCopy();

eV_geff_cut_RT_Maintenance->SetLineColor(color_rainbow(0,6));
eV_geff_cut_RT_Maintenance->SetMarkerColor(color_rainbow(0,6));
eV_geff_cut_RT_Maintenance->Draw("LPsame");
eV_geff_cut_Bank->SetLineColor(color_rainbow(1,6));
eV_geff_cut_Bank->SetMarkerColor(color_rainbow(1,6));
eV_geff_cut_Bank->Draw("LPsame");
eV_geff_Trigger->SetLineColor(color_rainbow(2,6));
eV_geff_Trigger->SetMarkerColor(color_rainbow(2,6));
eV_geff_Trigger->Draw("LPsame");
eV_geff_cut_Chi2NormION->SetLineColor(color_rainbow(5,6));
eV_geff_cut_Chi2NormION->SetMarkerColor(color_rainbow(5,6));
eV_geff_cut_Chi2NormION->Draw("LPsame");
eV_geff_cut_Eion->SetLineColor(color_rainbow(1,6));
eV_geff_cut_Eion->SetMarkerColor(color_rainbow(1,6));
eV_geff_cut_Eion->Draw("LPsame");
eV_geff_cut_chalfid->SetLineColor(color_rainbow(0,6));
eV_geff_cut_chalfid->SetMarkerColor(color_rainbow(0,6));
eV_geff_cut_chalfid->Draw("LPsame");
eV_geff_cut_Chi2Norm->SetLineColor(color_rainbow(3,6));
eV_geff_cut_Chi2Norm->SetMarkerColor(color_rainbow(3,6));
eV_geff_cut_Chi2Norm->Draw("LPsame");
eV_geff_cut_Chi2BNorm->SetLineColor(color_rainbow(5,6));
eV_geff_cut_Chi2BNorm->SetMarkerColor(color_rainbow(5,6));
eV_geff_cut_Chi2BNorm->Draw("LPsame");
eV_geff_cut_DeltaChi2->SetLineColor(color_rainbow(4,6));
eV_geff_cut_DeltaChi2->SetMarkerColor(color_rainbow(4,6));
//eV_geff_cut_DeltaChi2->Draw("LPsame");
TFile* output_efficiency = new TFile(Form("Efficiencies_%3.2f_%s_cutFid.root",cutIon*1000,common_name.c_str()),"RECREATE");
eV_geff_cut_Chi2BNorm->Write(Form("Efficienciy_cutChihB_%s",common_name.c_str()));
eV_geff_cut_Chi2Norm ->Write(Form("Efficienciy_cutChihA_%s",common_name.c_str()));
output_efficiency->Close();


gPad->SetGridx();
gPad->SetGridy();
TLegend *legneweff=new TLegend(0.5,0.1,0.9,0.4);
legneweff->AddEntry(eV_geff_cut_RT_Maintenance,LegendcutRT_Maintenance.c_str(),"lep");
legneweff->AddEntry(eV_geff_cut_Bank,"Cut 10 keV Bank","lep");
legneweff->AddEntry(eV_geff_Trigger,"Trigger (coinc Simu anticoinc Data)","lep");
legneweff->AddEntry(eV_geff_cut_Chi2NormION,"Cut Chi2 ION","lep");
legneweff->AddEntry(eV_geff_cut_Eion,Form("Cut Eion < %3.0f",cutIon*1000),"lep");
legneweff->AddEntry(eV_geff_cut_chalfid,"Fid heat cut");
legneweff->AddEntry(eV_geff_cut_Chi2Norm,"Cut Chi2 Normal","lep");
legneweff->AddEntry(eV_geff_cut_Chi2BNorm,"Cut Chi2 B Normal","lep");
//legneweff->AddEntry(eV_geff_cut_DeltaChi2,Form("Cut Delta Chi2 < %.3f" ,par0_DELTA),"lep");
legneweff->DrawClone();
if(verbose > 0) cout<<"end calculate"<<endl;
eV_cneweff->Update();
eV_cneweff->SaveAs(Form("effcuts_eV_%s_%3.2feVion_cutFid.png",common_name.c_str(),cutIon*1000));
//toto->Run();
///////////////////////////////////// PLOT ///////////////////////////
TCanvas *c1=new TCanvas("c1","c1",67,55,1853,1025);
c1->Divide(2,2);
TLegend *leg=new TLegend(0.1,0.1,0.9,0.9);

for(int tol_index=0;tol_index<N_tolerances;tol_index++){	
	string option;
	if(tol_index==0){option="E";}else{option="Esame";}
	int color=color_rainbow(tol_index,N_tolerances);	
	c1->cd(1);
	he_simu_coinc_data[tol_index]->SetTitle("(Simu) coinc Data;Energy [keV];Efficiency");	
	he_simu_coinc_data[tol_index]->SetLineColor(color);
	he_simu_coinc_data[tol_index]->DrawCopy(option.c_str());
	gPad->SetLogx();
	c1->cd(2);
	he_simu_coinc_datasimu[tol_index]->SetTitle("(Simu) coinc (Data+Simu);Energy [keV];Efficiency");	
	he_simu_coinc_datasimu[tol_index]->SetLineColor(color);
	he_simu_coinc_datasimu[tol_index]->DrawCopy(option.c_str());
	gPad->SetLogx();
	c1->cd(3);
	he_simu_coinc_datasimu_anticoinc_data[tol_index]->SetTitle("(Simu) coinc (Data+Simu) anticoinc Data;Energy [keV];Efficiency");	
	he_simu_coinc_datasimu_anticoinc_data[tol_index]->SetLineColor(color);
	he_simu_coinc_datasimu_anticoinc_data[tol_index]->DrawCopy(option.c_str());
	gPad->SetLogx();
	//c1->cd(4);
	leg->AddEntry(he_simu_coinc_data[tol_index],Form("tol = %1.1f timestamp",vectolerances[tol_index]),"l");

}

c1->cd(4);
leg->DrawClone();
c1->Update();
//c1->SaveAs(Form("%s/c1-%s_%s.eps",Figurespath.c_str(),study_index.c_str(),prefixrootfilesaveLines.c_str()));
//c1->SaveAs(Form("%s/c1-%s_%s.png",Figurespath.c_str(),study_index.c_str(),prefixrootfilesaveLines.c_str()));
if(verbose > 0) cout<<"canvas c2"<<endl;
TCanvas *c2=new TCanvas("c2","c2",67,55,1853,1025);
c2->Divide(2,2);
c2->cd(1);
TLegend *leg2=new TLegend(0.1,0.7,0.5,0.9);
/*he_simu_coinc_datasimu_anticoinc_data[0]->SetTitle(Form("(Simu) coinc (Data+Simu) anticoinc Data tol = %1.1f ts;Energy [keV];Efficiency",vectolerances[0])); 
he_simu_coinc_datasimu_anticoinc_data_cut1[0]->SetTitle(Form("(Simu) coinc (Data+Simu) anticoinc Data with Cut tol = %1.1f ts;Energy [keV];Efficiency",vectolerances[0])); 
CorrectedE_he_simu_coinc_datasimu_anticoinc_data[0]->SetTitle(Form("(Simu) coinc (Data+Simu) anticoinc Data tol = %1.1f ts;Energy [keV];Efficiency",vectolerances[0])); 
CorrectedE_h2eff_simu_coinc_datasimu_anticoinc_data_cut1[0]->SetTitle(Form("(Simu) coinc (Data+Simu) anticoinc Data tol = %1.1f ts;Energy [keV];Efficiency",vectolerances[0])); 
CorrectedE_he_simu_coinc_datasimu_anticoinc_data[0]->SetLineColor(color_rainbow(0,4));
CorrectedE_h2eff_simu_coinc_datasimu_anticoinc_data_cut1[0]->SetLineColor(color_rainbow(3,4));
he_simu_coinc_datasimu_anticoinc_data[0]->SetLineColor(color_rainbow(1,4));
he_simu_coinc_datasimu_anticoinc_data_cut1[0]->SetLineColor(color_rainbow(2,4));
leg2->AddEntry(CorrectedE_he_simu_coinc_datasimu_anticoinc_data[0],"no cut Energy_reconstructed","l");
leg2->AddEntry(CorrectedE_h2eff_simu_coinc_datasimu_anticoinc_data_cut1[0],"all cuts Energy_reconstructed","l");
leg2->AddEntry(he_simu_coinc_datasimu_anticoinc_data[0],"no cut Energy_In","l");
leg2->AddEntry(he_simu_coinc_datasimu_anticoinc_data_cut1[0],"all cuts Energy_In","l");
he_simu_coinc_datasimu_anticoinc_data[0]->DrawCopy("E");
he_simu_coinc_datasimu_anticoinc_data_cut1[0]->DrawCopy("Esame");
CorrectedE_he_simu_coinc_datasimu_anticoinc_data[0]->DrawCopy("Esame");
CorrectedE_h2eff_simu_coinc_datasimu_anticoinc_data_cut1[0]->DrawCopy("Esame");
leg2->DrawClone();*/
int markerstyle=20;
double markersize=0.3;
h2D_0_NORMALchi2BvsE[0]->SetTitle(Form("Heat B (Simu) coinc (Data+Simu) anticoinc Data tol = %1.1f ts;Energy reconstructed [keV];Chi2 Normal",vectolerances[0])); 
h2D_0_NORMALchi2BvsE[0]->SetMarkerStyle(markerstyle);
h2D_0_NORMALchi2BvsE[0]->SetMarkerSize(markersize);
h2D_0_NORMALchi2BvsE[0]->SetMarkerColor(color_rainbow(0,3));
h2D_0_NORMALchi2BvsE[0]->DrawCopy("scat");
leg2->AddEntry(h2D_0_NORMALchi2BvsE[0],"not in coincidence","p");

h2D_1_NORMALchi2BvsE[0]->SetMarkerStyle(markerstyle);
h2D_1_NORMALchi2BvsE[0]->SetMarkerSize(markersize);
h2D_1_NORMALchi2BvsE[0]->SetMarkerColor(color_rainbow(1,3));
h2D_1_NORMALchi2BvsE[0]->DrawCopy("scatsame");
leg2->AddEntry(h2D_1_NORMALchi2BvsE[0],"coincidence but fail cuts","p");

h2D_2_NORMALchi2BvsE[0]->SetMarkerStyle(markerstyle);
h2D_2_NORMALchi2BvsE[0]->SetMarkerSize(markersize);
h2D_2_NORMALchi2BvsE[0]->SetMarkerColor(color_rainbow(2,3)); 
h2D_2_NORMALchi2BvsE[0]->DrawCopy("scatsame");
leg2->AddEntry(h2D_2_NORMALchi2BvsE[0],"coincidence and pass cuts","p");
fchi2B_Normal->DrawCopy("same");
leg2->DrawClone();

if(optionlogmonitoring){
gPad->SetLogx();
}
gPad->SetGridx();
gPad->SetGridy();
if(saveLines){
		foutsaveLines->cd(); 
		he_simu_coinc_datasimu_anticoinc_data[0]->Write("h_eff_no_cut");
		he_simu_coinc_datasimu_anticoinc_data_cut1[0]->Write("h_eff_cut");
		keV_geff_cut_DeltaChi2->Write("g_eff_Ein_finalcut");
}

TLegend *leg3=new TLegend(0.1,0.7,0.5,0.9);
c2->cd(2);
h2D_0_NORMALchi2vsE[0]->SetTitle(Form("Heat A (Simu) coinc (Data+Simu) anticoinc Data tol = %1.1f ts;Energy reconstructed [keV];Chi2 Normal",vectolerances[0])); 
h2D_0_NORMALchi2vsE[0]->SetMarkerStyle(markerstyle);
h2D_0_NORMALchi2vsE[0]->SetMarkerSize(markersize);
h2D_0_NORMALchi2vsE[0]->SetMarkerColor(color_rainbow(0,3));
h2D_0_NORMALchi2vsE[0]->DrawCopy("scat");
leg3->AddEntry(h2D_0_NORMALchi2vsE[0],"not in coincidence","p");

h2D_1_NORMALchi2vsE[0]->SetMarkerStyle(markerstyle);
h2D_1_NORMALchi2vsE[0]->SetMarkerSize(markersize);
h2D_1_NORMALchi2vsE[0]->SetMarkerColor(color_rainbow(1,3));
h2D_1_NORMALchi2vsE[0]->DrawCopy("scatsame");
leg3->AddEntry(h2D_1_NORMALchi2vsE[0],"coincidence but fail cuts","p");

h2D_2_NORMALchi2vsE[0]->SetMarkerStyle(markerstyle);
h2D_2_NORMALchi2vsE[0]->SetMarkerSize(markersize);
h2D_2_NORMALchi2vsE[0]->SetMarkerColor(color_rainbow(2,3)); 
h2D_2_NORMALchi2vsE[0]->DrawCopy("scatsame");
leg3->AddEntry(h2D_2_NORMALchi2vsE[0],"coincidence and pass cuts","p");
leg3->DrawClone();

fchi2_Normal->DrawCopy("same");

//gPad->SetLogx();
if(optionlogmonitoring){
gPad->SetLogy();
gPad->SetLogx();
}
c2->cd(3);
TLegend *leg3bis=new TLegend(0.1,0.7,0.5,0.9);
h2D_0_DELTAchi2vsE[0]->SetTitle(Form("(Simu) coinc (Data+Simu) anticoinc Data tol = %1.1f ts;Energy reconstructed [keV]; ION Chi2",vectolerances[0])); 
h2D_0_DELTAchi2vsE[0]->SetMarkerStyle(markerstyle);
h2D_0_DELTAchi2vsE[0]->SetMarkerSize(markersize);
h2D_0_DELTAchi2vsE[0]->SetMarkerColor(color_rainbow(0,3));
h2D_0_DELTAchi2vsE[0]->DrawCopy("scat");
leg3bis->AddEntry(h2D_0_DELTAchi2vsE[0],"not in coincidence","p");

h2D_1_DELTAchi2vsE[0]->SetMarkerStyle(markerstyle);
h2D_1_DELTAchi2vsE[0]->SetMarkerSize(markersize);
h2D_1_DELTAchi2vsE[0]->SetMarkerColor(color_rainbow(1,3));
h2D_1_DELTAchi2vsE[0]->DrawCopy("scatsame");
leg3bis->AddEntry(h2D_1_DELTAchi2vsE[0],"coincidence but fail cuts","p");

h2D_2_DELTAchi2vsE[0]->SetMarkerStyle(markerstyle);
h2D_2_DELTAchi2vsE[0]->SetMarkerSize(markersize);
h2D_2_DELTAchi2vsE[0]->SetMarkerColor(color_rainbow(2,3)); 
h2D_2_DELTAchi2vsE[0]->DrawCopy("scatsame");
leg3bis->AddEntry(h2D_2_DELTAchi2vsE[0],"coincidence and pass cuts","p");
leg3bis->DrawClone();
fchi2_ION->DrawCopy("same");
if(optionlogmonitoring){
   gPad->SetLogx();
}
TLegend *leg4=new TLegend(0.1,0.7,0.5,0.9);
c2->cd(4);
h2D_0_EoutvsE[0]->SetTitle(Form("(Simu) coinc (Data+Simu) anticoinc Data tol = %1.1f ts;Energy_In [keV];Energy reconstructed [keV]",vectolerances[0])); 
h2D_0_EoutvsE[0]->SetMarkerStyle(markerstyle);
h2D_0_EoutvsE[0]->SetMarkerSize(markersize);
h2D_0_EoutvsE[0]->SetMarkerColor(color_rainbow(0,3));
h2D_0_EoutvsE[0]->DrawCopy("scat");
leg4->AddEntry(h2D_0_EoutvsE[0],"not in coincidence","p");
h2D_1_EoutvsE[0]->SetMarkerStyle(markerstyle);
h2D_1_EoutvsE[0]->SetMarkerSize(markersize);
h2D_1_EoutvsE[0]->SetMarkerColor(color_rainbow(1,3));
h2D_1_EoutvsE[0]->DrawCopy("scatsame");
leg4->AddEntry(h2D_1_EoutvsE[0],"coincidence but fail cuts","p");
h2D_2_EoutvsE[0]->SetMarkerStyle(markerstyle);
h2D_2_EoutvsE[0]->SetMarkerSize(markersize);
h2D_2_EoutvsE[0]->SetMarkerColor(color_rainbow(2,3)); 
h2D_2_EoutvsE[0]->DrawCopy("scatsame");
leg4->AddEntry(h2D_2_EoutvsE[0],"coincidence and pass cuts","p");
TF1 *fx=new TF1("fx","[0]*x",0.0001,10);
fx->SetParameter(0,conversionkeV2ADU);
fx->SetNpx(10000);
fx->SetLineColor(kGray+2);
fx->DrawCopy("same");
leg4->DrawClone();
if(optionlogmonitoring){
gPad->SetLogx();
gPad->SetLogy();
}
c2->SaveAs("chi2_vs_energy_both_channel.png");
if(verbose > 0) cout<<"canvas cst"<<endl;
c2->Update();
TCanvas *cst=new TCanvas("cst","cst",67,55,1853,1025);
cst->Divide(2,2);
THStack *hs_0 = new THStack("hs_0",Form("All ;Energy_Reconstructed [keV];Counts"));
THStack *hs_1 = new THStack("hs_1",Form("Coinc anticoinc Data tol = %1.1f ts;Energy_Reconstructed [keV] ;Counts",vectolerances[0]));
THStack *hs_2 = new THStack("hs_2",Form("Coinc anticoinc Data tol = %1.1f ts and Pass Cuts ;Energy_Reconstructed [keV] ;Counts",vectolerances[0]));
THStack *hs_N = new THStack("hs_N",Form("Coinc anticoinc Data tol = %1.1f ts and Pass Cuts Normalized;Energy_Reconstructed [keV];Counts / bin",vectolerances[0]));
THStack *hs_0B = new THStack("hs_0B",Form("All ;Energy_Reconstructed [keV];Counts"));
THStack *hs_1B = new THStack("hs_1B",Form("Coinc anticoinc Data tol = %1.1f ts;Energy_Reconstructed [keV] ;Counts",vectolerances[0]));
THStack *hs_2B = new THStack("hs_2B",Form("Coinc anticoinc Data tol = %1.1f ts and Pass Cuts ;Energy_Reconstructed [keV] ;Counts",vectolerances[0]));
THStack *hs_NB = new THStack("hs_NB",Form("Coinc anticoinc Data tol = %1.1f ts and Pass Cuts Normalized;Energy_Reconstructed [keV];Counts / bin",vectolerances[0]));
TLegend *leghs=new TLegend(0.1,0.1,0.9,0.9);
if(saveLines){ 
   foutsaveLines->cd();
}
if(verbose > 0) cout<<Nlines<<endl;

TFile* output_chi2 = new TFile(Form("output_chi2_%s_%3.2f_%s.root",Runs_to_analyse.c_str(),cutIon*1000,common_name.c_str()),"RECREATE");
output_chi2->cd();

for(int l=0;l<Nlines;l++){

   for(int itchi =  0 ; itchi < 6 ; itchi ++){
      continue;
      h1Chi2Line_2_EoutvsChi2[0][l][itchi]->Fit("gaus","","",0.1,5.);
      h1Chi2Line_2_EoutvsChi2[0][l][itchi]->Write();
      h1Chi2Line_0_EoutvsChi2[0][l][itchi]->Fit("gaus","","",0.1,5.);
      h1Chi2Line_0_EoutvsChi2[0][l][itchi]->Write();
      h1Chi2Line_1_EoutvsChi2[0][l][itchi]->Fit("gaus","","",0.1,5.);
      h1Chi2Line_1_EoutvsChi2[0][l][itchi]->Write();
                  
      chi2_vs_Ein_2[itchi]->SetPoint(l,Eline[l],h1Chi2Line_2_EoutvsChi2[0][l][itchi]->GetFunction("gaus")->GetParameter(1));
      chi2_vs_Ein_1[itchi]->SetPoint(l,Eline[l],h1Chi2Line_1_EoutvsChi2[0][l][itchi]->GetFunction("gaus")->GetParameter(1));
      chi2_vs_Ein_0[itchi]->SetPoint(l,Eline[l],h1Chi2Line_0_EoutvsChi2[0][l][itchi]->GetFunction("gaus")->GetParameter(1));
      sigchi2_vs_Ein_2[itchi]->SetPoint(l,Eline[l],h1Chi2Line_2_EoutvsChi2[0][l][itchi]->GetFunction("gaus")->GetParameter(2));
      sigchi2_vs_Ein_1[itchi]->SetPoint(l,Eline[l],h1Chi2Line_1_EoutvsChi2[0][l][itchi]->GetFunction("gaus")->GetParameter(2));
   }
   bool fit = true;
   if(fit == true){
   h1EoutLine_0_EoutvsE[0][l]->Fit("gaus","","",Eline[l]-0.1,Eline[l]+0.1);
   h1EoutLine_1_EoutvsE[0][l]->Fit("gaus","","",Eline[l]-0.1,Eline[l]+0.1);
   h1EoutLine_2_EoutvsE[0][l]->Fit("gaus","","",Eline[l]-0.1,Eline[l]+0.1);
   
   
   h1EoutLine_0_EoutvsEB[0][l]->Fit("gaus","","",Eline[l]-0.1,Eline[l]+0.1);
   h1EoutLine_1_EoutvsEB[0][l]->Fit("gaus","","",Eline[l]-0.1,Eline[l]+0.1);
   h1EoutLine_2_EoutvsEB[0][l]->Fit("gaus","","",Eline[l]-0.1,Eline[l]+0.1);
   
   h1Line_0_EoutvsEion[0][l][0]->Fit("gaus","","",Eline[l]-0.2,Eline[l]+0.2);
   h1Line_0_EoutvsEion[0][l][1]->Fit("gaus","","",Eline[l]-0.2,Eline[l]+0.2);
   
   h1Line_1_EoutvsEion[0][l][0]->Fit("gaus","","",Eline[l]-0.2,Eline[l]+0.2);
   h1Line_1_EoutvsEion[0][l][0]->Write();
   h1Line_1_EoutvsEion[0][l][1]->Fit("gaus","","",Eline[l]-0.2,Eline[l]+0.2);
   h1Line_1_EoutvsEion[0][l][1]->Write();
   h1Line_2_EoutvsEion[0][l][0]->Fit("gaus","","",Eline[l]-0.2,Eline[l]+0.2);
   h1Line_2_EoutvsEion[0][l][1]->Fit("gaus","","",Eline[l]-0.2,Eline[l]+0.2);
   h1Line_2_EoutvsEion[0][l][0]->Write();
   h1Line_2_EoutvsEion[0][l][1]->Write();
    Ea_vs_Ein_0->SetPoint(l,Eline[l],h1EoutLine_0_EoutvsE[0][l]->GetFunction("gaus")->GetParameter(1));
	Eb_vs_Ein_0->SetPoint(l,Eline[l],h1EoutLine_0_EoutvsEB[0][l]->GetFunction("gaus")->GetParameter(1)); 
	Ea_vs_Ein_1->SetPoint(l,Eline[l],h1EoutLine_1_EoutvsE[0][l]->GetFunction("gaus")->GetParameter(1));
	Eb_vs_Ein_1->SetPoint(l,Eline[l],h1EoutLine_1_EoutvsEB[0][l]->GetFunction("gaus")->GetParameter(1));
	Ea_vs_Ein_2->SetPoint(l,Eline[l],h1EoutLine_2_EoutvsE[0][l]->GetFunction("gaus")->GetParameter(1));
	Eb_vs_Ein_2->SetPoint(l,Eline[l],h1EoutLine_2_EoutvsEB[0][l]->GetFunction("gaus")->GetParameter(1));
	
	sigEa_vs_Ein_1->SetPoint(l,Eline[l],h1EoutLine_1_EoutvsE[0][l]->GetFunction("gaus")->GetParameter(2));
	sigEb_vs_Ein_1->SetPoint(l,Eline[l],h1EoutLine_1_EoutvsEB[0][l]->GetFunction("gaus")->GetParameter(2));
	sigEa_vs_Ein_2->SetPoint(l,Eline[l],h1EoutLine_2_EoutvsE[0][l]->GetFunction("gaus")->GetParameter(2));
	sigEb_vs_Ein_2->SetPoint(l,Eline[l],h1EoutLine_2_EoutvsEB[0][l]->GetFunction("gaus")->GetParameter(2));
	
	Et_vs_EionB_2->SetPoint(l,h1EoutLine_2_EoutvsE[0][l]->GetFunction("gaus")->GetParameter(1),h1Line_2_EoutvsEion[0][l][0]->GetFunction("gaus")->GetParameter(1));
	Et_vs_EionB_1->SetPoint(l,h1EoutLine_1_EoutvsE[0][l]->GetFunction("gaus")->GetParameter(1),h1Line_1_EoutvsEion[0][l][0]->GetFunction("gaus")->GetParameter(1));
	
	Et_vs_EionFid_2->SetPoint(l,h1EoutLine_2_EoutvsE[0][l]->GetFunction("gaus")->GetParameter(1),h1Line_2_EoutvsEion[0][l][0]->GetFunction("gaus")->GetParameter(1));
	Et_vs_EionFid_1->SetPoint(l,h1EoutLine_1_EoutvsE[0][l]->GetFunction("gaus")->GetParameter(1),h1Line_1_EoutvsEion[0][l][0]->GetFunction("gaus")->GetParameter(1));
	
	}
	
	int color=color_rainbow(l,Nlines);
	h1EoutLine_0_EoutvsE[0][l]->SetLineColor(color);
	h1EoutLine_0_EoutvsE[0][l]->SetFillColor(color);
	h1EoutLine_0_EoutvsE[0][l]->SetLineWidth(2);
	hs_0->Add(h1EoutLine_0_EoutvsE[0][l],"HIST");
	h1EoutLine_0_EoutvsEB[0][l]->SetLineColor(color);
	h1EoutLine_0_EoutvsEB[0][l]->SetFillColor(color);
	h1EoutLine_0_EoutvsEB[0][l]->SetLineWidth(2);
	hs_0B->Add(h1EoutLine_0_EoutvsEB[0][l],"HIST");
	int Nall=h1EoutLine_0_EoutvsE[0][l]->GetEntries();
	
	h1EoutLine_1_EoutvsE[0][l]->SetLineColor(color);
	h1EoutLine_1_EoutvsE[0][l]->SetFillColor(color);
	h1EoutLine_1_EoutvsE[0][l]->SetLineWidth(2);
	hs_1->Add(h1EoutLine_1_EoutvsE[0][l],"HIST");
	h1EoutLine_1_EoutvsEB[0][l]->SetLineColor(color);
	h1EoutLine_1_EoutvsEB[0][l]->SetFillColor(color);
	h1EoutLine_1_EoutvsEB[0][l]->SetLineWidth(2);
	hs_1B->Add(h1EoutLine_1_EoutvsEB[0][l],"HIST");
	int Ncoinc=h1EoutLine_1_EoutvsE[0][l]->GetEntries();
	double eff_coinc;
	if(Nall>0){
	   eff_coinc=(double)Ncoinc/Nall;
	}else{
	   eff_coinc=0.;
	} // ajout minuit
	h1EoutLine_2_EoutvsE[0][l]->SetLineColor(color);
	h1EoutLine_2_EoutvsE[0][l]->SetFillColor(color);
	h1EoutLine_2_EoutvsE[0][l]->SetLineWidth(2);
	hs_2->Add(h1EoutLine_2_EoutvsE[0][l],"HIST");
	h1EoutLine_2_EoutvsEB[0][l]->SetLineColor(color);
	h1EoutLine_2_EoutvsEB[0][l]->SetFillColor(color);
	h1EoutLine_2_EoutvsEB[0][l]->SetLineWidth(2);
	hs_2B->Add(h1EoutLine_2_EoutvsEB[0][l],"HIST");
	
	int Ncoinc_cut=h1EoutLine_2_EoutvsE[0][l]->GetEntries();
	double eff_coinc_cut;
	if(Nall>0){
	   eff_coinc_cut=(double)Ncoinc_cut/Nall; 
	}else{
	   eff_coinc_cut=0.; 
	} // ajout minuit
	h1EoutLine_2Normalized_EoutvsE[0][l]->SetLineColor(color);
	h1EoutLine_2Normalized_EoutvsE[0][l]->SetLineWidth(2);
	h1EoutLine_2Normalized_EoutvsEB[0][l]->SetLineColor(color);
	h1EoutLine_2Normalized_EoutvsEB[0][l]->SetLineWidth(2);
	if(Nall>0){ // ajout minuit
		h1EoutLine_2Normalized_EoutvsE[0][l]->Scale(1./Nall);
		h1EoutLine_2Normalized_EoutvsEB[0][l]->Scale(1./Nall);
	}
	hs_N->Add(h1EoutLine_2Normalized_EoutvsE[0][l],"HIST");
	hs_NB->Add(h1EoutLine_2Normalized_EoutvsEB[0][l],"HIST");
	double IntNorm=h1EoutLine_2Normalized_EoutvsE[0][l]->Integral();	
	if(SIMU_Type=="Line"){
		h1EoutLine_0_EoutvsE[0][l]->SetTitle(Form("%s",vec_Lines[l].c_str()));
		h1EoutLine_2Normalized_EoutvsE[0][l]->SetTitle(Form("%s",vec_Lines[l].c_str()));
		leghs->AddEntry(h1EoutLine_2_EoutvsE[0][l],Form("%s keV, Nall=%i, Ncoinc=%i (#epsilon=%.2f), Ncoinc_cut=%i (#epsilon=%.2f) IntNorm=%.3f",vec_Lines[l].c_str(),Nall,Ncoinc,eff_coinc,Ncoinc_cut,eff_coinc_cut,IntNorm),"l");  // they are all same color so give any of the 3 histo
		h1EoutLine_0_EoutvsEB[0][l]->SetTitle(Form("%s",vec_Lines[l].c_str()));
		h1EoutLine_2Normalized_EoutvsEB[0][l]->SetTitle(Form("%s",vec_Lines[l].c_str()));
		leghs->AddEntry(h1EoutLine_2_EoutvsEB[0][l],Form("%s keV, Nall=%i, Ncoinc=%i (#epsilon=%.2f), Ncoinc_cut=%i (#epsilon=%.2f) IntNorm=%.3f",vec_Lines[l].c_str(),Nall,Ncoinc,eff_coinc,Ncoinc_cut,eff_coinc_cut,IntNorm),"l");
	}	
	if(saveLines){
		h1EoutLine_0_EoutvsE[0][l]->Write(Form("h_All_Line_%i",l));
		h1EoutLine_2Normalized_EoutvsE[0][l]->Write(Form("h_Line_%i",l));
		h1EoutLine_0_EoutvsEB[0][l]->Write(Form("h_AllB_Line_%i",l));
		h1EoutLine_2Normalized_EoutvsEB[0][l]->Write(Form("h_LineB_%i",l));
	}		
}

for(int itchi = 0 ; itchi< 6 ; itchi++){
      break;
      chi2_vs_Ein_2[itchi]->Write(Form("chi_passcut_coinc_%i",itchi));
      chi2_vs_Ein_1[itchi]->Write(Form("chi_coinc_%i",itchi));
      chi2_vs_Ein_0[itchi]->Write(Form("chi_notcoinc_%i",itchi));
      sigchi2_vs_Ein_2[itchi]->Write(Form("sigma_chi_passcut_coinc_%i",itchi));
      sigchi2_vs_Ein_1[itchi]->Write(Form("sigma_chi_coinc_%i",itchi));
}
Ea_vs_Ein_0->Write("EAout_vsEin_notcoinc");
Eb_vs_Ein_0->Write("EBout_vsEin_notcoinc");
Ea_vs_Ein_1->Write("EAout_vsEin_coinc");
Eb_vs_Ein_1->Write("EBout_vsEin_coinc");
Ea_vs_Ein_2->Write("EAout_vsEin_passcut_coinc");
Eb_vs_Ein_2->Write("EBout_vsEin_passcut_coinc");
Et_vs_EionB_2->Write("EHout_vsEionB_passcut_coinc");
Et_vs_EionB_1->Write("EHout_vsEionB_coinc");
Et_vs_EionFid_2->Write("EHout_vsEionFid_passcut_coinc");
Et_vs_EionFid_1->Write("EHout_vsEionFid_coinc");
sigEa_vs_Ein_1->Write("sigma_EAout_vsEin_coinc");
sigEb_vs_Ein_1->Write("sigma_EBout_vsEin_coinc");
sigEa_vs_Ein_2->Write("sigma_EAout_vsEin_passcut_coinc");
sigEb_vs_Ein_2->Write("sigma_EBout_vsEin_passcut_coinc");
output_chi2->Close();
cst->cd(1);
hs_0->Draw("nostack");
hs_0->GetXaxis()->SetMoreLogLabels();
gPad->SetLogx();
gPad->RedrawAxis();
cst->cd(2);
hs_1->Draw("nostack");
hs_1->GetXaxis()->SetMoreLogLabels();
gPad->SetLogx();
gPad->RedrawAxis();
cst->cd(3);
hs_2->Draw("nostack");
hs_2->GetXaxis()->SetMoreLogLabels();
gPad->SetLogx();
gPad->RedrawAxis();
cst->cd(4);
leghs->Draw();
cst->SaveAs("Testing_output_2.png");
cst->Update();




TCanvas *cstB=new TCanvas("cstB","cstB",67,55,1853,1025);
cstB->Divide(2,2);

cstB->cd(1);
hs_0B->Draw("nostack");
hs_0B->GetXaxis()->SetMoreLogLabels();
gPad->SetLogx();
gPad->RedrawAxis();
cstB->cd(2);
hs_1B->Draw("nostack");
hs_1B->GetXaxis()->SetMoreLogLabels();
gPad->SetLogx();
gPad->RedrawAxis();
cstB->cd(3);
hs_2B->Draw("nostack");
hs_2B->GetXaxis()->SetMoreLogLabels();
gPad->SetLogx();
gPad->RedrawAxis();
cstB->cd(4);
leghs->Draw();
cstB->SaveAs("Testing_output_3.png");
cstB->Update();

if(verbose > 0) cout<<"canvas cnorm"<<endl;
TCanvas *cnorm=new TCanvas("cnorm","cnorm",67,55,1853,1025);
cnorm->cd(1);
hs_N->Draw("nostack");
hs_N->GetXaxis()->SetMoreLogLabels();
cnorm->SaveAs("Testing_output_1.png");
gPad->SetLogx();
gPad->RedrawAxis();
cnorm->Update();

if(saveLines){
   foutsaveLines->Close();
}
if(verbose > 0) cout<<"NFILES bugged :"<<NFILES_BUG<<endl;
for(int i=0;i<NFILES_BUG;i++){
	if(verbose > 0) cout<<filesbugged[i]<<endl;
}
//toto->Run();  
	return 0;
}   


//////////// check if file exists ////////////
bool FileExists(std::string const & filename)
{
    std::ifstream f(filename.c_str());
    if(f.good())
    {
        f.close();
        return true;
    }
    else
    {
        f.close();
        return false;
    }
} 

int color_rainbow(int const & i,int  const & nbval)
{
	int color;
    //100,798,845,64,52
    // 100 red
    // 64 bleu clair
    // 798 orane
    // 845 vert
    // (2 vert
    if(nbval==1){return 1;}
    double tempcolor = 51+ (100-51.)*i/(nbval-1);
	color = (int)tempcolor;
	
	
   if(nbval==3){
    	if(i==0){color=100;}
    	if(i==1){color=52;} 
    	if(i==2){color=64;}	
    }	
    
     if(nbval==4){
    	if(i==0){color=100;}
    	if(i==1){color=52;} 
    	if(i==2){color=64;}	
    	if(i==3){color=798;}
    }	
	
    return color;
} 


bool IsItaCoincidence(double const &  t,std::vector<double> const & vec,double const & tolerance)
{
	if(coincidencemethod=="slow"){
		int vectorsize=(int)vec.size();
		bool coincbool=false;
		
		// cout<<"size = "<<vectorsize<<endl;
		for(int i=0;i<vectorsize;i++){
			
			double timevec=vec[i];
			
			if(abs(t-timevec)<=tolerance){
				coincbool=true;
				break; // get out of loop if coinc found
			}else{
				if(timevec>t){break;} // get out of loop since time can only get higher so no chance of finding a coinc 
			}
			
		}

		return coincbool; // true means there is a coincidence 
	}else if(coincidencemethod=="fast"){
	
		auto first=lower_bound(vec.begin(),vec.end(),t-tolerance);
		auto last=upper_bound(vec.begin(),vec.end(),t+tolerance);
	
		if(first!=last){
			//cout<<"coinc found"<<endl; 
			return true;
		}else{ return false;}
	
	}else{cout<<"PROBLEM"<<endl; }
	
	return false;
}

bool IsItaCoincidenceAndAntiCoincidence(double const & t,std::vector<double>  const & veccoinc,std::vector<double> const & vecanticoinc,double const & tolerance)
{
	
	
	auto firstcoinc=lower_bound(veccoinc.begin(),veccoinc.end(),t-tolerance);
	auto  lastcoinc=upper_bound(veccoinc.begin(),veccoinc.end(),t+tolerance);
	bool coinc;
	if(firstcoinc!=lastcoinc){
		// cout<<"coinc found"<<endl; 
		coinc=true;
	
	}else{ coinc=false;}
	
	auto firstanticoinc=lower_bound(vecanticoinc.begin(),vecanticoinc.end(),t-tolerance);
	auto  lastanticoinc=upper_bound(vecanticoinc.begin(),vecanticoinc.end(),t+tolerance);
	bool anticoinc;
	if(firstanticoinc!=lastanticoinc){
		//if(verbose > 0) cout<<"coinc found"<<endl; 
		anticoinc=false;
	
	}else{ anticoinc=true;}
	
	bool coincanticoinc=coinc&&anticoinc;
	
	return coincanticoinc;
}

double GetBinomialUncertainty(int const &  ntrue,int const & nfalse)
{
	int N=ntrue+nfalse;
	double p=(double)ntrue/(double)N;
	double q=(1.-p);
	
	double uncertainty=TMath::Sqrt(N*p*q)/N;  // see binomialsimple.C on Desktop as a proof
	return uncertainty;
	
}
void linspace(double  const & min,double  const & max,int const & N,std::vector<double> &vec)
{
	// min always included 
	// if N>=2 both min and max always included 
	if(N==1){
		vec.push_back(min);
	}else if(N==2){ 
		vec.push_back(min);
		vec.push_back(max);
	}else{
		vec.push_back(min);
		int Nprime=N-2;
		double delta=(max-min)/(double)(Nprime+1);
		
		for(int i=1;i<=Nprime;i++){
			double value=min+delta*i;
			vec.push_back(value);
		}
		vec.push_back(max);
	}		
}

void logspace(double const &  min,double const & max,int const & N,std::vector<double> &vec)
{
	// min always included 
	// if N>=2 both min and max always included 
	if(N==1){
		vec.push_back(min);
	}else if(N==2){ 
		vec.push_back(min);
		vec.push_back(max);
	}else{
		vec.push_back(min);
		int Nprime=N-2;
		double logmax=TMath::Log10(max);
		double logmin=TMath::Log10(min);
		double deltalog=(logmax-logmin)/(double)(Nprime+1);
	
		for(int i=1;i<=Nprime;i++){
			double logvalue=logmin+deltalog*i;
			double value=TMath::Power(10,logvalue);
			vec.push_back(value);
		}
		vec.push_back(max);
	}		
}

bool AreTheseFilesDifferent(string const & file1,string const & file2)
{
	string command=Form("cmp %s %s | wc -l",file1.c_str(),file2.c_str());
	int N = gSystem->GetFromPipe(command.c_str()).Atoi(); // the bash command returns a TString that is converted into int with Atoi();
	if(N==1){
		return true; // Files are different	
	}else if(N==0){
		return false; // Files are the same
	}else{
		 cout<<"This shouldn't be possible"<<endl; exit(1);
	}
}

double CalculateDeltaT(double const & value,std::vector<double> const & vec)
{
	auto first=lower_bound(vec.begin(),vec.end(),value);
	//auto last=upper_bound(vec.begin(),vec.end(),value);

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


void buildefficiency(TH1D *h1,TH1D *h0,TH1D *he,int const & NBINS){

	int content_true=0;
	int content_false=0;
	int totalcontent=0;
	double error=0;
	double efficiency=0;
	
		for(int bin=1;bin<=NBINS;bin++){
			content_true = (int)h1->GetBinContent(bin);
			content_false= (int)h0->GetBinContent(bin);
			
			totalcontent=content_true+content_false;
			if(totalcontent==0){
				efficiency=0;
				error=0;
			}else{
				efficiency=content_true/(double)totalcontent;
				error=GetBinomialUncertainty(content_true,content_false);
			}
			he->SetBinContent(bin,efficiency);	
			he->SetBinError(bin,error);
			//if(efficiency>0){if(verbose > 0) cout<<"bin= "<<bin<<" eff="<<efficiency<<"	error="<<error<<endl;}
		}
}


void LoadBinary(string const & binaryfilepath,string const & run,int RunPart,string const & Channel,std::vector <string> const & Detector_List,double *amplitudes,int const & detID,int const & NB_Samples_to_read_per_channel,int const & stampstart)
{

	const int NDET=(int)Detector_List.size();
	int dictionnary[NDET];
	// example GeCo2 is abritrarily the first detector in TimingInformation, the first element of dictionnary will tell which channel it does correspond to in the binary file

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
	if (!fichierlog){if(verbose > 0) cout<< << " ==> ERROR: log file "<<logfile<<" does not exist !!" << endl;}
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
	
	// find the different channels 
	ifstream fichier;
	fichier.open(binaryfile, ios::binary);
	if (!fichier){ cout<< " ==> ERROR: binary file "<<binaryfile<<" does not exist !!" << endl;}
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
 
 
 	///////////////////// check dictionnary ////////////
 	// for(unsigned int p=0;p<Detector_List.size();p++){ if(verbose > 0) cout<<Detector_List[p]<<" "<<Channel<<" is "<<p<<" in Timing list but is "<<dictionnary[p]<<" in binary file"<<endl; }
 	
	fichier.seekg(0,ios::end);  //vec_octet_offset[RunPart]
	//long int end = fichier.tellg();	// fichier.tellg() returns the number of bytes between offset and end of file
	
	int Nchannels=(int)vec_Channels.size();
	
	int single_length = sizeof(short)*Nchannels;	// pour chaque coup d'horloge, nombre d'octets enregistres
    short *buffer = new short[single_length/sizeof(short)]; // tableau de taille Nchannels

	long int octet_start=octet_offset+(stampstart-NB_Samples_to_read_per_channel/2)*single_length;
   
	fichier.seekg(octet_start, std::ios::beg);


	for(int s=0;s<NB_Samples_to_read_per_channel;s++){
	
		fichier.read((char*)buffer,single_length);
		
		if(!memcmp(buffer,"* Arret ",8) || !memcmp(buffer,"* Arret ",3)){cout<<"End of Binary"<<endl; break;}
	
		amplitudes[s]=buffer[dictionnary[detID]];
	}

 	fichier.close();
}

void LoadFilteredBinary(string const & fileandpath,double *amplitudes,int const & NB_Samples_to_read_per_channel,int const & stampstart)
{
	ifstream fichier;
	fichier.open(fileandpath.c_str(), ios::binary);
	if (!fichier){ cout<< " ==> ERROR: binary file "<<fileandpath<<" does not exist !!" << endl;}
 	string line;
 
 	
	fichier.seekg(0,ios::end);  //vec_octet_offset[RunPart]
	//long int end = fichier.tellg();	// fichier.tellg() returns the number of bytes between offset and end of file
	
	int Nchannels=1;
	
	int single_length = sizeof(double)*Nchannels;	// pour chaque coup d'horloge, nombre d'octets enregistres
	//if(verbose > 0) cout<<"single_length="<<single_length<<endl;
    	//short *buffer = new short[single_length/sizeof(short)]; // tableau de taille Nchannels
    	double *buffer = new double[single_length/sizeof(double)]; // tableau de taille Nchannels
	//if(verbose > 0) cout<<"single_length/sizeof(short)="<<single_length/sizeof(short)<<endl;
	long int octet_offset=0;
	long int octet_start=octet_offset+(stampstart-NB_Samples_to_read_per_channel/2)*single_length;
   	//if(verbose > 0) cout<<"octet_offset ="<<octet_offset<<endl;
   	//if(verbose > 0) cout<<"octet_start ="<<octet_start<<" vs (stampstart-NB_Samples_to_read_per_channel/2)*single_length="<< (stampstart-NB_Samples_to_read_per_channel/2)*single_length<<" for stampstart="<<stampstart<<endl;
   	
	fichier.seekg(octet_start, std::ios::beg);
	 cout<<"where am I = "<<fichier.tellg()<<endl;

	for(int s=0;s<NB_Samples_to_read_per_channel;s++){
	
		fichier.read((char*)buffer,single_length);
		
		if(!memcmp(buffer,"* Arret ",8) || !memcmp(buffer,"* Arret ",3)){ cout<<"End of Binary"<<endl; break;}
	
		amplitudes[s]=buffer[0];
		//if(verbose > 0) cout<<"s="<<s<<"	ampl="<<amplitudes[s]<<endl;
	}

 	fichier.close();
	
	
}


void itx(double const & x=0.7,double const & y=0.7,TString const & mess="",int const & color=1, double const & size=0.06,int const & opt=0)
{
	TText *t=new TText(x,y,mess);
	if(opt==0){
	t->SetNDC(kTRUE);
	}
	t->SetTextColor(color);
	t->SetTextSize(size);
	t->Draw();
}

TCanvas* createcanvas(string const & namecanvas,string  const & style)
{
		TCanvas *c; //
	 c=new TCanvas(namecanvas.c_str(),namecanvas.c_str(),67,55,900,1025);
	 /*
    string commandheight="xrandr | awk '$0 ~ \"*\" {print $1}' | awk -F\"x\" '{print $2}'";
    string commandwidth="xrandr | awk '$0 ~ \"*\" {print $1}' | awk -F\"x\" '{print $1}'";
    int gScreenHeight = gSystem->GetFromPipe(commandheight.c_str()).Atoi(); // the bash command returns a
    int gScreenWidth = gSystem->GetFromPipe(commandwidth.c_str()).Atoi(); // the bash command returns a
    //if(verbose > 0) cout<<"screen :"<<gScreenWidth<<"   x   "<<gScreenHeight<<endl;

    if(style=="halfscreen"){
        int w=0.45*gScreenWidth;
        int h=0.95*gScreenHeight;
        c=new TCanvas(namecanvas.c_str(),namecanvas.c_str(),w,h);
    }else if(style=="fullscreen"){     
        int w=0.95*gScreenWidth;
        int h=0.95*gScreenHeight;
        c=new TCanvas(namecanvas.c_str(),namecanvas.c_str(),w,h);
    }else{
        if(verbose > 0) cout<<"this mode is not implemented"<<endl;
        c=new TCanvas(namecanvas.c_str(),namecanvas.c_str(),67,55,900,1025);
    }
*/
    return c;
}

void BuildEfficiencyIterativeCuts(TGraphErrors *geff,double *N_cut,double *N_all,double *Eline,int const & Nlines,string const & scale)
{
	geff->SetLineWidth(2);
	
	 const int N=Nlines;
	 double eff;
	 double eff_error;
	 double energy;
	 double energy_error;
	 for(int l=0;l<Nlines;l++){
	 	
	 	int Ntot=(int)N_all[l];
	 	int N_true=(int)N_cut[l];
	 	int N_false=Ntot-N_true;
	 	if(N_all[l]==0){
				eff=0;
				eff_error=0;
		}else{
				eff=N_cut[l]/N_all[l];
				eff_error=GetBinomialUncertainty(N_true,N_false);
		}
	 	
	 	
	 	if(scale=="eV"){
	 	energy=Eline[l]*1000.;
	 	}else{
	 	energy=Eline[l];
	 	}
	 	energy_error=0.;
	 	
	 	geff->SetPoint(l,energy,eff);
	 	geff->SetPointError(l,energy_error,eff_error);
	 	//if(verbose > 0) cout<<"E="<<energy<<" +- "<<energy_error<<"	Eff="<<eff<<" +- "<<eff_error<<endl;

	 }
	


}

void BuildEfficiencyIterativeCutsfromTH1D(TGraphErrors *geff_Trigger,TH1D *hcut,TH1D *hall)
{
	 double eff;
	 double eff_error;
	 double energy;
	 double energy_error;
	 geff_Trigger->SetLineWidth(2);
	 //if(verbose > 0) cout<<"HELLO IN"<<endl;
	for(int bin=1;bin<=hcut->GetNbinsX();bin++){
		energy=hcut->GetBinCenter(bin);
		double energy_check=hall->GetBinCenter(bin);
		energy_error=0.;
		if(energy!=energy_check){cout<<"serious problem"<<endl; exit(1);}
		
		int pt=bin-1;
		
		
		int Ntot=hall->GetBinContent(bin);
	 	int N_true=hcut->GetBinContent(bin);
	 	int N_false=Ntot-N_true;
	 	
	 	if(Ntot==0){
				eff=0;
				eff_error=0;
		}else{
				eff=(double)N_true/Ntot;
				eff_error=GetBinomialUncertainty(N_true,N_false);
		}
		geff_Trigger->SetPoint(pt,energy,eff);
		geff_Trigger->SetPointError(pt,energy_error,eff_error);
		//if(verbose > 0) cout<<"point "<<pt<<"	energy="<<energy<<" eff="<<eff<<endl;
	}
 //if(verbose > 0) cout<<"HELLO OUT"<<endl;

}

