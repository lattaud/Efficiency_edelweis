#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <istream>
#include <time.h>
#include <vector>
#include "TROOT.h"
#include "TSystem.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixTSym.h"
#include "TDecompLU.h"
#include "TApplication.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TGraph2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TArrow.h"
#include "TGraph.h"
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
#include "TEfficiency.h"
#include "TVirtualFFT.h"
#include "TGraphErrors.h"
#include "TComplex.h"
#include "TVirtualFFT.h"

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer.h"
#include "Math/GSLSimAnMinimizer.h"
#include "TColor.h"


using namespace std;

TColor *Color = new TColor();
Float_t R_vec[9] = {0.89,0.21,0.30,0.59,1.00,1.00,0.65,0.96,0.};
Float_t G_vec[9] = {0.10,0.49,0.68,0.30,0.50,0.80,0.34,0.51,0.};
Float_t B_vec[9] = {0.11,0.72,0.29,0.64,0.00,0.00,0.16,0.75,0.};

const int N_Max_Channel = 6;
const int N_Max_Ion = 4;
const int N_Max_Heat = 2;
int NumberOfPartition;
string DataListName;
int NFiles = 0;
vector<string> Filelist;

string path_processed, path_trigtraces, select_detector, Input_ProcessedData, Output_TemplateTraces, run, path_selectedtemplates;
int trigger_channel;

void GetInfoFromRunTree(int *Channel_ON, int &N_samples_Heat, int &N_samples_Ion, double &f_max_heat, double &f_max_ion, double &TimeWindow_Heat, double &TimeWindow_Ion, double *Polarity, int &N_tot_channel, int &N_Heat_Channel, int &N_Ion_Channel);
void EventSelection(int *Channel_ON, double *Polarity, int N_tot_channel, int N_Heat_Channel, int N_Ion_Channel, int N_samples_Heat, int N_samples_Ion, double f_max_heat);
int GetRootFileNames(string DataList);
double GetHeatCalibration(double GainA, double GainB, double *NonLinearityA, double *NonLinearityB, double FiducialVoltage, double *Amplitude_OF, string HeatChannel);
void GetIonCalibration(double GainIonA, double GainIonB, double GainIonC, double GainIonD, double CrossTalkAB, double CrossTalkBA, double CrossTalkAC, double CrossTalkCA, double CrossTalkAD, double CrossTalkDA, double CrossTalkBC, double CrossTalkCB, double CrossTalkBD, double CrossTalkDB, double CrossTalkCD, double CrossTalkDC, double *Amplitude_OF, double *Energy_Ion);
int Halfchi2_analysis;

int main(int argc, char *argv[]){
    
    TApplication toto("toto",0,0);
    gStyle->SetFillColor(10);
    gStyle->SetFrameFillColor(10);
    gStyle->SetCanvasColor(10);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetNdivisions(515, "xyz");
    gStyle->SetLabelOffset(0.005, "xyz");
    gStyle->SetLabelSize(0.05, "xyz");
    gStyle->SetTitleSize(0.05, "xyz");
    gStyle->SetOptStat(0);
    //gStyle->SetOptTitle(kFALSE);
    gStyle->SetTickLength(0.05, "xyz");
    gRandom->SetSeed(0);
    gStyle->SetCanvasDefH(800);
    gStyle->SetCanvasDefW(1600);
    
    path_processed = argv[1];                                 // /Users/billard/EDELWEISS/NEPAL_Stream/Processing_CC/ProcessedData
    path_trigtraces = argv[2];                                // /Users/billard/EDELWEISS/NEPAL_Stream/Processing_CC/Output_Trig
    path_selectedtemplates = argv[3];                         // /Users/billard/EDELWEISS/NEPAL_Stream/Processing_CC/SelectedTemplateEvent
    select_detector = argv[4];
    trigger_channel = atoi(argv[5]);
    DataListName = argv[6];
    Halfchi2_analysis = atoi(argv[7]);
    
    NFiles = GetRootFileNames(DataListName);
    
    system(Form("mkdir -p %s/%s",path_selectedtemplates.c_str(),DataListName.c_str()));
    
    int Channel_ON[N_Max_Channel], N_samples_Heat, N_samples_Ion, N_tot_channel, N_Heat_Channel, N_Ion_Channel;
    double f_max_heat, f_max_ion, TimeWindow_Heat, TimeWindow_Ion;
    double Polarity[N_Max_Ion];
    
    GetInfoFromRunTree(Channel_ON, N_samples_Heat, N_samples_Ion, f_max_heat, f_max_ion, TimeWindow_Heat, TimeWindow_Ion, Polarity, N_tot_channel, N_Heat_Channel, N_Ion_Channel);
    
    EventSelection(Channel_ON, Polarity, N_tot_channel, N_Heat_Channel, N_Ion_Channel, N_samples_Heat, N_samples_Ion, f_max_heat);
    
    cout<<"done"<<endl;
    toto.Run();
    return 0;
}


void GetInfoFromRunTree(int *Channel_ON, int &N_samples_Heat, int &N_samples_Ion, double &f_max_heat, double &f_max_ion, double &TimeWindow_Heat, double &TimeWindow_Ion, double *Polarity, int &N_tot_channel, int &N_Heat_Channel, int &N_Ion_Channel)
{
    string TemplateType = "Normal";

    int Chan_ON[N_Max_Channel], NumPart;
    double Polar_Ion[N_Max_Ion];

    TChain * RunTreeIn = new TChain(Form("RunTree_%s",TemplateType.c_str()));
    for(int i = 0; i<NFiles; i++)
    {
        RunTreeIn->Add(Form("%s/%s/%s/ProcessedData_%s_*_%s_ChanTrig%i.root",path_processed.c_str(),Filelist[i].c_str(),select_detector.c_str(),Filelist[i].c_str(),select_detector.c_str(),  trigger_channel));
        std::cout << " Adding file "<< Form("%s/%s/%s/ProcessedData_%s_*_%s_ChanTrig%i.root",path_processed.c_str(),Filelist[i].c_str(),select_detector.c_str(),Filelist[i].c_str(),select_detector.c_str(),  trigger_channel)<< std::endl;
    }
    
    RunTreeIn->SetBranchAddress("Chan_ON",&Chan_ON);
    RunTreeIn->SetBranchAddress("Polar_Ion",&Polar_Ion);
    RunTreeIn->SetBranchAddress("N_samples_Heat",&N_samples_Heat);
    RunTreeIn->SetBranchAddress("N_samples_Ion",&N_samples_Ion);
    RunTreeIn->SetBranchAddress("f_max_heat",&f_max_heat);
    RunTreeIn->SetBranchAddress("f_max_ion",&f_max_ion);
    RunTreeIn->SetBranchAddress("TimeWindow_Heat",&TimeWindow_Heat);
    RunTreeIn->SetBranchAddress("TimeWindow_Ion",&TimeWindow_Ion);
    RunTreeIn->SetBranchAddress("TimeWindow_Ion",&TimeWindow_Ion);
    RunTreeIn->SetBranchAddress("NumPart",&NumPart);
    
    NumberOfPartition = RunTreeIn->GetEntries();
    RunTreeIn->GetEntry(0);
    N_tot_channel = 0;
    N_Heat_Channel = 0;
    N_Ion_Channel = 0;
    for(int i = 0; i<N_Max_Channel; i++)
    {
        Channel_ON[i] = Chan_ON[i];
        if(i<N_Max_Ion)Polarity[i] = Polar_Ion[i];
        cout<<"\t"<<Channel_ON[i];
        N_tot_channel = N_tot_channel + Channel_ON[i];
        if(i<N_Max_Heat){N_Heat_Channel = N_Heat_Channel + Channel_ON[i];}
        else{N_Ion_Channel = N_Ion_Channel + Channel_ON[i];}
    }
    //if(select_detector == "NbSi209") N_Ion_Channel = 2 ;
    cout<<endl;
    cout<<"N_tot_channel = "<<N_tot_channel<<"\t"<<"N_Heat_Channel = "<<N_Heat_Channel<<"\t"<<"N_Ion_Channel = "<<N_Ion_Channel<<endl;
    
    cout<<"Heat: "<<N_samples_Heat<<"\t"<<f_max_heat<<"\t"<<TimeWindow_Heat<<endl;
    cout<<"Ion: "<<N_samples_Ion<<"\t"<<f_max_ion<<"\t"<<TimeWindow_Ion<<endl;
    cout<<"Total number of partitions = "<<NumberOfPartition<<endl;
    
    delete RunTreeIn;
    
    return;
}

void EventSelection(int *Channel_ON, double *Polarity, int N_tot_channel, int N_Heat_Channel, int N_Ion_Channel, int N_samples_Heat, int N_samples_Ion, double f_max_heat)
{
    string TemplateType = "Normal";
    string EventType = "trig";
    string filteringType;
    string filteringType_everest;
    string EventType_everest = "Trig";
    if(N_Ion_Channel>2){filteringType = "filt_decor"; filteringType_everest = "Filt_Decor";}
    else{filteringType = "filt"; filteringType_everest = "Filt";}
    int N_Max_Trace = 30;
    
    double Time_ind[N_samples_Heat], Time_Heat[N_samples_Heat];
    for(int i = 0;i<N_samples_Heat;i++){Time_ind[i] = i;}
    for(int i = 0; i<N_samples_Heat; i++){Time_Heat[i] = i*1./f_max_heat + 0.5/f_max_heat;}

    
    double Energy_OF[N_Max_Channel], chi2_OF[N_Max_Channel], chi2_OF_raw[N_Max_Channel], Amplitude_OF[N_Max_Channel], Amplitude_NEPAL_OF[N_Max_Channel], chi2_half[6];
    double Trace_Heat_A_Filt[N_samples_Heat], Trace_Heat_B_Filt[N_samples_Heat], Trace_Ion_A_Filt[N_samples_Heat], Trace_Ion_B_Filt[N_samples_Heat], Trace_Ion_C_Filt[N_samples_Heat], Trace_Ion_D_Filt[N_samples_Heat];
    double Trace_Heat_A_Filt_Decor[N_samples_Heat], Trace_Heat_B_Filt_Decor[N_samples_Heat], Trace_Ion_A_Filt_Decor[N_samples_Heat], Trace_Ion_B_Filt_Decor[N_samples_Heat], Trace_Ion_C_Filt_Decor[N_samples_Heat], Trace_Ion_D_Filt_Decor[N_samples_Heat];
    int Event_Number, MegaStp, MicroStp, NumRun, NumPart;
    double Trace_Max[N_Max_Channel], Trace_Min[N_Max_Channel];
    double Off[N_Max_Channel];
    
    
    double EhA, EhB, EiA, EiB, EiC, EiD, EiFID, EiVETO;
    double heatA, heatB, ionA, ionB, ionC, ionD;
    double heatChi2A, heatChi2B, ionChi2A, ionChi2B, ionChi2C, ionChi2D;
    
    int indexIonCalibData, indexHeatCalibData, indexParamCalibData, indexStrucCalibData;
    
    double GainA, GainB;
    
    double NonLinearityA[5], NonLinearityB[5];
    double FiducialVoltage;
    
    double GainIonA, GainIonB, GainIonC, GainIonD;
    double CrossTalkAB, CrossTalkBA, CrossTalkAC, CrossTalkCA, CrossTalkAD, CrossTalkDA, CrossTalkBC, CrossTalkCB, CrossTalkBD, CrossTalkDB, CrossTalkCD, CrossTalkDC;

    TChain * TreeIn = new TChain(Form("EventTree_%s_%s_%s",EventType.c_str(),TemplateType.c_str(),filteringType.c_str()));
    TChain * TreeRawIn = new TChain(Form("EventTree_%s_%s_%s",EventType.c_str(),TemplateType.c_str(),"raw"));
    TChain * TrigTrace = new TChain("tree");
    
    TChain * TreeEverestEnergy = new TChain(Form("Energies_%s_%s",EventType_everest.c_str(),filteringType_everest.c_str()));
    TChain * TreeEverestAmplitudes = new TChain(Form("Amplitudes_%s_%s",EventType_everest.c_str(),filteringType_everest.c_str()));
    TChain * TreeBoloHeader = new TChain("boloHeader");
    TChain * TreeIonCalibData = new TChain("ionCalibData");
    TChain * TreeHeatCalibData = new TChain("heatCalibData");
    TChain * TreeParamCalibData = new TChain("paramCalibData");

    for(int i = 0; i<NFiles; i++)
    {
        TreeIn->Add(Form("%s/%s/%s/ProcessedData_%s_*_%s_ChanTrig%i.root",path_processed.c_str(),Filelist[i].c_str(),select_detector.c_str(),Filelist[i].c_str(),select_detector.c_str(),  trigger_channel));
        TrigTrace->Add(Form("%s/%s/%s/RootFiles/TriggerData_%s_*_%s_ChanTrig%i.root",path_trigtraces.c_str(),Filelist[i].c_str(),select_detector.c_str(),Filelist[i].c_str(),select_detector.c_str(),  trigger_channel));
        TreeRawIn->Add(Form("%s/%s/%s/ProcessedData_%s_*_%s_ChanTrig%i.root",path_processed.c_str(),Filelist[i].c_str(),select_detector.c_str(),Filelist[i].c_str(),select_detector.c_str(),  trigger_channel));

        TreeEverestEnergy->Add(Form("%s/%s/%s/CalibratedData/CalibratedData_%s_*_%s_ChanTrig%i.root",path_processed.c_str(),Filelist[i].c_str(),select_detector.c_str(),Filelist[i].c_str(),select_detector.c_str(),  trigger_channel));
        TreeEverestAmplitudes->Add(Form("%s/%s/%s/CalibratedData/CalibratedData_%s_*_%s_ChanTrig%i.root",path_processed.c_str(),Filelist[i].c_str(),select_detector.c_str(),Filelist[i].c_str(),select_detector.c_str(),  trigger_channel));
        TreeBoloHeader->Add(Form("%s/%s/%s/CalibratedData/CalibratedData_%s_*_%s_ChanTrig%i.root",path_processed.c_str(),Filelist[i].c_str(),select_detector.c_str(),Filelist[i].c_str(),select_detector.c_str(),  trigger_channel));
        
        if(i==0)
        {
            TreeIonCalibData->Add(Form("%s/%s/%s/CalibratedData/CalibratedData_%s_S00_%s_ChanTrig%i.root",path_processed.c_str(),Filelist[i].c_str(),select_detector.c_str(),Filelist[i].c_str(),select_detector.c_str(),  trigger_channel));
            TreeHeatCalibData->Add(Form("%s/%s/%s/CalibratedData/CalibratedData_%s_S00_%s_ChanTrig%i.root",path_processed.c_str(),Filelist[i].c_str(),select_detector.c_str(),Filelist[i].c_str(),select_detector.c_str(),  trigger_channel));
            TreeParamCalibData->Add(Form("%s/%s/%s/CalibratedData/CalibratedData_%s_S00_%s_ChanTrig%i.root",path_processed.c_str(),Filelist[i].c_str(),select_detector.c_str(),Filelist[i].c_str(),select_detector.c_str(),  trigger_channel));
        }
        cout<<Form("%s/%s/%s/RootFiles/TriggerData_%s_*_%s_ChanTrig%i.root",path_trigtraces.c_str(),Filelist[i].c_str(),select_detector.c_str(),Filelist[i].c_str(),select_detector.c_str(),  trigger_channel)<<endl;
        cout<<Form("%s/%s/%s/CalibratedData/CalibratedData_%s_*_%s_ChanTrig%i.root",path_processed.c_str(),Filelist[i].c_str(),select_detector.c_str(),Filelist[i].c_str(),select_detector.c_str(),  trigger_channel)<<endl;
    }
    

    TreeIn->SetBranchAddress("Energy_OF",&Amplitude_NEPAL_OF);
    //if(select_detector == "NbSi209"){
   //   TreeIn -> SetBranchAddress ("chi2_OF_half",&chi2_OF);
   // }else{
      TreeIn->SetBranchAddress("chi2_OF",&chi2_OF);
    
   // }
    TreeIn->SetBranchAddress("Event_Number",&Event_Number);
    TreeIn->SetBranchAddress("MegaStp",&MegaStp);
    TreeIn->SetBranchAddress("MicroStp",&MicroStp);
    TreeIn->SetBranchAddress("NumRun",&NumRun);
    TreeIn->SetBranchAddress("NumPart",&NumPart);
    
    TreeRawIn->SetBranchAddress("Off",&Off);
    if(select_detector == "NbSi209"){
      TreeRawIn -> SetBranchAddress ("chi2_OF_half",&chi2_half);
    }
    TreeRawIn->SetBranchAddress("chi2_OF",&chi2_OF_raw);
    
    
      
    TrigTrace->SetBranchAddress("Trace_Heat_A_Filt",&Trace_Heat_A_Filt);
    TrigTrace->SetBranchAddress("Trace_Heat_B_Filt",&Trace_Heat_B_Filt);
    TrigTrace->SetBranchAddress("Trace_Ion_A_Filt",&Trace_Ion_A_Filt);
    TrigTrace->SetBranchAddress("Trace_Ion_B_Filt",&Trace_Ion_B_Filt);
    TrigTrace->SetBranchAddress("Trace_Ion_C_Filt",&Trace_Ion_C_Filt);
    TrigTrace->SetBranchAddress("Trace_Ion_D_Filt",&Trace_Ion_D_Filt);
    
    TrigTrace->SetBranchAddress("Trace_Heat_A_Filt_Decor",&Trace_Heat_A_Filt_Decor);
    TrigTrace->SetBranchAddress("Trace_Heat_B_Filt_Decor",&Trace_Heat_B_Filt_Decor);
    TrigTrace->SetBranchAddress("Trace_Ion_A_Filt_Decor",&Trace_Ion_A_Filt_Decor);
    TrigTrace->SetBranchAddress("Trace_Ion_B_Filt_Decor",&Trace_Ion_B_Filt_Decor);
    TrigTrace->SetBranchAddress("Trace_Ion_C_Filt_Decor",&Trace_Ion_C_Filt_Decor);
    TrigTrace->SetBranchAddress("Trace_Ion_D_Filt_Decor",&Trace_Ion_D_Filt_Decor);
    
    TreeEverestEnergy->SetBranchAddress("EhA",&EhA);
    TreeEverestEnergy->SetBranchAddress("EhB",&EhB);
    TreeEverestEnergy->SetBranchAddress("EiA",&EiA);
    TreeEverestEnergy->SetBranchAddress("EiB",&EiB);
    TreeEverestEnergy->SetBranchAddress("EiC",&EiC);
    TreeEverestEnergy->SetBranchAddress("EiD",&EiD);
    TreeEverestEnergy->SetBranchAddress("EiFID",&EiFID);
    TreeEverestEnergy->SetBranchAddress("EiVETO",&EiVETO);
    
    TreeEverestAmplitudes->SetBranchAddress("heatA",&heatA);
    TreeEverestAmplitudes->SetBranchAddress("heatB",&heatB);
    TreeEverestAmplitudes->SetBranchAddress("ionA",&ionA);
    TreeEverestAmplitudes->SetBranchAddress("ionB",&ionB);
    TreeEverestAmplitudes->SetBranchAddress("ionC",&ionC);
    TreeEverestAmplitudes->SetBranchAddress("ionD",&ionD);
    TreeEverestAmplitudes->SetBranchAddress("heatChi2A",&heatChi2A);
    TreeEverestAmplitudes->SetBranchAddress("heatChi2B",&heatChi2B);
    TreeEverestAmplitudes->SetBranchAddress("ionChi2A",&ionChi2A);
    TreeEverestAmplitudes->SetBranchAddress("ionChi2B",&ionChi2B);
    TreeEverestAmplitudes->SetBranchAddress("ionChi2C",&ionChi2C);
    TreeEverestAmplitudes->SetBranchAddress("ionChi2D",&ionChi2D);
    
    TreeBoloHeader->SetBranchAddress("indexIonCalibData",&indexIonCalibData);
    TreeBoloHeader->SetBranchAddress("indexHeatCalibData",&indexHeatCalibData);
    TreeBoloHeader->SetBranchAddress("indexParamCalibData",&indexParamCalibData);
    TreeBoloHeader->SetBranchAddress("indexStrucCalibData",&indexStrucCalibData);
    
    TreeHeatCalibData->SetBranchAddress("GainA",&GainA);
    TreeHeatCalibData->SetBranchAddress("GainB",&GainB);
    
    TreeParamCalibData->SetBranchAddress("NonLinearityA",&NonLinearityA);
    TreeParamCalibData->SetBranchAddress("NonLinearityB",&NonLinearityB);
    TreeParamCalibData->SetBranchAddress("FiducialVoltage",&FiducialVoltage);
    
    TreeIonCalibData->SetBranchAddress("GainA",&GainIonA);
    TreeIonCalibData->SetBranchAddress("GainB",&GainIonB);
    TreeIonCalibData->SetBranchAddress("GainC",&GainIonC);
    TreeIonCalibData->SetBranchAddress("GainD",&GainIonD);
    TreeIonCalibData->SetBranchAddress("CrossTalkAB",&CrossTalkAB);
    TreeIonCalibData->SetBranchAddress("CrossTalkBA",&CrossTalkBA);
    TreeIonCalibData->SetBranchAddress("CrossTalkAC",&CrossTalkAC);
    TreeIonCalibData->SetBranchAddress("CrossTalkCA",&CrossTalkCA);
    TreeIonCalibData->SetBranchAddress("CrossTalkAD",&CrossTalkAD);
    TreeIonCalibData->SetBranchAddress("CrossTalkDA",&CrossTalkDA);
    TreeIonCalibData->SetBranchAddress("CrossTalkBC",&CrossTalkBC);
    TreeIonCalibData->SetBranchAddress("CrossTalkCB",&CrossTalkCB);
    TreeIonCalibData->SetBranchAddress("CrossTalkBD",&CrossTalkBD);
    TreeIonCalibData->SetBranchAddress("CrossTalkDB",&CrossTalkDB);
    TreeIonCalibData->SetBranchAddress("CrossTalkCD",&CrossTalkCD);
    TreeIonCalibData->SetBranchAddress("CrossTalkDC",&CrossTalkDC);
    
    int N_event = TreeIn->GetEntries();
    cout<<"test N entries = "<<N_event<<"\t"<<TreeEverestEnergy->GetEntries()<<"\t"<<TreeEverestAmplitudes->GetEntries()<<endl;
    int N_trace = 0;
    int N_select = 0;
    int N_chi2 = 0;
    int col;
    double Param_x0[6] ;
    double Param_x1[6] ;
    double Param_x2[6] ;
    double NDOF_chi2 = N_samples_Heat;
    if(Halfchi2_analysis ==1)NDOF_chi2 = N_samples_Heat/2.;
    if(select_detector == "RED30"){
        
        Param_x0[0] = NDOF_chi2*1.15;
        Param_x1[0] = 22;  
        Param_x2[0] = 2;        
        Param_x0[1] = NDOF_chi2*1.15;
        Param_x1[1] = 0;
        Param_x2[1] = 0; 
        Param_x0[2] = N_samples_Ion*1.15;
        Param_x1[2] = 0;
        Param_x2[2] = 0; 
        Param_x0[3] = N_samples_Ion*1.4;
        Param_x1[3] = 0;
        Param_x2[3] = 0; 
        Param_x0[4] = N_samples_Ion*1.2;
        Param_x1[4] = 0;
        Param_x2[4] = 0; 
        Param_x0[5] = N_samples_Ion*1.4;
        Param_x1[5] = 0;
        Param_x2[5] = 0; 
    }
    
    if(select_detector == "NbSi209"){
        
        Param_x0[0] = NDOF_chi2*1.2;
        Param_x1[0] = 15;  
        Param_x2[0] = 2;        
        Param_x0[1] = NDOF_chi2*1.2;
        Param_x1[1] = 5;
        Param_x2[1] = 2; 
        Param_x0[2] = N_samples_Ion*1.2;
        Param_x1[2] = 0;
        Param_x2[2] = 0; 
        Param_x0[3] = N_samples_Ion*1.2;
        Param_x1[3] = 0;
        Param_x2[3] = 0; 
        Param_x0[4] = N_samples_Ion*10.;
        Param_x1[4] = 0;
        Param_x2[4] = 0; 
        Param_x0[5] = N_samples_Ion*10.;
        Param_x1[5] = 0;
        Param_x2[5] = 0; 
    }
    
    if(select_detector == "FID848"){
        
        Param_x0[0] = NDOF_chi2*1.15;
        Param_x1[0] = 22.;  
        Param_x2[0] = 5;        
        Param_x0[1] = NDOF_chi2*1.15;
        Param_x1[1] = 20.;
        Param_x2[1] = 5; 
        Param_x0[2] = N_samples_Ion*1.2;
        Param_x1[2] = 0;
        Param_x2[2] = 0; 
        Param_x0[3] = N_samples_Ion*1.2;
        Param_x1[3] = 0;
        Param_x2[3] = 0; 
        Param_x0[4] = N_samples_Ion*1.2;
        Param_x1[4] = 0;
        Param_x2[4] = 0; 
        Param_x0[5] = N_samples_Ion*1.2;
        Param_x1[5] = 0;
        Param_x2[5] = 0; 
    }
   

    TF1 * fun_chi2[N_Max_Channel];
    TGraph * gr_chi2VsAmp[N_Max_Channel];
    TGraph * gr_chi2VsAmp_select[N_Max_Channel];
    TGraph * gr_chi2VsAmp_goodioni[N_Max_Channel];
    TGraph * gr_chi2VsAmp_partioni[N_Max_Channel];
    TGraph * gr_chi2VsAmp_ionitail[N_Max_Channel];
    
    TGraph * gr_IonAVsIonB = new TGraph();
    TGraph * gr_IonCVsIonD = new TGraph();
    TGraph * gr_IonBVsIonD = new TGraph();
    TGraph * gr_BottomVsTop = new TGraph();
    TGraph * gr_IonAVsIonB_select = new TGraph();
    TGraph * gr_IonCVsIonD_select = new TGraph();
    TGraph * gr_IonBVsIonD_select = new TGraph();
    TGraph * gr_BottomVsTop_select = new TGraph();
    TGraph * gr_HeatVsBD = new TGraph();
    TGraph * gr_HeatVsBD_select = new TGraph();
    TGraph * gr_HeatVsBD_select_goodioni = new TGraph();
    TGraph * gr_HeatVsBD_select_partioni = new TGraph();
    TGraph * gr_HeatVsBD_select_ionitail = new TGraph();
    
    
    TGraph * gr_HeatAVsHeatB = new TGraph();
    TGraph * gr_HeatAVsHeatB_select = new TGraph();
    
        
    TF1 * fun_EionVsEheat_inf = new TF1("fun_EionVsEheat_inf","[0]+[1]*x",0,2000);
    TF1 * fun_EionVsEheat_sup = new TF1("fun_EionVsEheat_sup","[0]+[1]*x",0,2000);
    //fun_EionVsEheat_inf->SetParameters(-150/5.3,5.33);
    //fun_EionVsEheat_sup->SetParameters(-220/5.3,5.33);
    fun_EionVsEheat_inf->SetParameters(-5,2.);
    fun_EionVsEheat_sup->SetParameters(5,2.);
    TH1F * h_Energy_All[N_Max_Channel];
    TH1F * h_Energy_select[N_Max_Channel];

    for(int i = 0; i<N_Max_Channel; i++)
    {
        gr_chi2VsAmp[i]          = new TGraph();
        gr_chi2VsAmp_select[i]   = new TGraph();
        gr_chi2VsAmp_goodioni[i] = new TGraph();
        gr_chi2VsAmp_partioni[i] = new TGraph();
        gr_chi2VsAmp_ionitail[i] = new TGraph();
        Trace_Max[i] = -1e6;
        Trace_Min[i] = 1e6;
        //fun_chi2[i] = new TF1(Form("fun_chi2_%i",i),"[0]+[1]*x*x",0,1e3);
        fun_chi2[i] = new TF1(Form("fun_chi2_%i",i),"[0]+100*[3]*TMath::Power((fabs(x)/[1]), [2])",0,1e3);
        double param_4  = NDOF_chi2 ;
        if(i > 1 ) param_4 = 0. ;
        fun_chi2[i]->SetParameters(Param_x0[i],Param_x1[i],Param_x2[i],param_4);
        h_Energy_All[i] = new TH1F(Form("h_Energy_All%i",i),Form("h_Energy_All%i",i),1000,0,15);
        h_Energy_select[i] = new TH1F(Form("h_Energy_select%i",i),Form("h_Energy_select%i",i),1000,0,15);

    }
    int Nchannel_trace = N_Max_Channel ;
   // if(select_detector == "NbSi209")  Nchannel_trace = 4 ;
    
    TGraph * gr_Trace[Nchannel_trace][N_Max_Trace];
    
    bool OK;
    string partNumber;
    double Gain_ADUperkeV[N_Max_Channel];
    
    TFile fileOut(Form("%s/%s/TemplateEvents_%s.root",path_selectedtemplates.c_str(),DataListName.c_str(),DataListName.c_str()),"recreate");
    TTree * tree = new TTree("tree","");

    tree->Branch("Trace_Heat_A_Filt",&Trace_Heat_A_Filt,Form("Trace_Heat_A_Filt[%i]/D",N_samples_Heat));
    tree->Branch("Trace_Heat_B_Filt",&Trace_Heat_B_Filt,Form("Trace_Heat_B_Filt[%i]/D",N_samples_Heat));
    tree->Branch("Trace_Ion_A_Filt",&Trace_Ion_A_Filt,Form("Trace_Ion_A_Filt[%i]/D",N_samples_Ion));
    tree->Branch("Trace_Ion_B_Filt",&Trace_Ion_B_Filt,Form("Trace_Ion_B_Filt[%i]/D",N_samples_Ion));
    tree->Branch("Trace_Ion_C_Filt",&Trace_Ion_C_Filt,Form("Trace_Ion_C_Filt[%i]/D",N_samples_Ion));
    tree->Branch("Trace_Ion_D_Filt",&Trace_Ion_D_Filt,Form("Trace_Ion_D_Filt[%i]/D",N_samples_Ion));
    tree->Branch("Trace_Heat_A_Filt_Decor",&Trace_Heat_A_Filt_Decor,Form("Trace_Heat_A_Filt_Decor[%i]/D",N_samples_Heat));
    tree->Branch("Trace_Heat_B_Filt_Decor",&Trace_Heat_B_Filt_Decor,Form("Trace_Heat_B_Filt_Decor[%i]/D",N_samples_Heat));
    tree->Branch("Trace_Ion_A_Filt_Decor",&Trace_Ion_A_Filt_Decor,Form("Trace_Ion_A_Filt_Decor[%i]/D",N_samples_Ion));
    tree->Branch("Trace_Ion_B_Filt_Decor",&Trace_Ion_B_Filt_Decor,Form("Trace_Ion_B_Filt_Decor[%i]/D",N_samples_Ion));
    tree->Branch("Trace_Ion_C_Filt_Decor",&Trace_Ion_C_Filt_Decor,Form("Trace_Ion_C_Filt_Decor[%i]/D",N_samples_Ion));
    tree->Branch("Trace_Ion_D_Filt_Decor",&Trace_Ion_D_Filt_Decor,Form("Trace_Ion_D_Filt_Decor[%i]/D",N_samples_Ion));
    tree->Branch("Gain_ADUperkeV",&Gain_ADUperkeV,Form("Gain_ADUperkeV[%i]/D",N_Max_Channel));
    tree->Branch("Amplitude_OF",&Amplitude_OF,Form("Amplitude_OF[%i]/D",N_Max_Channel));
    tree->Branch("Energy_OF",&Energy_OF,Form("Energy_OF[%i]/D",N_Max_Channel));
    tree->Branch("NumRun",&NumRun,"NumRun/I");
    tree->Branch("NumPart",&NumPart,"NumPart/I");
    
    double Energy_heat_A, Energy_heat_B;
    double Energy_Ion_test[4];

    for(int i = 0; i<N_event; i++)
    {
        TreeIn->GetEntry(i);
        TreeRawIn->GetEntry(i);

        TreeEverestAmplitudes->GetEntry(i);
        Amplitude_OF[0] = heatA;
        Amplitude_OF[1] = heatB;
        Amplitude_OF[2] = ionA;
        Amplitude_OF[3] = ionB;
        Amplitude_OF[4] = ionC;
        Amplitude_OF[5] = ionD;
        if(Halfchi2_analysis == 1 ){
         chi2_OF_raw[0] = chi2_half[0];//heatChi2A;
         chi2_OF_raw[1] = chi2_half[1];//heatChi2B;
        }else{
         chi2_OF_raw[0] = heatChi2A;
         chi2_OF_raw[1] = heatChi2B;        
        }
        chi2_OF_raw[2] = ionChi2A;
        chi2_OF_raw[3] = ionChi2B;
        chi2_OF_raw[4] = ionChi2C;
        chi2_OF_raw[5] = ionChi2D;
        
       // std::cout << "chi2 half A "<<  heatChi2A << " B "<< heatChi2B << "Ndof " <<N_samples_Heat << std::endl;
        
        
        TreeEverestEnergy->GetEntry(i);
        Energy_OF[0] = EhA;
        Energy_OF[1] = EhB;
        Energy_OF[2] = EiA;
        Energy_OF[3] = EiB;
        Energy_OF[4] = EiC;
        Energy_OF[5] = EiD;
        
        TreeBoloHeader->GetEntry(i);
        TreeParamCalibData->GetEntry(indexParamCalibData);
        TreeHeatCalibData->GetEntry(indexHeatCalibData);
        TreeIonCalibData->GetEntry(indexIonCalibData);
        
        Energy_heat_A = GetHeatCalibration(GainA, GainB, NonLinearityA, NonLinearityB, FiducialVoltage, Amplitude_OF, "HeatA");
        Energy_heat_B = GetHeatCalibration(GainA, GainB, NonLinearityA, NonLinearityB, FiducialVoltage, Amplitude_OF, "HeatB");
        GetIonCalibration(GainIonA, GainIonB, GainIonC, GainIonD, CrossTalkAB, CrossTalkBA, CrossTalkAC, CrossTalkCA, CrossTalkAD, CrossTalkDA, CrossTalkBC, CrossTalkCB, CrossTalkBD, CrossTalkDB, CrossTalkCD, CrossTalkDC, Amplitude_OF, Energy_Ion_test);

        for(int j = 0;j<N_Max_Channel;j++)
        {
            Gain_ADUperkeV[j] = 1.;
            if(Channel_ON[j] == 1)
            {
                if(j>=N_Max_Heat){if(Polarity[j-N_Max_Heat]>0){Energy_OF[j] = Energy_OF[j];}}
                gr_chi2VsAmp[j]->SetPoint(i,Energy_OF[j],chi2_OF_raw[j]);
                Gain_ADUperkeV[j] = Amplitude_OF[j]/Energy_OF[j];
            }
        }

        OK = true;
        for(int i = 0;i<N_Max_Channel;i++){if(Channel_ON[i]==1 && chi2_OF_raw[i]>fun_chi2[i]->Eval(Energy_OF[i])){OK=false;}}
        for(int j = 0;j<N_Max_Channel;j++){if(Channel_ON[j]==1 && fabs(Off[j])>15000){OK=false;}}
        if(Event_Number == 92  && select_detector == "FID848" )OK=false;
        if(select_detector == "FID848" && fabs(Energy_OF[0]- Energy_OF[1]) > 0.5 )OK=false;
        if(OK)
        {
            gr_IonAVsIonB->SetPoint(N_chi2,Energy_OF[2],Energy_OF[3]);
            gr_IonCVsIonD->SetPoint(N_chi2,Energy_OF[4],Energy_OF[5]);
            gr_IonBVsIonD->SetPoint(N_chi2,Energy_OF[3],Energy_OF[5]);
            gr_BottomVsTop->SetPoint(N_chi2,Energy_OF[2]+Energy_OF[3],Energy_OF[4]+Energy_OF[5]);
            gr_HeatVsBD    ->SetPoint(N_chi2,Energy_OF[0],EiA+EiB+EiC+EiD);
            gr_HeatAVsHeatB->SetPoint(N_chi2,Energy_OF[0],Energy_OF[1]);
            for(int j = 0;j<N_Max_Channel;j++){if(Channel_ON[j]==1){h_Energy_All[j]->Fill(Energy_OF[j]);}}
            N_chi2++;
        }
        /*
        if(N_Ion_Channel>0)
        {
            if(!((EiA+EiB+EiC+EiD)<fun_EionVsEheat_inf->Eval(EhA)) && ((EiA+EiB+EiC+EiD)>fun_EionVsEheat_sup->Eval(EhA))){OK = false;}
        }
        */
       

        double cut_E_up   = 12.6;
        double cut_E_down =  1.5;
        if(select_detector == "NbSi209"){
            cut_E_up   = 12.6;
            cut_E_down =  2.0;
        }else if (select_detector == "RED30"){
        
            cut_E_up   = 11.;
            cut_E_down =  1.5;
        
        }else if (select_detector == "FID848"){
            
            cut_E_up   = 11.;
            cut_E_down =  2.;
        
        } 
        
        if(OK && Energy_OF[0]>cut_E_down && Energy_OF[0]<cut_E_up)
        {
            TrigTrace->GetEntry(i);
            
            for(int j = 0;j<N_Max_Channel;j++)
            {
                if(Channel_ON[j] == 1){
                     gr_chi2VsAmp_select[j]   ->SetPoint(N_select,Energy_OF[j],chi2_OF_raw[j]);
                     if(Energy_OF[0] >= 9. && Energy_OF[0] <= 12.7 && EiA+EiB+EiC+EiD <= 22. && EiA+EiB+EiC+EiD >= 18.) gr_chi2VsAmp_goodioni[j] ->SetPoint(N_select,Energy_OF[j],chi2_OF_raw[j]);
                     if(Energy_OF[0] >= 9. && Energy_OF[0] <= 12.7 && EiA+EiB+EiC+EiD <= 12. && EiA+EiB+EiC+EiD >= 9.)  gr_chi2VsAmp_partioni[j] ->SetPoint(N_select,Energy_OF[j],chi2_OF_raw[j]);
                     if(Energy_OF[0] >= 2. && Energy_OF[0] <= 9. && EiA+EiB+EiC+EiD <= 22. && EiA+EiB+EiC+EiD >= 1.5)   gr_chi2VsAmp_ionitail[j] ->SetPoint(N_select,Energy_OF[j],chi2_OF_raw[j]);
                }
            }
            gr_IonAVsIonB_select->SetPoint(N_select,Energy_OF[2],Energy_OF[3]);
            gr_IonCVsIonD_select->SetPoint(N_select,Energy_OF[4],Energy_OF[5]);
            gr_IonBVsIonD_select->SetPoint(N_select,Energy_OF[3],Energy_OF[5]);
            gr_BottomVsTop_select->SetPoint(N_select,Energy_OF[2]+Energy_OF[3],Energy_OF[4]+Energy_OF[5]);
            gr_HeatVsBD_select->SetPoint(N_select,Energy_OF[0],EiA+EiB+EiC+EiD);
            gr_HeatAVsHeatB_select->SetPoint(N_chi2,Energy_OF[0],Energy_OF[1]);
            if(Energy_OF[0] >= 9. && Energy_OF[0] <= 12.7 && EiA+EiB+EiC+EiD <= 22. && EiA+EiB+EiC+EiD >= 18.) gr_HeatVsBD_select_goodioni->SetPoint(N_select,Energy_OF[0],EiA+EiB+EiC+EiD);
            if(Energy_OF[0] >= 9. && Energy_OF[0] <= 12.7 && EiA+EiB+EiC+EiD <= 12. && EiA+EiB+EiC+EiD >= 9.)  gr_HeatVsBD_select_partioni->SetPoint(N_select,Energy_OF[0],EiA+EiB+EiC+EiD);
            if(Energy_OF[0] >= 2. && Energy_OF[0] <= 9. && EiA+EiB+EiC+EiD <= 22. && EiA+EiB+EiC+EiD >= 1.5)   gr_HeatVsBD_select_ionitail->SetPoint(N_select,Energy_OF[0],EiA+EiB+EiC+EiD);
            
            for(int j = 0;j<N_Max_Channel;j++){if(Channel_ON[j]==1){h_Energy_select[j]->Fill(Energy_OF[j]);}}
            N_select++;
            
            
            if(N_trace<N_Max_Trace)
            {
                if(N_Ion_Channel>2)
                {
                    if(Channel_ON[0] == 1)gr_Trace[0][N_trace] = new TGraph(N_samples_Heat,Time_ind,Trace_Heat_A_Filt_Decor);
                    if( Trace_Heat_A_Filt_Decor[50] >= 40.  ) std::cout<< " Event "<< i <<" Energy A "<<Energy_OF[0]<< " chi2 A "<< chi2_OF_raw[0]/1024. <<" Energy B "<<Energy_OF[1]<< " chi2 B "<< chi2_OF_raw[1]/1024.<<" event number "<< Event_Number<<" run "<< NumRun<<" partition "<<NumPart<<std::endl;
                    if(Channel_ON[1] == 1)gr_Trace[1][N_trace] = new TGraph(N_samples_Heat,Time_ind,Trace_Heat_B_Filt_Decor);
                    if(Channel_ON[2] == 1)gr_Trace[2][N_trace] = new TGraph(N_samples_Heat,Time_ind,Trace_Ion_A_Filt_Decor);
                    if(Channel_ON[3] == 1)gr_Trace[3][N_trace] = new TGraph(N_samples_Heat,Time_ind,Trace_Ion_B_Filt_Decor);
                    if(Channel_ON[4] == 1)gr_Trace[4][N_trace] = new TGraph(N_samples_Heat,Time_ind,Trace_Ion_C_Filt_Decor);
                    if(Channel_ON[5] == 1 )gr_Trace[5][N_trace] = new TGraph(N_samples_Heat,Time_ind,Trace_Ion_D_Filt_Decor);
                    
                    if(TMath::MaxElement(N_samples_Heat,Trace_Heat_A_Filt_Decor)>=Trace_Max[0]){Trace_Max[0] = TMath::MaxElement(N_samples_Heat,Trace_Heat_A_Filt_Decor);}
                    if(TMath::MinElement(N_samples_Heat,Trace_Heat_A_Filt_Decor)<=Trace_Min[0]){Trace_Min[0] = TMath::MinElement(N_samples_Heat,Trace_Heat_A_Filt_Decor);}
                    if(TMath::MaxElement(N_samples_Heat,Trace_Heat_B_Filt_Decor)>=Trace_Max[1]){Trace_Max[1] = TMath::MaxElement(N_samples_Heat,Trace_Heat_B_Filt_Decor);}
                    if(TMath::MinElement(N_samples_Heat,Trace_Heat_B_Filt_Decor)<=Trace_Min[1]){Trace_Min[1] = TMath::MinElement(N_samples_Heat,Trace_Heat_B_Filt_Decor);}
                    if(TMath::MaxElement(N_samples_Heat,Trace_Ion_A_Filt_Decor)>=Trace_Max[2]){Trace_Max[2] = TMath::MaxElement(N_samples_Heat,Trace_Ion_A_Filt_Decor);}
                    if(TMath::MinElement(N_samples_Heat,Trace_Ion_A_Filt_Decor)<=Trace_Min[2]){Trace_Min[2] = TMath::MinElement(N_samples_Heat,Trace_Ion_A_Filt_Decor);}
                    if(TMath::MaxElement(N_samples_Heat,Trace_Ion_B_Filt_Decor)>=Trace_Max[3]){Trace_Max[3] = TMath::MaxElement(N_samples_Heat,Trace_Ion_B_Filt_Decor);}
                    if(TMath::MinElement(N_samples_Heat,Trace_Ion_B_Filt_Decor)<=Trace_Min[3]){Trace_Min[3] = TMath::MinElement(N_samples_Heat,Trace_Ion_B_Filt_Decor);}
                    if(TMath::MaxElement(N_samples_Heat,Trace_Ion_C_Filt_Decor)>=Trace_Max[4]){Trace_Max[4] = TMath::MaxElement(N_samples_Heat,Trace_Ion_C_Filt_Decor);}
                    if(TMath::MinElement(N_samples_Heat,Trace_Ion_C_Filt_Decor)<=Trace_Min[4]){Trace_Min[4] = TMath::MinElement(N_samples_Heat,Trace_Ion_C_Filt_Decor);}
                    if(TMath::MaxElement(N_samples_Heat,Trace_Ion_D_Filt_Decor)>=Trace_Max[5]){Trace_Max[5] = TMath::MaxElement(N_samples_Heat,Trace_Ion_D_Filt_Decor);}
                    if(TMath::MinElement(N_samples_Heat,Trace_Ion_D_Filt_Decor)<=Trace_Min[5]){Trace_Min[5] = TMath::MinElement(N_samples_Heat,Trace_Ion_D_Filt_Decor);}
                }
                else
                {
                    if(Channel_ON[0] == 1)gr_Trace[0][N_trace] = new TGraph(N_samples_Heat,Time_ind,Trace_Heat_A_Filt);
                    if(Channel_ON[1] == 1)gr_Trace[1][N_trace] = new TGraph(N_samples_Heat,Time_ind,Trace_Heat_B_Filt);
                    if(Channel_ON[2] == 1)gr_Trace[2][N_trace] = new TGraph(N_samples_Heat,Time_ind,Trace_Ion_A_Filt);
                    if(Channel_ON[3] == 1)gr_Trace[3][N_trace] = new TGraph(N_samples_Heat,Time_ind,Trace_Ion_B_Filt);
                    if(Channel_ON[4] == 1 /*&& select_detector != "NbSi209"*/)gr_Trace[4][N_trace] = new TGraph(N_samples_Heat,Time_ind,Trace_Ion_C_Filt);
                    if(Channel_ON[5] == 1 /*&& select_detector != "NbSi209"*/)gr_Trace[5][N_trace] = new TGraph(N_samples_Heat,Time_ind,Trace_Ion_D_Filt);
                    
                    if(TMath::MaxElement(N_samples_Heat,Trace_Heat_A_Filt)>=Trace_Max[0]){Trace_Max[0] = TMath::MaxElement(N_samples_Heat,Trace_Heat_A_Filt);}
                    if(TMath::MinElement(N_samples_Heat,Trace_Heat_A_Filt)<=Trace_Min[0]){Trace_Min[0] = TMath::MinElement(N_samples_Heat,Trace_Heat_A_Filt);}
                    if(TMath::MaxElement(N_samples_Heat,Trace_Heat_B_Filt)>=Trace_Max[1]){Trace_Max[1] = TMath::MaxElement(N_samples_Heat,Trace_Heat_B_Filt);}
                    if(TMath::MinElement(N_samples_Heat,Trace_Heat_B_Filt)<=Trace_Min[1]){Trace_Min[1] = TMath::MinElement(N_samples_Heat,Trace_Heat_B_Filt);}
                    if(TMath::MaxElement(N_samples_Heat,Trace_Ion_A_Filt)>=Trace_Max[2]){Trace_Max[2] = TMath::MaxElement(N_samples_Heat,Trace_Ion_A_Filt);}
                    if(TMath::MinElement(N_samples_Heat,Trace_Ion_A_Filt)<=Trace_Min[2]){Trace_Min[2] = TMath::MinElement(N_samples_Heat,Trace_Ion_A_Filt);}
                    if(TMath::MaxElement(N_samples_Heat,Trace_Ion_B_Filt)>=Trace_Max[3]){Trace_Max[3] = TMath::MaxElement(N_samples_Heat,Trace_Ion_B_Filt);}
                    if(TMath::MinElement(N_samples_Heat,Trace_Ion_B_Filt)<=Trace_Min[3]){Trace_Min[3] = TMath::MinElement(N_samples_Heat,Trace_Ion_B_Filt);}
                    if(TMath::MaxElement(N_samples_Heat,Trace_Ion_C_Filt)>=Trace_Max[4]){Trace_Max[4] = TMath::MaxElement(N_samples_Heat,Trace_Ion_C_Filt);}
                    if(TMath::MinElement(N_samples_Heat,Trace_Ion_C_Filt)<=Trace_Min[4]){Trace_Min[4] = TMath::MinElement(N_samples_Heat,Trace_Ion_C_Filt);}
                    if(TMath::MaxElement(N_samples_Heat,Trace_Ion_D_Filt)>=Trace_Max[5]){Trace_Max[5] = TMath::MaxElement(N_samples_Heat,Trace_Ion_D_Filt);}
                    if(TMath::MinElement(N_samples_Heat,Trace_Ion_D_Filt)<=Trace_Min[5]){Trace_Min[5] = TMath::MinElement(N_samples_Heat,Trace_Ion_D_Filt);}
                }
                
                N_trace++;
            }
            tree->Fill();
        }
    }
    
    for(int i = 0; i<N_trace;i++)
    {
        col = (int) (51. + 49./(N_trace-1)*i);
        for(int j = 0; j<Nchannel_trace; j++){if(Channel_ON[j] == 1)gr_Trace[j][i]->SetLineColor(col);}
    }
    cout<<"Total number of trig events = "<<N_event<<endl;
    cout<<"Number of good trig events = "<<N_chi2<<endl;
    cout<<"total number of template events = "<<N_select<<endl;

    
    TH2F * h2_Chi2VsAmp[N_Max_Channel];
    string chanName;
    TCanvas * cc = new TCanvas();
    cc->Divide(2,3);
    for(int i = 0; i<N_Max_Channel; i++)
    {
        if(i==0)chanName = "Chal A";
        if(i==1)chanName = "Chal B";
        if(i==2)chanName = "Ion A";
        if(i==3)chanName = "Ion B";
        if(i==4)chanName = "Ion C";
        if(i==5)chanName = "Ion D";
        h2_Chi2VsAmp[i] = new TH2F(Form("%s",chanName.c_str()),Form("%s",chanName.c_str()),1000,0.01,300,1000,5e2,1e6);
        h2_Chi2VsAmp[i]->GetXaxis()->SetTitle("Energy [keVee]");
        h2_Chi2VsAmp[i]->GetYaxis()->SetTitle("Chi2 [a.u.]");
        
        cc->cd(i+1);
        h2_Chi2VsAmp[i]->Draw();
        if(Channel_ON[i] == 1)
        {
            gr_chi2VsAmp[i]->Draw("samep");
            gr_chi2VsAmp[i]->SetMarkerStyle(6);
            gr_chi2VsAmp_select[i]->Draw("samep");
            gr_chi2VsAmp_select[i]->SetMarkerStyle(6);
            gr_chi2VsAmp_select[i]->SetMarkerColor(2);
            fun_chi2[i]->Draw("samel");
            fun_chi2[i]->SetLineColor(4);
        }
        cc->cd(i+1)->SetGrid();
        cc->cd(i+1)->SetLogy();
        cc->cd(i+1)->SetLogx();
        //cc->cd(i+1)->Update();
        //cc->cd(i+1)->Update();
    }
 //   cc->Update();
    cc->SaveAs(Form("%s/%s/Chi2Selection_%s.png",path_selectedtemplates.c_str(),DataListName.c_str(),DataListName.c_str()));
  /*  
    for(int i = 0; i<N_Max_Channel; i++)
    {
        if(i==0)chanName = "Chal A";
        if(i==1)chanName = "Chal B";
        if(i==2)chanName = "Ion A";
        if(i==3)chanName = "Ion B";
        if(i==4)chanName = "Ion C";
        if(i==5)chanName = "Ion D";
        h2_Chi2VsAmp[i] = new TH2F(Form("%s",chanName.c_str()),Form("%s",chanName.c_str()),1000,0.01,300,1000,5e2,1e6);
        h2_Chi2VsAmp[i]->GetXaxis()->SetTitle("Energy [keVee]");
        h2_Chi2VsAmp[i]->GetYaxis()->SetTitle("Chi2 [a.u.]");
        
        cc->cd(i+1);
        h2_Chi2VsAmp[i]->Draw();
        if(Channel_ON[i] == 1)
        {
            gr_chi2VsAmp_goodioni[i] ->SetMarkerColor(kRed+1);
            gr_chi2VsAmp_partioni[i] ->SetMarkerColor(kBlue+1);
            gr_chi2VsAmp_ionitail[i] ->SetMarkerColor(kGreen +1);

            gr_chi2VsAmp_goodioni[i] ->SetMarkerStyle(6);
            gr_chi2VsAmp_partioni[i] ->SetMarkerStyle(6);
            gr_chi2VsAmp_ionitail[i] ->SetMarkerStyle(6);


            gr_chi2VsAmp_goodioni[i] ->Draw("samep");
            gr_chi2VsAmp_partioni[i] ->Draw("samep");
            gr_chi2VsAmp_ionitail[i] ->Draw("samep");

            fun_chi2[i]->Draw("samel");
            fun_chi2[i]->SetLineColor(4);
        }
        cc->cd(i+1)->SetGrid();
        cc->cd(i+1)->SetLogy();
        cc->cd(i+1)->SetLogx();
        cc->cd(i+1)->Update();
        cc->cd(i+1)->Update();
    }
    */
    //cc->Update();
    //cc->SaveAs(Form("%s/%s/Chi2Selected_population_%s.png",path_selectedtemplates.c_str(),DataListName.c_str(),DataListName.c_str()));

    TH2F * h2_Trace[N_Max_Channel];
    TCanvas * ccTrace = new TCanvas();
    ccTrace->Divide(2,3);
    for(int i = 0; i<Nchannel_trace; i++)
    {
        if(i==0)chanName = "Chal A";
        if(i==1)chanName = "Chal B";
        if(i==2)chanName = "Ion A";
        if(i==3)chanName = "Ion B";
        if(i==4)chanName = "Ion C";
        if(i==5)chanName = "Ion D";
        h2_Trace[i] = new TH2F(Form("Trace %s",chanName.c_str()),Form("%s",chanName.c_str()),100,0,N_samples_Heat,10000,-30000,30000);
        h2_Trace[i]->GetXaxis()->SetTitle("Time [ADU]");
        h2_Trace[i]->GetYaxis()->SetTitle("Amplitude [ADU]");
        
        ccTrace->cd(i+1);
        h2_Trace[i]->Draw();
        if(Channel_ON[i] == 1)
        {
            for(int j=0;j<N_trace;j++)
            {
                gr_Trace[i][j]->Draw("samel");
                //std::cout<<" test "<<std::endl;
            }
            h2_Trace[i]->GetYaxis()->SetRangeUser(Trace_Min[i] - (Trace_Max[i] - Trace_Min[i])*0.1,Trace_Max[i] + (Trace_Max[i] - Trace_Min[i])*0.1);
        }
        ccTrace->cd(i+1)->SetGrid();
      //  ccTrace->cd(i+1)->Update();
    }
   // ccTrace->Update();
    ccTrace->SaveAs(Form("%s/%s/TraceSelection_%s.png",path_selectedtemplates.c_str(),DataListName.c_str(),DataListName.c_str()));
    std::cout<<Form("%s/%s/TraceSelection_%s.png",path_selectedtemplates.c_str(),DataListName.c_str(),DataListName.c_str())<<std::endl;
    TCanvas * ccEnergy = new TCanvas();
    ccEnergy->Divide(2,3);
    for(int i = 0; i<N_Max_Channel; i++)
    {
        if(i==0)chanName = "Chal A";
        if(i==1)chanName = "Chal B";
        if(i==2)chanName = "Ion A";
        if(i==3)chanName = "Ion B";
        if(i==4)chanName = "Ion C";
        if(i==5)chanName = "Ion D";
        
        ccEnergy->cd(i+1);
        h_Energy_All[i]->Draw();
        h_Energy_All[i]->SetLineColor(1);
        h_Energy_select[i]->Draw("same");
        h_Energy_select[i]->SetLineColor(2);
        h_Energy_All[i]->SetTitle(Form("%s",chanName.c_str()));
        h_Energy_All[i]->GetXaxis()->SetTitle("Energy [keVee]");
        h_Energy_All[i]->GetYaxis()->SetTitle("Number of counts");
        ccEnergy->cd(i+1)->SetLogy();
        ccEnergy->cd(i+1)->SetGrid();
     //   ccEnergy->cd(i+1)->Update();
    }
   // ccEnergy->Update();
    ccEnergy->SaveAs(Form("%s/%s/EnergyDistribution_%s.png",path_selectedtemplates.c_str(),DataListName.c_str(),DataListName.c_str()));

    if(N_Ion_Channel>0)
    {
    TCanvas * ccCharge = new TCanvas();
    ccCharge->Divide(2,2);
    ccCharge->cd(1);
    gr_IonAVsIonB->Draw("ap");
    gr_IonAVsIonB->SetMarkerStyle(6);
    gr_IonAVsIonB_select->Draw("samep");
    gr_IonAVsIonB_select->SetMarkerStyle(6);
    gr_IonAVsIonB_select->SetMarkerColor(2);
    gr_IonAVsIonB->GetXaxis()->SetRangeUser(-5,14);
    gr_IonAVsIonB->GetYaxis()->SetRangeUser(-5,14);
    gr_IonAVsIonB->GetXaxis()->SetTitle("Decor Ion A [keVee]");
    gr_IonAVsIonB->GetYaxis()->SetTitle("Decor Ion B [keVee]");
    //ccCharge->cd(1)->Update();
    
    ccCharge->cd(2);
    gr_IonCVsIonD->Draw("ap");
    gr_IonCVsIonD->SetMarkerStyle(6);
    gr_IonCVsIonD_select->Draw("samep");
    gr_IonCVsIonD_select->SetMarkerStyle(6);
    gr_IonCVsIonD_select->SetMarkerColor(2);
    gr_IonCVsIonD->GetXaxis()->SetRangeUser(-5,14);
    gr_IonCVsIonD->GetYaxis()->SetRangeUser(-5,14);
    gr_IonCVsIonD->GetXaxis()->SetTitle("Decor Ion C [keVee]");
    gr_IonCVsIonD->GetYaxis()->SetTitle("Decor Ion D [keVee]");
    //ccCharge->cd(2)->Update();

    ccCharge->cd(3);
    gr_IonBVsIonD->Draw("ap");
    gr_IonBVsIonD->SetMarkerStyle(6);
    gr_IonBVsIonD_select->Draw("samep");
    gr_IonBVsIonD_select->SetMarkerStyle(6);
    gr_IonBVsIonD_select->SetMarkerColor(2);
    gr_IonBVsIonD->GetXaxis()->SetRangeUser(-5,14);
    gr_IonBVsIonD->GetYaxis()->SetRangeUser(-5,14);
    gr_IonBVsIonD->GetXaxis()->SetTitle("Decor Ion B [keVee]");
    gr_IonBVsIonD->GetYaxis()->SetTitle("Decor Ion D [keVee]");
    //ccCharge->cd(3)->Update();

    ccCharge->cd(4);
    gr_BottomVsTop->Draw("ap");
    gr_BottomVsTop->SetMarkerStyle(6);
    gr_BottomVsTop_select->Draw("samep");
    gr_BottomVsTop_select->SetMarkerStyle(6);
    gr_BottomVsTop_select->SetMarkerColor(2);
    gr_BottomVsTop->GetXaxis()->SetRangeUser(-5,14);
    gr_BottomVsTop->GetYaxis()->SetRangeUser(-5,14);
    gr_BottomVsTop->GetXaxis()->SetTitle("Decor A+B [keVee]");
    gr_BottomVsTop->GetYaxis()->SetTitle("Decor C+D [keVee]");
    //ccCharge->cd(4)->Update();
    //ccCharge->Update();
    ccCharge->SaveAs(Form("%s/%s/ChargeSelection_%s.png",path_selectedtemplates.c_str(),DataListName.c_str(),DataListName.c_str()));

    
    TCanvas * ccHeatVsIon = new TCanvas();
    gr_HeatVsBD->Draw("ap");
    gr_HeatVsBD->SetMarkerStyle(6);
    gr_HeatVsBD_select->Draw("samep");
    gr_HeatVsBD_select->SetMarkerStyle(6);
    gr_HeatVsBD_select->SetMarkerColor(2);
  //  fun_EionVsEheat_inf->Draw("same");
  //  fun_EionVsEheat_sup->Draw("same");
    fun_EionVsEheat_inf->SetLineColor(4);
    fun_EionVsEheat_sup->SetLineColor(4);
    gr_HeatVsBD->GetXaxis()->SetRangeUser(-2,15);
    gr_HeatVsBD->GetYaxis()->SetRangeUser(-4,30);
    gr_HeatVsBD->GetXaxis()->SetTitle("Heat [keVee]");
    gr_HeatVsBD->GetYaxis()->SetTitle(" EiFID+EiVETO [keVee]");
    //ccHeatVsIon->Update();
    ccHeatVsIon->SaveAs(Form("%s/%s/HeatVsIonFIDSelection_%s.png",path_selectedtemplates.c_str(),DataListName.c_str(),DataListName.c_str()));
    
    
    
    gr_HeatVsBD_select_goodioni ->SetMarkerColor(kRed+1);
    gr_HeatVsBD_select_partioni ->SetMarkerColor(kBlue+1);
    gr_HeatVsBD_select_ionitail ->SetMarkerColor(kGreen +1);

    gr_HeatVsBD_select_goodioni ->SetMarkerStyle(6);
    gr_HeatVsBD_select_partioni ->SetMarkerStyle(6);
    gr_HeatVsBD_select_ionitail ->SetMarkerStyle(6);


   gr_HeatVsBD_select_goodioni ->Draw("p");
   gr_HeatVsBD_select_partioni ->Draw("samep");
   gr_HeatVsBD_select_ionitail ->Draw("samep");
   //ccHeatVsIon->Update();
   ccHeatVsIon->SaveAs(Form("%s/%s/HeatVsIonFIDSelection_population_%s.png",path_selectedtemplates.c_str(),DataListName.c_str(),DataListName.c_str()));
   
   
   TCanvas * ccHeatVsHeat = new TCanvas();
   TH2D * axes = new TH2D("axe","axe",1000,0.,15.,1000,0.,15.);
   axes->GetXaxis()->SetTitle(" Heat A (keV)");
	axes->GetXaxis()->SetTitleOffset(1.25);
	axes->GetYaxis()->SetTitle("Heat B (keV)");
	axes->GetYaxis()->SetTitleOffset(1.3);
	axes->SetStats(kFALSE);
	axes->Draw();
   gr_HeatAVsHeatB ->SetMarkerStyle(6);
   gr_HeatAVsHeatB->Draw("PSAME");
   gr_HeatAVsHeatB_select ->SetMarkerStyle(6);
   gr_HeatAVsHeatB_select ->SetMarkerColor(2);
   gr_HeatAVsHeatB_select->Draw("PSAME");
   ccHeatVsHeat->SaveAs(Form("%s/%s/HeatAVsHeatB_%s.png",path_selectedtemplates.c_str(),DataListName.c_str(),DataListName.c_str()));
    }


    fileOut.cd();
    tree->Write();
    fileOut.Close();
    
    
    return;
}

int GetRootFileNames(string DataList)
{
    string RQfilename = Form("%s.list",DataList.c_str());
    int NumofFiles;
    string anRQ;
    Filelist.clear();
    ifstream infile (RQfilename.c_str(), ios_base::in);
    if(infile.fail())
    {
        cout <<"Error reading data file List !" << endl;
        exit(1);
    }
    
    while (getline(infile, anRQ, '\n'))
        Filelist.push_back( anRQ );
    NumofFiles = Filelist.size();
    cout<<"Num of data files considered = "<<NumofFiles<<endl;
    
    return NumofFiles;
}

double GetHeatCalibration(double GainA, double GainB, double *NonLinearityA, double *NonLinearityB, double FiducialVoltage, double *Amplitude_OF, string HeatChannel)
{
    double Amplitude, Energy;
    double Gain;
    double NonLinearity[5];
    double HeatLinear, correctedHeat, rlog, arg, tmp, gcor;
    
    if(HeatChannel == "HeatA")
    {
        Amplitude = Amplitude_OF[0];
        Gain = GainA;
        for(int i = 0; i<5; i++){NonLinearity[i] = NonLinearityA[i];}
    }
    if(HeatChannel == "HeatB")
    {
        Amplitude = Amplitude_OF[1];
        Gain = GainB;
        for(int i = 0; i<5; i++){NonLinearity[i] = NonLinearityB[i];}
    }

    HeatLinear = Amplitude * Gain;
    if(HeatLinear>0 && NonLinearity[1] != 0)
    {
        rlog = log( HeatLinear * (3.+fabs(FiducialVoltage))/11.);
        arg = (NonLinearity[2] - rlog)/NonLinearity[3];
        tmp = NonLinearity[1] * (rlog - NonLinearity[2]) *  (1. - NonLinearity[4] * (rlog - NonLinearity[2]));
        gcor = NonLinearity[0] + tmp / ( 1. + exp(arg));
        correctedHeat = HeatLinear / gcor;
        
    } else {
        
        correctedHeat = HeatLinear / NonLinearity[0];
    }
    
    Energy = correctedHeat;
    
    
    return Energy;
}

void GetIonCalibration(double GainIonA, double GainIonB, double GainIonC, double GainIonD, double CrossTalkAB, double CrossTalkBA, double CrossTalkAC, double CrossTalkCA, double CrossTalkAD, double CrossTalkDA, double CrossTalkBC, double CrossTalkCB, double CrossTalkBD, double CrossTalkDB, double CrossTalkCD, double CrossTalkDC, double *Amplitude_OF, double *Energy_Ion)
{
    TMatrixD Xtalk(4,4);
    Xtalk(0,0) = -1;
    Xtalk(0,1) = CrossTalkAB;
    Xtalk(0,2) = CrossTalkAC;
    Xtalk(0,3) = CrossTalkAD;
    Xtalk(1,0) = CrossTalkBA;
    Xtalk(1,1) = -1;
    Xtalk(1,2) = CrossTalkBC;
    Xtalk(1,3) = CrossTalkBD;
    Xtalk(2,0) = CrossTalkCA;
    Xtalk(2,1) = CrossTalkCB;
    Xtalk(2,2) = -1;
    Xtalk(2,3) = CrossTalkCD;
    Xtalk(3,0) = CrossTalkDA;
    Xtalk(3,1) = CrossTalkDB;
    Xtalk(3,2) = CrossTalkDC;
    Xtalk(3,3) = -1;
    
    double Amplitude_Xtalk[4];
    double Amplitude_Raw[4];
    double Gain_vec[4];
    Gain_vec[0] = GainIonA;
    Gain_vec[1] = GainIonB;
    Gain_vec[2] = GainIonC;
    Gain_vec[3] = GainIonD;
    
    for(int i = 0; i<4;i++)
    {
        Amplitude_Xtalk[i] = 0;
        Amplitude_Raw[i] = Amplitude_OF[2+i];
        Energy_Ion[i] = 0;
    }
    
    for(int i = 0; i<4; i++)
    {
        for(int j = 0; j<4; j++)
        {
            Amplitude_Xtalk[i] = Amplitude_Xtalk[i] + Xtalk(i,j)*Amplitude_Raw[j];
        }
    }
    
    for(int i = 0; i<4; i++){Energy_Ion[i] = -1. * Gain_vec[i] * Amplitude_Xtalk[i];}
    
    
    return;
}

