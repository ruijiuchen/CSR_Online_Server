#include <stdlib.h>
#include <stdio.h>
#include <sstream> 
#include <iostream> 
#include <sys/ioctl.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <signal.h>
#include <math.h>

#include <assert.h>
#include <TApplication.h>
#include <TKey.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TSocket.h>
#include <TServerSocket.h>
#include <TMessage.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TFolder.h>
#include <TRint.h>
#include <TCutG.h>
#include <TThread.h>
#include <TChain.h>
#include <set> // Include the set header
#define POINTER_T long int
#define INT int
#define BOOL bool
#define TRUE true
#define FALSE false
#define SUCCESS 1
#define CM_SUCCESS 1
#include <fstream>
#include <algorithm> // Include the algorithm header
#include <string>

/* Our own ROOT global objects */
TApplication *manaApp;
TFolder *gManaHistosFolder               = NULL;      // Container for all histograms
TFolder *histo_folder                    = NULL; // histo folder
TFile *gManaOutputFile                   = NULL;  // MIDAS output file
TObjArray *gHistoFolderStack             = NULL;    //
int InjectionRange                       = 300;
int InjectionRange_Last16Hours           = NULL;
TFile *f = NULL;
TH2F *FFT = NULL;
TH1D *projectionX = NULL;
// Define a struct to store your data
struct MyData_iq {
  double time;
  double FileID;
  double FrameID;
  double I_avr;
  double Q_avr;
  double IQ_avr;
  double Area;
  double Width;
  double Mean;
  
  std::string file_path;
};
struct MyData_sc {
  double time;
  double sc_ch1;
  double sc_ch2;
  double sc_ch3;
  double sc_ch4;
  double sc_ch5;
  double sc_ch6;
  double sc_ch7;
  double sc_ch8;
};
std::vector<MyData_iq> data_iq; // Vector to store the accumulated data
std::vector<MyData_iq> data_iq_injection; // Vector to store the accumulated data
TH1F *hI_avr_time = NULL;
TH1F *hQ_avr_time = NULL;
TH1F *hIQ_avr_time = NULL;
TH1F *hArea_time = NULL;
TH1F *hWidth_time = NULL;
TH1F *hMean_time = NULL;
std::vector<MyData_sc> data_sc; // Vector to store the accumulated data
TH1F *hsc_ch1_time = NULL;
TH1F *hsc_ch2_time = NULL;
TH1F *hsc_ch3_time = NULL;
TH1F *hsc_ch4_time = NULL;
TH1F *hsc_ch5_time = NULL;
TH1F *hsc_ch6_time = NULL;
TH1F *hsc_ch7_time = NULL;
TH1F *hsc_ch8_time = NULL;
TH1F *hDCCT_time = NULL;
//TFile *f;
//TH2F *FFT;
TH1D *hAccumulation = NULL;
int frameID_offset_number=50;
std::vector<TH1D *> hAccumulation_vector;
std::vector<std::string> analyzed_files_sc;
std::vector<std::string> analyzed_files_iq;
//double range_sc = nx_sc;
//int nx_sc = 3600*24;
double range_sc = 262144*1./20e6*4096*1000;
int nx_sc = 4096*1000;

double time_offset = 0.060; //s
double range_iq = 262144*1./20e6*4096*1000;
int nx_iq = 4096*1000;
double frameID_offset=-10;
double frameID_range=1;

int injection = 0;
//double time = 0,sc_ch5=0; 
/* Global mutex */
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
/*book functions put a root object in a suitable folder
 *for histos, there are a lot of types, so we use templates.
 *for other objects we have one function per object
 */
int sum, Injection,Ions_0=0,IonsZoomed_0=0,Ions_1=0,IonsZoomed_1=0,IonsWithin;
double REVTIME_EffectiveIons_Min=660,REVTIME_EffectiveIons_Max=680,REVTIME_101In_Min=669.8,REVTIME_101In_Max=670.6;
double Hour_Range=16,Hour_Range1=2;
TDatime Time0(1995,1,1,8,0,0);
int Time_X0 = Time0.Convert(),IonNumberPerInjection=0,EffectiveIonNumberPerInjection=0;

using namespace std;
template < typename TH1X >
TH1X * h1_book(const char *name, const char *title,
	       int bins, double min, double max)
{
  TH1X *hist;
  /* check if histo already exists */
  if (!gHistoFolderStack->Last())
    hist = (TH1X *) gManaHistosFolder->FindObjectAny(name);
  else
    hist = (TH1X *) ((TFolder *) gHistoFolderStack->Last())->FindObjectAny(name);
  
  if (hist == NULL) {
    hist = new TH1X(name, title, bins, min, max);
    if (!gHistoFolderStack->Last())
      gManaHistosFolder->Add(hist);
    else
      ((TFolder *) gHistoFolderStack->Last())->Add(hist);
  }
  return hist;
}
template < typename TGraph >TGraph * g1_book(const char *name, const char *title)
{
  TGraph *graph;
  /* check if histo already exists */
  if (!gHistoFolderStack->Last())
    graph = (TGraph *) gManaHistosFolder->FindObjectAny(name);
  else
    graph = (TGraph *) ((TFolder *) gHistoFolderStack->Last())->FindObjectAny(name);
  if (graph == NULL) 
    {
      graph = new TGraph();
      graph->SetName(name);
      graph->SetTitle(title);
      if (!gHistoFolderStack->Last())
	gManaHistosFolder->Add(graph);
      else
	((TFolder *) gHistoFolderStack->Last())->Add(graph);
    }
  return graph;
}

bool ss_kbhit()
/********************************************************************\

  Routine: ss_kbhit

  Purpose: Returns TRUE if a key is pressed

  Input:
    none

  Output:
    none

  Function value:
    FALSE                 No key has been pressed
    TRUE                  Key has been pressed

\********************************************************************/
{
  int n;
  ioctl(0, FIONREAD, &n);
  return (n > 0);
}
// Define a custom comparison function for sorting based on time
bool compareByTime_iq(const MyData_iq& a, const MyData_iq& b) {
  return a.time < b.time;
}
void load_data_from_files_list_sc(string file_list_path, std::vector<std::string>& analyzed_files, std::vector<MyData_sc>& data)
{
  // Create a set to store unique file paths
  std::set<std::string> unique_files;
  for (const std::string& file_path : analyzed_files) {
    unique_files.insert(file_path); 
  }
  // Read file paths from "analyzed_files_sc.txt"
  std::ifstream file_list(file_list_path);
  if (file_list.is_open()) {
    std::string line,file_path;
    while (std::getline(file_list, line)) {
      size_t lastSpace = line.find_last_of(" \t");
      if (lastSpace != std::string::npos) {
	std::string file_path = line.substr(lastSpace + 1); // Extract the last part
	// Check if the file path is not a duplicate before adding it
	if (unique_files.find(file_path) == unique_files.end()) {
	  analyzed_files.push_back(file_path);
	  unique_files.insert(file_path); // Add to the set to mark as seen
	  cout<<file_path<<endl;
	  // read data from the file_path to data.
	  TFile *f = new TFile(file_path.c_str());
	  TTree *t;
	  f->GetObject("t",t);
	  double time = 0;
	  double sc_ch1=0;
	  double sc_ch2=0;
	  double sc_ch3=0;
	  double sc_ch4=0;
	  double sc_ch5=0;
	  double sc_ch6=0;
	  double sc_ch7=0;
	  double sc_ch8=0;
	  
	  t->SetBranchAddress("time",   &time);
	  t->SetBranchAddress("sc_ch1", &sc_ch1);
	  t->SetBranchAddress("sc_ch2", &sc_ch2);
	  t->SetBranchAddress("sc_ch3", &sc_ch3);
	  t->SetBranchAddress("sc_ch4", &sc_ch4);
	  t->SetBranchAddress("sc_ch5", &sc_ch5);
	  t->SetBranchAddress("sc_ch6", &sc_ch6);
	  t->SetBranchAddress("sc_ch7", &sc_ch7);
	  t->SetBranchAddress("sc_ch8", &sc_ch8);
	  int NEntries = t->GetEntries();
	  for(int i=0; i < NEntries; i++)
	    {
	      t->GetEntry(i);
	      // Store the data in the vector
	      MyData_sc entry;
	      entry.time = time;
	      entry.sc_ch1 = sc_ch1;
	      entry.sc_ch2 = sc_ch2;
	      entry.sc_ch3 = sc_ch3;
	      entry.sc_ch4 = sc_ch4;
	      entry.sc_ch5 = sc_ch5;
	      entry.sc_ch6 = sc_ch6;
	      entry.sc_ch7 = sc_ch7;
	      entry.sc_ch8 = sc_ch8;
	      data.push_back(entry);
	    }
	  f->Close();
	}
      }
    }
    file_list.close();
  }
  else{
      std::cerr << "Unable to open file: analyzed_files_sc.txt" << std::endl;
    }  
}
void load_data_from_files_list_iq(string file_list_path, std::vector<std::string>& analyzed_files, std::vector<MyData_iq>& data)
{
  // Create a set to store unique file paths
  std::set<std::string> unique_files;
  for (const std::string& file_path : analyzed_files) {
    unique_files.insert(file_path); 
  }
  // Read file paths from "analyzed_files_iq.txt"
  std::ifstream file_list(file_list_path);
  if (file_list.is_open()) {
    std::string line,file_path;
    while (std::getline(file_list, line)) {
      size_t lastSpace = line.find_last_of(" \t");
      if (lastSpace != std::string::npos) {
	std::string file_path = line.substr(lastSpace + 1); // Extract the last part
	// Check if the file path is not a duplicate before adding it
	if (unique_files.find(file_path) == unique_files.end()) {
	  analyzed_files.push_back(file_path);
	  unique_files.insert(file_path); // Add to the set to mark as seen
	  //cout<<file_path<<endl;
	  //// read data from the file_path to data.
	  TFile *f = new TFile(file_path.c_str());
	  TTree *t;
	  TH2F *FFT = NULL;
	  f->GetObject("t_avr",t);
	  double time = 0,FileID,FrameID,I_avr,Q_avr,IQ_avr,Area=0,Width=0,Mean=0,FreRange[10][2];
	  int NFreRange[10][2];
	  t->SetBranchAddress("time",   &time);
	  t->SetBranchAddress("FileID", &FileID);
	  t->SetBranchAddress("FrameID",&FrameID);
	  t->SetBranchAddress("I_avr",  &I_avr);
	  t->SetBranchAddress("Q_avr",  &Q_avr);
	  t->SetBranchAddress("IQ_avr", &IQ_avr);
	  int NEntries = t->GetEntries();
	  f->GetObject("FFT", FFT);
	  FreRange[0][0] = -1910; // min of freqnency range
	  FreRange[0][1] = -1880; // max 
	  for(int i=0; i < NEntries; i++)
	    {
	      t->GetEntry(i);
	      // Store the data in the vector
	      MyData_iq entry;
	      entry.time = time;
	      entry.FileID = FileID;
	      entry.FrameID = FrameID;
	      entry.I_avr = I_avr;
	      entry.Q_avr = Q_avr;
	      entry.IQ_avr = IQ_avr;
	      entry.file_path=file_path;
	      
	      TH1D *projectionX = FFT->ProjectionX(Form("ProjectionX_File%d_Frame%d", FileID, FrameID), FrameID, FrameID);
	      //projectionX->GetXaxis()->SetRangeUser(FreRange[0][1], FreRange[0][1]);
	      
	      Mean  = projectionX->GetXaxis()->GetBinCenter(projectionX->GetMaximumBin());
              Width = projectionX->GetRMS();
	      NFreRange[0][0] = projectionX->GetXaxis()->FindBin(Mean-1);
              NFreRange[0][1] = projectionX->GetXaxis()->FindBin(Mean+1);
              //for(int k=0;k<2;k++)  NFreRange[0][k] = projectionX->GetXaxis()->FindBin(FreRange[0][k]);
	      Area  = projectionX->Integral(NFreRange[0][0],NFreRange[0][1]);
	      //cout<<" FileID = "<<FileID<<" FrameID = "<<FrameID<<" Mean = "<<Mean<<endl;
	      entry.Area  = Area;
	      entry.Width = Width;
	      entry.Mean  = Mean;
	      delete projectionX;
	      data.push_back(entry);
	    }
	  //////Calculate area, width and mean values of peak.
	  
	  //////
	  f->Close();
	}
      }
    }
    file_list.close();
  }
  else{
      std::cerr << "Unable to open file: analyzed_files_sc.txt" << std::endl;
    }  
}

std::vector<MyData_iq>  find_injection(
				       std::vector<MyData_sc>& data_sc,
				       std::vector<MyData_iq>& data_iq,
				       std::vector<MyData_iq>& data_iq_injection,
				       double time_offset
				       )
{
  cout<<"chenrj ... find_injection 1"<<endl;
  std::vector<MyData_iq> data_iq_injection_new;
  for (const auto& entry_sc : data_sc) {
    if (entry_sc.sc_ch4 == 1) {
      // Add a condition to check if entry_sc.time is equal to data_iq_injection.time - time_offset
      bool skipEntry = false;
      //cout<<" find_injection "<< injection<<endl;
      for (const auto& entry_iq_inj : data_iq_injection) {
	//if (entry_sc.time == entry_iq_inj.time) {
	if (fabs(entry_sc.time - (entry_iq_inj.time - time_offset))<0.1) {
	  skipEntry = true;
	  break; // Exit the loop as we don't need to process this entry_sc
	}
      }
      if (skipEntry) {
	continue; // Skip the rest of the loop for this entry_sc
      }
      cout<<"chenrj ... find_injection 2"<<endl;
      double closest_time_difference = std::numeric_limits<double>::max();
      MyData_iq* closest_entry_iq = nullptr;
      for (auto& entry_iq : data_iq) {
	double time_difference = fabs(entry_sc.time - (entry_iq.time - time_offset));
	//double time_difference = fabs(entry_sc.time - entry_iq.time);
	if (time_difference < closest_time_difference && time_difference < 0.01) {
	  closest_time_difference = time_difference;
	  closest_entry_iq = &entry_iq;
	}
      }
      cout<<"chenrj ... find_injection 3"<<endl;
      if (closest_entry_iq) {
	//cout<<" chenrj2 "<< std::fixed<<entry_sc.time <<"  "<<closest_entry_iq->time - time_offset<<" closest_time_difference "<<closest_time_difference<<" time_offset "<<time_offset<<"FileID = "<<closest_entry_iq->FileID<<" FrameID = "<<closest_entry_iq->FrameID<<endl;
	//cout<<" chenrj2 "<< std::fixed<<entry_sc.time <<"  "<<entry_iq.time<<" "<<entry_sc.time - (entry_iq.time)<<" time_offset "<<time_offset<<"FileID = "<<entry_iq.FileID<<" FrameID = "<<entry_iq.FrameID<<endl;
	
	data_iq_injection.    push_back(*closest_entry_iq);
	data_iq_injection_new.push_back(*closest_entry_iq);
      }
    }
  }
  cout<<"chenrj ... find_injection 4"<<endl;
  return data_iq_injection_new;
}

void create_accumulated_spectrum(
				 std::vector<MyData_iq>& data_iq_injection,
				 double frameID_offset,
				 double frameID_range
				 ) {
  cout<<"chenrj ... create_accumulated_spectrum 1"<<endl;
  std::string previousFilePath;
  //TFile *f = NULL;
  //TH2F *FFT = NULL;
  double xMin=-2000;
  double xMax= 2000;
  double IQSamplingRate = 20e6; //Samples/s
  double IQRecordPerFrame=pow(2,18);
  double hstep = IQSamplingRate/IQRecordPerFrame/1000;
  int nBinsX = (xMax - xMin)/hstep;
  cout<<"chenrj ... create_accumulated_spectrum 1"<<endl;
  for (const auto& entry_iq_injection : data_iq_injection) {
    if (entry_iq_injection.file_path != previousFilePath) {
      previousFilePath = entry_iq_injection.file_path;
      cout << "injection = " << injection << " " << entry_iq_injection.file_path << " FrameID " << entry_iq_injection.FrameID << endl;
      if (f) {
	f->Close();
	f = nullptr;
      }
      f = new TFile(entry_iq_injection.file_path.c_str());
      f->GetObject("FFT", FFT);
    }
    cout<<"chenrj ... create_accumulated_spectrum 2"<<endl;
    if (FFT) {
      for(int i=0;i<frameID_offset_number;i++){
	frameID_offset = i;
	int fileID = static_cast<int>(entry_iq_injection.FileID);
	int frameID = static_cast<int>(entry_iq_injection.FrameID);
	int frameID_min = frameID - frameID_offset;
	int frameID_max = frameID - frameID_offset + frameID_range;
	projectionX = FFT->ProjectionX(Form("ProjectionX_%d_File%d_Frame%d_%d", injection, fileID, frameID_min, frameID_max), frameID_min, frameID_max);
	
	if (injection == 0) {
	  nBinsX = projectionX->GetNbinsX();
	  xMin = projectionX->GetXaxis()->GetXmin();
	  xMax = projectionX->GetXaxis()->GetXmax();
	  hAccumulation_vector[i]->SetBins(nBinsX, xMin, xMax);
	}
	
	for (int binX = 1; binX <= nBinsX; ++binX) {
	  double binContent = projectionX->GetBinContent(binX);
	  double binCenterX = projectionX->GetXaxis()->GetBinCenter(binX);
	  int binX_accumulation = hAccumulation_vector[i]->GetXaxis()->FindBin(binCenterX);
	  double accumulationValue = hAccumulation_vector[i]->GetBinContent(binX_accumulation);
	  hAccumulation_vector[i]->SetBinContent(binX_accumulation, accumulationValue + binContent);
	}
	delete projectionX;
      }
    }
    injection = injection + 1;
  }
  cout<<"chenrj ... create_accumulated_spectrum 3"<<endl;
  //f->Close();
  cout<<"chenrj ... create_accumulated_spectrum 4"<<endl;
}


/* Create histograms */
void create_hist()
{
  double xmax_sc = 1000;
  // scaler data
  load_data_from_files_list_sc("/lustre/astrum/rchen/OnlineDataAnalysisSystem/analyze_iq_sc/analyzed_files_sc.txt",analyzed_files_sc, data_sc);
  load_data_from_files_list_iq("/lustre/astrum/rchen/OnlineDataAnalysisSystem/analyze_iq_sc/analyzed_files_iq.txt",analyzed_files_iq, data_iq);
  // Check if data_sc is not empty before updating xmax_sc
  if (!data_sc.empty()) {
    xmax_sc = data_sc.back().time;
  }
  double xmin_sc = xmax_sc - range_sc;
  hsc_ch1_time    = h1_book<TH1F>("hsc_ch1_time", "sc_ch1 vs time", nx_sc, xmin_sc, xmax_sc);
  hsc_ch2_time    = h1_book<TH1F>("hsc_ch2_time", "sc_ch2 vs time", nx_sc, xmin_sc, xmax_sc);
  hsc_ch3_time    = h1_book<TH1F>("hsc_ch3_time", "sc_ch3 vs time", nx_sc, xmin_sc, xmax_sc);
  hsc_ch4_time    = h1_book<TH1F>("hsc_ch4_time", "sc_ch4 vs time", nx_sc, xmin_sc, xmax_sc);
  hsc_ch5_time    = h1_book<TH1F>("hsc_ch5_time", "sc_ch5 vs time", nx_sc, xmin_sc, xmax_sc);
  hsc_ch6_time    = h1_book<TH1F>("hsc_ch6_time", "sc_ch6 vs time", nx_sc, xmin_sc, xmax_sc);
  hsc_ch7_time    = h1_book<TH1F>("hsc_ch7_time", "sc_ch7 vs time", nx_sc, xmin_sc, xmax_sc);
  hsc_ch8_time    = h1_book<TH1F>("hsc_ch8_time", "sc_ch8 vs time", nx_sc, xmin_sc, xmax_sc);
  hDCCT_time      = h1_book<TH1F>("hDCCT_time",   "DCCT vs time",   nx_sc, xmin_sc, xmax_sc);

  
  hsc_ch1_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
  hsc_ch2_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
  hsc_ch3_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
  hsc_ch4_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
  hsc_ch5_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
  hsc_ch6_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
  hsc_ch7_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
  hsc_ch8_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
  hDCCT_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
  
  //Fill data from data_sc
  for (const auto& entry : data_sc) {
    hsc_ch1_time->Fill(entry.time, entry.sc_ch1);
    hsc_ch2_time->Fill(entry.time, entry.sc_ch2);
    hsc_ch3_time->Fill(entry.time, entry.sc_ch3);
    hsc_ch4_time->Fill(entry.time, entry.sc_ch4*-1899.5);
    hsc_ch5_time->Fill(entry.time, entry.sc_ch5);
    hsc_ch6_time->Fill(entry.time, entry.sc_ch6);
    hsc_ch7_time->Fill(entry.time, entry.sc_ch7);
    hsc_ch8_time->Fill(entry.time, entry.sc_ch8);
    hDCCT_time  ->Fill(entry.time, entry.sc_ch5);
  }
  
  // iq
  double min_time_iq = 0;
  if (!data_iq.empty()) {
    // Find the minimum time value in data_iq
    min_time_iq = data_iq[0].time;
    for (const auto& entry_iq : data_iq) {
      if (entry_iq.time < min_time_iq) {
	min_time_iq = entry_iq.time;
      }
    }
    // Calculate time_offset
    //time_offset = min_time_iq - data_sc[0].time;
    cout << "time_offset " << time_offset << endl;
  }
  double xmax_iq = 1000;
  if (!data_iq.empty()) {
    xmax_iq = data_iq.back().time;
  }
  double xmin_iq = xmax_iq - range_iq;
  hI_avr_time     = h1_book<TH1F>("hI_avr_time",  "I_avr vs time",  nx_iq, xmin_iq, xmax_iq);
  hQ_avr_time     = h1_book<TH1F>("hQ_avr_time",  "Q_avr vs time",  nx_iq, xmin_iq, xmax_iq);
  hIQ_avr_time    = h1_book<TH1F>("hIQ_avr_time", "IQ_avr vs time", nx_iq, xmin_iq, xmax_iq);
  hArea_time      = h1_book<TH1F>("hArea_time",   "Area vs time",   nx_iq, xmin_iq, xmax_iq);
  hWidth_time     = h1_book<TH1F>("hWidth_time",  "Width vs time",  nx_iq, xmin_iq, xmax_iq);
  hMean_time      = h1_book<TH1F>("hMean_time",   "Mean vs time",   nx_iq, xmin_iq, xmax_iq);
  
  // Fill data from data_iq
  for (const auto& entry : data_iq) {
    hI_avr_time ->Fill(entry.time - time_offset, entry.I_avr);
    hQ_avr_time ->Fill(entry.time - time_offset, entry.Q_avr);
    hIQ_avr_time->Fill(entry.time - time_offset, entry.IQ_avr);
    hArea_time  ->Fill(entry.time - time_offset, entry.Area);
    hWidth_time ->Fill(entry.time - time_offset, entry.Width);
    hMean_time  ->Fill(entry.time - time_offset, entry.Mean);
  }

  hAccumulation = h1_book<TH1D>("hAccumulation", "hAccumulation", 1000, -1000, 1000);
  for (int i=0;i<frameID_offset_number;i++)
    {
      int frameID_offset=i;
      int frameID_range =1;
      TH1D *h = h1_book<TH1D>(Form("hAccumulation_offset_%d_range_%d",frameID_offset,frameID_range), Form("accumulated spectrum after injection, offset=%d/frame, peojection range=%d/frame.",frameID_offset,frameID_range), 1000, -1000, 1000);
      hAccumulation_vector.push_back(h);
    }
  std::vector<MyData_iq> data_iq_injection_new =  find_injection(data_sc, data_iq, data_iq_injection, time_offset);
  create_accumulated_spectrum(data_iq_injection_new, frameID_offset, frameID_range);
  TFile *fout = new TFile("result.root","recreate");
  hsc_ch1_time->Write();
  hsc_ch2_time->Write();
  hsc_ch3_time->Write();
  hsc_ch4_time->Write();
  hsc_ch5_time->Write();
  hsc_ch6_time->Write();
  hsc_ch7_time->Write();
  hsc_ch8_time->Write();
  hDCCT_time  ->Write();
  hI_avr_time ->Write();
  hQ_avr_time ->Write();
  hIQ_avr_time->Write();
  hArea_time  ->Write();
  hWidth_time ->Write();
  hMean_time  ->Write();
  for (int i=0;i<frameID_offset_number;i++)
    hAccumulation_vector[i]->Write();

  fout->Close();
}
/* update histograms */
void update_hists()
{
  double Start  = clock();
  cout<<"update_hists"<<endl;  
  /* To be implemented */
  load_data_from_files_list_sc("/lustre/astrum/rchen/OnlineDataAnalysisSystem/analyze_iq_sc/analyzed_files_sc.txt",analyzed_files_sc,data_sc);
  load_data_from_files_list_iq("/lustre/astrum/rchen/OnlineDataAnalysisSystem/analyze_iq_sc/analyzed_files_iq.txt",analyzed_files_iq,data_iq);
  cout<<"chenrj 1..."<<endl;
  // sc
  // Check if data_sc is not empty before updating xmax
  double xmax_sc = 1000;
  if (!data_sc.empty()) {
    int i = data_sc.size() - 1;
    xmax_sc =  data_sc[i].time;
    double xmin_sc = xmax_sc  - range_sc;
    //cout<<"xmax_sc = "<<xmax_sc<<endl;
    hsc_ch1_time->Reset();
    hsc_ch2_time->Reset();
    hsc_ch3_time->Reset();
    hsc_ch4_time->Reset();
    hsc_ch5_time->Reset();
    hsc_ch6_time->Reset();
    hsc_ch7_time->Reset();
    hsc_ch8_time->Reset();
    hDCCT_time  ->Reset();
    hsc_ch1_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
    hsc_ch2_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
    hsc_ch3_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
    hsc_ch4_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
    hsc_ch5_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
    hsc_ch6_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
    hsc_ch7_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
    hsc_ch8_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
    hDCCT_time->GetXaxis()->SetLimits(xmin_sc, xmax_sc);
    // Fill data from data_sc into hDCCT_time starting from the last element
    for (int i = data_sc.size() - 1; i >= 0; i--) {
      const auto& entry = data_sc[i];
      if (entry.time < xmin_sc) {
	// Stop filling when entry.time is less than xmax_sc - nx_sc
	break;
      }
      hsc_ch1_time->Fill(entry.time, entry.sc_ch1);
      hsc_ch2_time->Fill(entry.time, entry.sc_ch2);
      hsc_ch3_time->Fill(entry.time, entry.sc_ch3);
      hsc_ch4_time->Fill(entry.time, entry.sc_ch4);
      hsc_ch5_time->Fill(entry.time, entry.sc_ch5);
      hsc_ch6_time->Fill(entry.time, entry.sc_ch6);
      hsc_ch7_time->Fill(entry.time, entry.sc_ch7);
      hsc_ch8_time->Fill(entry.time, entry.sc_ch8);
      hDCCT_time  ->Fill(entry.time, entry.sc_ch5);
    }
  }
  //
  cout<<"chenrj 2..."<<endl;
  // iq
  // Check if data_iq is not empty before updating xmax
  double xmax_iq = 1000;
  if (!data_iq.empty()) {
    int i = data_iq.size() - 1;
    xmax_iq =  data_iq[i].time;
    double xmin_iq = xmax_iq  - range_iq;
    //cout<<"xmax_iq = "<<xmax_iq<<endl;
    hI_avr_time  ->Reset();
    hQ_avr_time  ->Reset();
    hIQ_avr_time ->Reset();
    hArea_time   ->Reset();
    hWidth_time  ->Reset();
    hMean_time   ->Reset();
    hI_avr_time  ->GetXaxis()->SetLimits(xmin_iq, xmax_iq);
    hQ_avr_time  ->GetXaxis()->SetLimits(xmin_iq, xmax_iq);
    hIQ_avr_time ->GetXaxis()->SetLimits(xmin_iq, xmax_iq);
    hArea_time   ->GetXaxis()->SetLimits(xmin_iq, xmax_iq);
    hWidth_time  ->GetXaxis()->SetLimits(xmin_iq, xmax_iq);
    hMean_time   ->GetXaxis()->SetLimits(xmin_iq, xmax_iq);
    // Fill data from data_sc into hDCCT_time starting from the last element
    for (int i = data_iq.size() - 1; i >= 0; i--) {
      const auto& entry = data_iq[i];
      if (entry.time < xmin_iq) {
	// Stop filling when entry.time is less than xmax_sc - nx_sc
	break;
      }
      hI_avr_time ->Fill(entry.time - time_offset, entry.I_avr);
      hQ_avr_time ->Fill(entry.time - time_offset, entry.Q_avr);
      hIQ_avr_time->Fill(entry.time - time_offset, entry.IQ_avr);
      hArea_time  ->Fill(entry.time - time_offset, entry.Area);
      hWidth_time ->Fill(entry.time - time_offset, entry.Width);
      hMean_time  ->Fill(entry.time - time_offset, entry.Mean);
    }
    cout<<"chenrj 3..."<<endl;
    std::vector<MyData_iq> data_iq_injection_new =  find_injection(data_sc, data_iq, data_iq_injection, time_offset);
    create_accumulated_spectrum(data_iq_injection_new, frameID_offset, frameID_range);
    cout<<"chenrj 3.1..."<<endl;
    cout<<" data_iq_injection.size() = "<<data_iq_injection.size()<<endl;
    cout<<"chenrj 4..."<<endl;
  }
  
}

TFolder *ReadFolderPointer(TSocket * fSocket)
{
  //read pointer to current folder
  TMessage *m = 0;
  fSocket->Recv(m);
  POINTER_T p;
  *m >> p;
  return (TFolder *) p;
}

/*------------------------------------------------------------------*/
void root_server_thread(void *arg)
/*
  Serve histograms over TCP/IP socket link
*/
{
  char request[256];
  TSocket *sock = (TSocket *) arg;
  //signal(SIGUSR1, SIG_IGN);

  do {
    /* close connection if client has disconnected */
    if (sock->Recv(request, sizeof(request)) <= 0) {
      // printf("Closed connection to %s\n", sock->GetInetAddress().GetHostName());
      sock->Close();
      delete sock;
      return;
    } else {
      TMessage *message = new TMessage(kMESS_OBJECT);
      if (strcmp(request, "GetListOfFolders") == 0) {
	TFolder *folder = ReadFolderPointer(sock);
	if (folder == NULL) {
	  message->Reset(kMESS_OBJECT);
	  message->WriteObject(NULL);
	  sock->Send(*message);
	  delete message;
	  continue;
	}
	//get folder names
	TObject *obj;
	TObjArray *names = new TObjArray(100);
	pthread_mutex_lock(&mutex1);
	TCollection *folders = folder->GetListOfFolders();
	TIterator *iterFolders = folders->MakeIterator();
	while ((obj = iterFolders->Next()) != NULL)
	  names->Add(new TObjString(obj->GetName()));
	//write folder names
	message->Reset(kMESS_OBJECT);
	message->WriteObject(names);
	sock->Send(*message);
	pthread_mutex_unlock(&mutex1);
	for (int i = 0; i < names->GetLast() + 1; i++)
	  delete(TObjString *) names->At(i);
	delete names;
	delete message;
      } else if (strncmp(request, "FindObject", 10) == 0) {
	TFolder *folder = ReadFolderPointer(sock);
	//get object
	TObject *obj;
	if (strncmp(request + 10, "Any", 3) == 0)
	  obj = folder->FindObjectAny(request + 14);
	else
	  obj = folder->FindObject(request + 11);
	//write object
	if (!obj)
	  sock->Send("Error");
	else {
	  message->Reset(kMESS_OBJECT);
	  pthread_mutex_lock(&mutex1);
	  message->WriteObject(obj);
	  sock->Send(*message);
	  pthread_mutex_unlock(&mutex1);
	}
	delete message;
      } else if (strncmp(request, "FindFullPathName", 16) == 0) {
	TFolder *folder = ReadFolderPointer(sock);
	//find path
	const char *path = folder->FindFullPathName(request + 17);
	//write path
	if (!path) {
	  sock->Send("Error");
	} else {
	  TObjString *obj = new TObjString(path);
	  message->Reset(kMESS_OBJECT);
	  message->WriteObject(obj);
	  sock->Send(*message);
	  delete obj;
	}
	delete message;
      } else if (strncmp(request, "Occurence", 9) == 0) {
	TFolder *folder = ReadFolderPointer(sock);
	//read object
	TMessage *m = 0;
	sock->Recv(m);
	TObject *obj = ((TObject *) m->ReadObject(m->GetClass()));
	//get occurence
	Int_t retValue = folder->Occurence(obj);
	//write occurence
	message->Reset(kMESS_OBJECT);
	*message << retValue;
	sock->Send(*message);
	delete message;
      } else if (strncmp(request, "GetPointer", 10) == 0) {
	//find object
	TObject *obj = gROOT->FindObjectAny(request + 11);
	//write pointer
	message->Reset(kMESS_ANY);
	POINTER_T p = (POINTER_T) obj;
	*message << p;
	sock->Send(*message);
	delete message;
      } else if (strncmp(request, "Command", 7) == 0) {
	char objName[100], method[100];
	sock->Recv(objName, sizeof(objName));
	sock->Recv(method, sizeof(method));
	TObject *object = gROOT->FindObjectAny(objName);
	if (object && object->InheritsFrom(TH1::Class())
	    && strcmp(method, "Reset") == 0)
	  static_cast < TH1 * >(object)->Reset();
      } else if (strncmp(request, "SetCut", 6) == 0) {
	//read new settings for a cut
	char name[256];
	sock->Recv(name, sizeof(name));
	TCutG *cut = (TCutG *) gManaHistosFolder->FindObjectAny(name);
	TMessage *m = 0;
	sock->Recv(m);
	TCutG *newc = ((TCutG *) m->ReadObject(m->GetClass()));
	if (cut) {
	  newc->TAttMarker::Copy(*cut);
	  newc->TAttFill::Copy(*cut);
	  newc->TAttLine::Copy(*cut);
	  newc->TNamed::Copy(*cut);
	  cut->Set(newc->GetN());
	  for (int i = 0; i < cut->GetN(); ++i) {
	    cut->SetPoint(i, newc->GetX()[i], newc->GetY()[i]);
	  }
	} 
	delete newc;
      } else
	printf("SocketServer: Received unknown command \"%s\"\n", request);
    }
  } while (1);
  return;
}

/*------------------------------------------------------------------*/

void root_socket_server(void *arg)
{
  // Server loop listening for incoming network connections on specified port.
  // Starts a searver_thread for each connection.
  int port;
  port = *(int *) arg;
  printf("Root server listening on port %d...\n", port);
  TServerSocket *lsock = new TServerSocket(port, kTRUE);
  do {
    TSocket *sock = lsock->Accept();
    // printf("Established connection to %s\n", sock->GetInetAddress().GetHostName());
    TThread *thread = new TThread("Server", root_server_thread, sock);
    thread->Run();
  } while (1);
  return;
}

/*------------------------------------------------------------------*/
void start_root_socket_server(int port)
{
  static int pport = port;
  TThread *thread = new TThread("server_loop", root_socket_server, &pport);
  thread->Run();
}

/*------------------------------------------------------------------*/
/*-- analyzer init routine -----------------------------------------*/
INT mana_init()
{
  // Create an ofstream object and open the file for writing
  std::ofstream outputFile("CSR_Online_Server_status.txt");
  if (!outputFile.is_open())
    {
      std::cerr << "Error opening file: CSR_Online_Server_status.txt " << std::endl;
      exit(1);
    }
  
  char str[256];
  const char *name = "onlie_analy";
  sprintf(str, "Histos for %s", name);
  histo_folder = (TFolder *) gROOT->FindObjectAny(name);
  if (!histo_folder)
    histo_folder =
      gManaHistosFolder->AddFolder(name, str);
  else if (strcmp(((TObject *) histo_folder)->ClassName(), "TFolder")
	   != 0) {
    exit(1);
  }
  gHistoFolderStack->Clear();
  gHistoFolderStack->Add((TObject *) histo_folder);

  create_hist();

  // Write the return value to the file
  outputFile << "Return value: " << SUCCESS << std::endl;
  outputFile.close();
  return SUCCESS;
}
/*------------------------------------------------------------------*/
INT loop_online()
{
  INT status = SUCCESS;
  int ch, sig;
  sigset_t set;

  sigemptyset(&set);
  sigaddset(&set, SIGUSR1);
  printf("Running analyzer online. Stop with \"!\"\n");
  /* main loop */
  do {
    /* check keyboard */
    ch = 0;
    if (ss_kbhit()) {
      ch = getchar();
      if ((char)ch == '!') {
	printf("exiting...\n");
	break;
      }
    }
    /* do some real work here*/
    /*************************/
    sigwait(&set, &sig);  
    pthread_mutex_lock(&mutex1);
    update_hists();
    pthread_mutex_unlock(&mutex1);
  } while (1);
  return status;
}

/*------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int argn = 1;
  char *argp = (char *) argv[0];
  int root_port = 9092;
  sigset_t set;

  sigemptyset(&set);
  sigaddset(&set, SIGUSR1);

  /* block signal */
  if(pthread_sigmask(SIG_BLOCK, &set, NULL) == -1) 
    perror("set process signal o be set failed\n");
  sum = 0;
  manaApp = new TRint("ranalyzer", &argn, &argp, NULL, 0, true);
  /* create the folder for analyzer histograms */
  gManaHistosFolder = gROOT->GetRootFolder()->AddFolder("histos", "MIDAS Analyzer Histograms");
  gHistoFolderStack = new TObjArray();
  gROOT->GetListOfBrowsables()->Add(gManaHistosFolder, "histos");

  /* start socket server */
  start_root_socket_server(root_port);

  /* analyzer init function */
  if (mana_init() != CM_SUCCESS) {
    fprintf(stderr, "initial failed! \n");
    return 1;
  }
  /*---- start main loop ----*/
  loop_online();

  return 0;
}
