#include <TSystem.h> // needed for AccessPathName()
#include <TFile.h>
#include <TBranch.h>
#include <TString.h>
#include <TH2.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

Long64_t process(TString);

Long64_t minTime=100000000000000;
Long64_t maxTime=0;
Long64_t totalEntries=0;

void dataRates(Int_t runNo)
{

	minTime=100000000000000;
	maxTime=0;
	totalEntries=0;

	TSystemDirectory dir("./", "./");
   TList *files = dir.GetListOfFiles();
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (!file->IsDirectory() && fname.EndsWith(".root") && fname.BeginsWith(Form("KCLrun_%i_",runNo))) {
         // if (!file->IsDirectory() && fname.EndsWith(".root") && fname.BeginsWith("KCLrun_"+runNo)) {
            cout << fname.Data() << endl;
						// readSpe2TH1F(fname.Data());
					 totalEntries+=process(fname);
         }
      }
   }

	 cout << "Start: " << minTime << endl;
	 cout << "End: " << maxTime << endl;

	 Float_t runLength = (Float_t)(maxTime-minTime)*1e-12;
	 cout << "Run length: " << runLength << " seconds." << endl;
	 cout << "Entries: " << totalEntries << endl;
	 cout << "Hit rate: " << (Float_t)totalEntries/runLength/1e6 << " M hits/s." << endl;
}

Long64_t process(TString name)
{

  TFile *f = new TFile(name);
  TTree *t1 = (TTree*)f->Get("data");

  Float_t         step1;	// value of parameter 1 (when doing parameter scans)
  Float_t         step2;	// value of parameter 2 (when doing parameter scans)
  Long64_t        time;	// time of event, in picoseconds
  UInt_t          channelID;
  Float_t         tot;		// In QDC mode this is the integration time in nanoseconds
  Float_t         energy;
  UShort_t        tacID;
  Int_t           xi;
  Int_t           yi;
  Float_t         x;
  Float_t         y;
  Float_t         z;
  Float_t         tqT;	// the fine timing of the chn (crossing time threshold), in TDC clock units (200 MHz => 5ns
  Float_t         tqE;

  // List of branches
  TBranch        *b_step1;   //!
  TBranch        *b_step2;   //!
  TBranch        *b_time;   //!
  TBranch        *b_channelID;   //!
  TBranch        *b_tot;   //!
  TBranch        *b_energy;   //!
  TBranch        *b_tacID;   //!
  TBranch        *b_xi;   //!
  TBranch        *b_yi;   //!
  TBranch        *b_x;   //!
  TBranch        *b_y;   //!
  TBranch        *b_z;   //!
  TBranch        *b_tqT;   //!
  TBranch        *b_tqE;   //!

  t1->SetBranchAddress("step1", &step1, &b_step1);
  t1->SetBranchAddress("step2", &step2, &b_step2);
  t1->SetBranchAddress("time", &time, &b_time);
  t1->SetBranchAddress("channelID", &channelID, &b_channelID);
  t1->SetBranchAddress("tot", &tot, &b_tot);
  t1->SetBranchAddress("energy", &energy, &b_energy);
  t1->SetBranchAddress("tacID", &tacID, &b_tacID);
  t1->SetBranchAddress("xi", &xi, &b_xi);
  t1->SetBranchAddress("yi", &yi, &b_yi);
  t1->SetBranchAddress("x", &x, &b_x);
  t1->SetBranchAddress("y", &y, &b_y);
  t1->SetBranchAddress("z", &z, &b_z);
  t1->SetBranchAddress("tqT", &tqT, &b_tqT);
  t1->SetBranchAddress("tqE", &tqE, &b_tqE);
////////////////////

  Long64_t nentries = t1->GetEntriesFast();
	Long64_t firstTime,lastTime;
	t1->GetEntry(0);
	firstTime=time;
	if(firstTime<minTime) minTime=firstTime;
	t1->GetEntry(nentries-1);
	lastTime=time;
	if(lastTime>minTime) maxTime=lastTime;
	// Float_t runLength = (Float_t)(lastTime-firstTime)*1e-12;

	// cout << "Run length: " << runLength << " seconds." << endl;
	// cout << "Entries: " << nentries << endl;
	// cout << "Hit rate: " << (Float_t)nentries/runLength/1e6 << " M hits/s." << endl;
	f->Close();

	return nentries;
} // End of Process()
