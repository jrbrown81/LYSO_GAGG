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

const Int_t nPixels=2048;
const Int_t pulserChan=500; // check this!!!
// Some cuts...
const Long64_t timeWindow = 22.; // in ns. This should get all hw_trig coincidences
const Int_t firstHit=1; // switch to choose default first hit. 0=highest energy, 1=second highest. // Julien's sims suggest largest hit is second
const Float_t eCutMin=150;
const Float_t eCutMax=1000;
const Float_t sumThresh=5;
///////////////////////

///////////////////////////////////////////////////////
// Keep all this global stuff so funtions aren't broken
Bool_t calibrated=0;

vector<Double_t> calCoef[nPixels];
Float_t timeShift[nPixels];
vector<Double_t> gainCorr[nPixels];
vector<Float_t> gainCoef[nPixels];
vector<Int_t> gainRuns;
TH2F *runVsEnergy_h[nPixels];
TH2F *runVsEnergy_all_h;
TH2F *runVsEnergy_TOP_h;
TH2F *runVsEnergy_BOT_h;

vector<Int_t> tempRuns_v;
vector<Float_t> temp_v[8][4];
vector<Float_t> calDatime_v;
vector<Float_t> photoPeakPos_v[nPixels];
Float_t timePeak[nPixels];

void Process(TString pathOut,TString pathIn,TString name, Long64_t toProcess);
Bool_t ReadCalFile(TString,Bool_t);
Bool_t ReadTimeCalFile(TString);
// Bool_t ReadGainCorrFile(TString);
Bool_t ReadGainCoefFile(TString);
Bool_t ReadPETsysTempFile(TString filename, int mod2Plot=-1, bool requireRunNo=1);
Float_t GetPhiAngleDeg(const Float_t, const Float_t);
Float_t GetPhiError(Float_t pitch, Float_t dist);
Float_t GetPhiError(Float_t pitch, Float_t dx, Float_t dy);
Float_t GetAngleFromEnergy(const Float_t, const Float_t);
Double_t CalculatePhotopeakStats(TH1F*);
Bool_t CalcGainCorrFactors();
Int_t GetRunNumber(TString,const char*,const char*);
void ClearVectors();
TH1F* AcceptanceCorrection(TH1F*,TH1F*,TH1F*);

// Keep all this global stuff so funtions aren't broken
///////////////////////////////////////////////////////

void LYSOonline(TString name="", TString pathIn="./", TString pathOut="./", Long64_t toProcess = 0)
{
	if(!gSystem->AccessPathName(pathIn+name)) {
  	cout << endl << "Processing single file: " << pathIn+name << endl;
		Process(pathOut,pathIn,name,toProcess);
	} else cout << "File: " << pathIn+name << " does not exists." << endl;
}

void Process(TString pathOut,TString pathIn,TString name, Long64_t toProcess)
{
	Bool_t goodEvent=kFALSE;
	Bool_t goodTOP=kFALSE;
	Bool_t goodBOT=kFALSE;
	Bool_t goodEnergyTOP=kFALSE;
	Bool_t goodEnergyBOT=kFALSE;
	Bool_t goodTimeTOP=kFALSE;
	Bool_t goodTimeBOT=kFALSE;

  TFile *f = new TFile(pathIn+name);
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

  Float_t eCal;
  Long64_t eventStart;
  Long64_t      prevTime=0;
  Int_t          hits=0;
  Int_t          hitsTOP=0;
  Int_t          hitsBOT=0;
	Int_t 				portID=0;

	vector<Int_t> portVector;
  vector<Float_t> energyVectorAll;
  vector<Float_t> energyVectorTOP;
  vector<Float_t> energyVectorBOT;
  // vector<Int_t> energyIndexTOP;
  // vector<Int_t> energyIndexBOT;
  vector<Int_t> channelVector;
  vector<Long64_t> timeVector;
  vector<Float_t> eventTime_v;

  Float_t energyTOP=0;
  Float_t energyBOT=0;
////////////////////

// Histograms
	Float_t eMin4Spec=-5;
	Float_t eMax4Spec=50;
	Int_t eBins1D=1100;
	Int_t eBins2D=110;

	TString outfile="output_" + name;
	cout << "Writing output to " << pathOut+outfile << endl;
  TFile *f2 = new TFile(pathOut+outfile,"recreate");
	f2->cd();
  TH1F* hits_h = new TH1F("hits","Number of hits per event",20,0,20);
  TH1F* hitsTOP_h = new TH1F("hitsTOP","Number of hits per event in top array",20,0,20);
  TH1F* hitsBOT_h = new TH1F("hitsBOT","Number of hits per event in bottom array",20,0,20);

	TH1F* hitPatternTOP_h = new TH1F("hitPatternTOP","Hit pattern for top array",1024,1024,nPixels);
	TH1F* hitPatternBOT_h = new TH1F("hitPatternBOT","Hit pattern for bottom array",1024,0,1024);

  TH1F* hitTime_h=new TH1F("hitTime","Time difference between hits;time difference (ns); counts",1000,-20,80);
  TH1F* eventTime_h=new TH1F("eventTime","Hit time relative to the first hit in the event;time difference (ns); counts",1000,-20,80);

  TH1F* energyTOP_h = new TH1F("energyTOP","Total energy in top array;Energy (ADC);Counts /bin",eBins1D,eMin4Spec,eMax4Spec);
  TH1F* energyBOT_h = new TH1F("energyBOT","Total energy in bottom array;Energy (ADC);Counts /bin",eBins1D,eMin4Spec,eMax4Spec);

	TH1F* pulserEnergy_h = new TH1F("pulserEnergy",Form("Pulser energy spectrum (channel %i)",pulserChan),eBins1D,eMin4Spec,eMax4Spec);
	TH1F* pulserTime_h = new TH1F("pulserTime",Form("Pulser time, relative to first hit in the event (channel %i)",pulserChan),100,-25,25);
  TH1F* energySpec_h[nPixels];
  TH1F* timeSpec_h[nPixels];
  for(int i=0;i<nPixels;i++) {
    energySpec_h[i]=new TH1F(Form("energySpec_%i",i),Form("Energy Spectrum for channel %i",i),eBins2D,eMin4Spec,eMax4Spec);
		timeSpec_h[i]=new TH1F(Form("timeSpec_%i",i),Form("Hit time relative to the first hit in the event for channel %i",i),100,-25,25);
  }

  TH2F* chnVsEnergy_h = new TH2F("chnVsEnergy","ChannelID vs. Energy (ADC);Energy (ADC);channelID",eBins2D,eMin4Spec,eMax4Spec,nPixels,0,nPixels);
  TH2F* chnVsEnergyTOP_h = new TH2F("chnVsEnergyTOP","ChannelID vs. Energy (ADC), Top array;Energy (ADC);channelID",eBins2D,eMin4Spec,eMax4Spec,1024,1024,nPixels);
  TH2F* chnVsEnergyBOT_h = new TH2F("chnVsEnergyBOT","ChannelID vs. Energy (ADC), Bottom array;Energy (ADC);channelID",eBins2D,eMin4Spec,eMax4Spec,1024,0,1024);
  TH2F* chnVsEventTime_h = new TH2F("chnVsEventTime","ChannelID vs. Time since start of event;Time Difference (ns);channelID",1000,-20,80,nPixels,0,nPixels);
////////////////////

  Long64_t nentries = t1->GetEntriesFast();
  cout << "Entries Found " << nentries << endl;
  if (toProcess == 0) {
    toProcess = nentries;
    cout << "Processing all entries... " << endl;
  } else cout << "Processing " << toProcess << " entries... " << endl;
  Long64_t jentry = 0;
  for (jentry=0; jentry<toProcess;jentry++) {
    if(jentry!=0 && jentry%1000000==0) cout << (double)jentry/(double)toProcess*100 << "% done       \r" << flush; // extra space added to clear previous text from the end of the line
    t1->GetEntry(jentry);
    if(jentry==0) eventStart=time; // So first hit will not trigger the filling of the tree etc.

// New event. Fill tree and clear vectors.
    if(time>(eventStart+timeWindow*1000))
		{
			hitsTOP=(int)energyVectorTOP.size();
      if(hitsTOP==2) goodTOP=kTRUE;
      hitsBOT=(int)energyVectorBOT.size();
      if(hitsBOT==2) goodBOT=kTRUE;
      // int indexTOP[hitsTOP];
      // int indexBOT[hitsBOT];

      for(int i=0;i<hitsTOP;i++) {
      	if(portVector[i]==2) {
        	energyVectorTOP.insert(energyVectorTOP.begin(),energyVectorAll[i]);
      	}
        if(energyVectorTOP[i]>sumThresh) {
          energyTOP+=energyVectorTOP[i];
        }
        // energyIndexTOP.push_back(indexTOP[i]);
      }
      for(int i=0;i<hitsBOT;i++) {
      	if(portVector[i]==1) {
        	energyVectorBOT.insert(energyVectorBOT.begin(),energyVectorAll[i]);
      	}
        if(energyVectorBOT[i]>sumThresh) {
          energyBOT+=energyVectorBOT[i];
        }
        // energyIndexBOT.push_back(indexBOT[i]);
      }

//////
      if(goodTOP && goodBOT && goodEnergyTOP && goodEnergyBOT) goodEvent=kTRUE;

      if(goodEvent){ // goodEvents will include all noScatEvents because they have the same conditions, just noScatEvents requires <7 keV in the scatterer

	    }

// Fill some event histograms
			hits_h->Fill(hits);
      hitsTOP_h->Fill(hitsTOP);
      hitsBOT_h->Fill(hitsBOT);

      for(int i=0;i<hits;i++) {
        chnVsEnergy_h->Fill(energyVectorAll[i],channelVector[i]);
				if(portVector[i]==2) {
					chnVsEnergyTOP_h->Fill(energyVectorAll[i],channelVector[i]);
					hitPatternTOP_h ->Fill(channelVector[i]);
				}
				if(portVector[i]==1) {
					chnVsEnergyBOT_h->Fill(energyVectorAll[i],channelVector[i]);
					hitPatternBOT_h ->Fill(channelVector[i]);
				}
        energySpec_h[channelVector[i]]->Fill(energyVectorAll[i]);
				timeSpec_h[channelVector[i]]->Fill(eventTime_v[i]);

				if(channelVector[i]==pulserChan) {
					pulserEnergy_h->Fill(energyVectorAll[i]);
					pulserTime_h->Fill(eventTime_v[i]);
				}

        if(i>0){ // skip the first hit as this defines t0
          chnVsEventTime_h->Fill(eventTime_v[i],channelVector[i]);
          eventTime_h->Fill(eventTime_v[i]);
        }

      }

      if(energyTOP!=0) {
        energyTOP_h->Fill(energyTOP);
      }
      if(energyBOT!=0) {
        energyBOT_h->Fill(energyBOT);
      }

// clear some variables
			portVector.clear();
      energyVectorAll.clear();
      energyVectorTOP.clear();
      energyVectorBOT.clear();
      channelVector.clear();
      // energyIndexTOP.clear();
      // energyIndexBOT.clear();
      timeVector.clear();
      eventTime_v.clear();
			portID=0;
      hits=0;
      hitsTOP=0;
      hitsBOT=0;
      energyTOP=0;
      energyBOT=0;
      eventStart=time;

      goodEvent=kFALSE;
      goodTOP=kFALSE;
      goodBOT=kFALSE;
      goodEnergyTOP=kFALSE;
      goodEnergyBOT=kFALSE;
      goodTimeTOP=kFALSE;
      goodTimeBOT=kFALSE;

    } // End of new event

		// adjust the channelIDs to a sensible range
		if(channelID>=131072 && channelID<262144) {
			channelID-=131072;
			portID=1; // bottom array
		} else if(channelID>=262144) {
			channelID-=(262144-1024);
			portID=2; // top array
		}
    // calibrate energy
    eCal=0; // reset calibrated energy first
    if(calCoef[channelID].size()!=0) for(int i=0;i<(int)calCoef[channelID].size();i++) eCal+=calCoef[channelID][i]*pow(energy,i);
		else eCal=energy;

		portVector.push_back(portID);
    energyVectorAll.push_back(eCal);
    channelVector.push_back(channelID);
    timeVector.push_back(time);
    eventTime_v.push_back((Float_t)(time-eventStart)/1000.);
    hits++;

    if(jentry>0) hitTime_h->Fill((time-prevTime)/1000.); // first entry has no previous time
    prevTime=time;
  } // End of loop over entries
  cout << (double)jentry/(double)toProcess*100 << "% done       \r" << endl;
	f->Close();

// Write histograms and tree to output file then close file
  hits_h->Write();
  hitsTOP_h->Write();
  hitsBOT_h->Write();

	hitTime_h->Write();
  eventTime_h->Write();
	chnVsEventTime_h->Write();

  energyTOP_h->Write();
  energyBOT_h->Write();

  hitPatternTOP_h->Write();
  hitPatternBOT_h->Write();

	chnVsEnergy_h->Write();
	chnVsEnergyTOP_h->Write();
	chnVsEnergyBOT_h->Write();

	pulserEnergy_h->Write();
	pulserTime_h->Write();

	for(int i=0;i<nPixels;i++) energySpec_h[i]->Write();
	for(int i=0;i<nPixels;i++) timeSpec_h[i]->Write();

	f2->Close();
} // End of Process()

Int_t GetRunNumber(TString fileName, const char* prefix, const char* suffix)
{
  fileName.ReplaceAll(prefix,"");
  fileName.ReplaceAll(suffix,"");
  Int_t runNo=fileName.Atoi();
  return runNo;
}

Bool_t ReadCalFile(TString calFileName, Bool_t recalibrate)
{
  if(recalibrate) {
    cout << "Forcing recalibration..." << endl;
    cout << "Calibraion warning will be suppressed." << endl;
  }

// Read energy calibration file
	string line, word;
	Int_t col=0;
	ifstream calFile(calFileName);

	Int_t chn=-1;
	Int_t nCalChns=0;

	if(calFile.is_open()){
		cout << "Reading calibration from file: " << calFileName << endl;
		while(getline(calFile,line)){
			nCalChns++;
			istringstream iss(line);
			while(getline(iss,word,',')){
				if(col==0) {
					chn=(Int_t)stoi(word);
					if(calCoef[chn].size()!=0) {
						if(!recalibrate) cerr << "Warning! Channel " << chn << " calibration is being overwritten." << endl;
						calCoef[chn].clear();
					}
				}
				else calCoef[chn].push_back((Double_t)stof(word));
				col++;
			}
			chn=-1;
			col=0;
		}
		cout << nCalChns << " calibration channels read." << endl;
	} else {
    cout << "Calibtaion file: " << calFileName << " not found!" << endl;
    return 0;
  }
  return 1;
}

Bool_t ReadTimeCalFile(TString calFileName)
{
// Read time calibration file
	string line, word;
	Int_t col=0;
	ifstream calFile(calFileName);

	Int_t chn=-1;
	Int_t nCalChns=0;

	if(calFile.is_open()){
		cout << endl << "Reading time calibration from file: " << calFileName << endl;
		while(getline(calFile,line)){
			istringstream iss(line);
			while(getline(iss,word,',')){
				if(col==0) chn=(Int_t)stoi(word);
				else timeShift[chn]=(Float_t)stof(word);
				col++;
			}
      if(timeShift[chn]!=0) nCalChns++;
			chn=-1;
			col=0;
		}
		cout << nCalChns << " time shifts read." << endl;
	} else {
    cout << endl << "Time calibtaion file: " << calFileName << " not found!" << endl;
    return 0;
  }
  return 1;
}

Bool_t ReadGainCoefFile(TString gainFileIn)
{
  string line, word;
	ifstream inFile(gainFileIn);

  Int_t chn=-1;
  Int_t col=0;
  Int_t nCalChns=0;

  if(inFile.is_open()){
    cout << endl << "Reading gain coefficients from file: " << gainFileIn << endl;

    while(getline(inFile,line)){
      istringstream iss(line);
      while(getline(iss,word,',')){
        if(col==0) chn=(Int_t)stoi(word);
        else gainCoef[chn].push_back((Float_t)stof(word));
        col++;
      }
      if(!(gainCoef[chn][0]==0 && gainCoef[chn][1]==1)) nCalChns++;
      chn=-1;
      col=0;
    }
    cout << "Gain coefficients read successfully for " << nCalChns << " channels." << endl;
  }
  else cout << endl << gainFileIn << " file not found!" << endl;

  if(nCalChns!=0) return 1;
  else return 0;
}

Bool_t CalcGainCorrFactors()
{
  for(int chn=0;chn<nPixels;chn++){
    Int_t module=chn/128; // choose module 0-8
    Int_t sensor=(chn%128)/64; // choose sensor, 0 for first half of module, 1 for second half
    for(Int_t run=0;run<(Int_t)temp_v[module][sensor].size();run++) gainCorr[chn].push_back(511./(gainCoef[chn][1]*temp_v[module][sensor][run]+gainCoef[chn][0]));
  }
  for(auto element : tempRuns_v) gainRuns.push_back(element); // copy the list of runs with temperature data to a list of runs with gain data

  if(gainCorr[21].size()==gainRuns.size()) return 1;
  else return 0;
}

void ClearVectors()
{
  for(int i=0;i<nPixels;i++) {
    calCoef[i].clear();
    gainCorr[i].clear();
    gainCoef[i].clear();
    photoPeakPos_v[i].clear();
  }

  gainRuns.clear();
  tempRuns_v.clear();
  calDatime_v.clear();

  for(int i=0;i<8;i++) for(int j=0;j<4;j++) temp_v[i][j].clear();
}

Bool_t ReadPETsysTempFile(TString filename, int mod2Plot, bool requireRunNo)
{
  int lineCount=0;
  int recordCount=0;
  int counter=0;
  vector<float> datime_v;
  vector<float> run_v;
  vector<bool> calFlag_v;
  vector<float> calRun_v;
  vector<float> calTemp_v[8][4];

  float temp;
  float tempBuffer[4];
  float run;
  int mod;
  string dummy;

  string line, word;
  ifstream infile(filename);
  TDatime datime;

  Bool_t newRecord=0;
  Bool_t goodRecord=0;
  Bool_t goodTime=0;
  Bool_t goodRun=0;
  Bool_t goodTemp=0;

  if(!infile.is_open())
    cout << "\nTemperature log file '" << filename << "' not found." << endl;
  else {
		cout << "\nReading temperature log from file: " << filename << endl;
    if(!requireRunNo) {
      cout << "----------------------------------------------------" << endl;
      cout << "!!! Warning, run numbers are not being required. !!!" << endl;
      cout << "!!! Any plots with run numbers in should be      !!!" << endl;
      cout << "!!! treated with suspicion                       !!!" << endl;
      cout << "----------------------------------------------------\n" << endl;
    } else cout << "Records without run numbers will be ignored" << endl;
		while(getline(infile,line)){
      lineCount++;
      if(line.find("---")!=string::npos) {
        line.erase(line.begin(),line.begin()+30);
        line.erase(line.end()-30,line.end());
        if(line.find("---")==string::npos) { // date/time line
          datime.Set(&line[0]);
          newRecord=1;
          goodTime=1;
          goodTemp=0;
          goodRecord=0;
        }
      } else if(line.rfind("0",0)==0) { // temperature lines
        istringstream iss(line);
        iss >> dummy >> dummy >> mod >> tempBuffer[0] >> dummy >> tempBuffer[1] >> dummy >> tempBuffer[2] >> dummy >> tempBuffer[3] >> dummy;
        goodTemp=1;
        if((!requireRunNo || goodRun) && goodTime && goodTemp && newRecord) goodRecord=1;
        if(goodRecord) {
          for(int s=0;s<4;s++) temp_v[mod][s].push_back(tempBuffer[s]);
          if(newRecord) {
            datime_v.push_back(datime.Convert());
            run_v.push_back(run);
            tempRuns_v.push_back(run);
            newRecord=0;
            goodRun=0;
            goodTime=0;
            goodTemp=0;
          }
        }
      } else {
        istringstream iss(line);
        while(getline(iss,word,' ')) {
          if(word.find("Run")!=string::npos) { // word contains "Run"
            run=stof(word.substr(3));
            calFlag_v.push_back(0);
            newRecord=1;
            goodRun=1;
            goodRecord=0;
            goodTime=0;
          } else if(word.find("Cal")!=string::npos) { // line contains "Cal"
            calFlag_v.back()=1;
          }
        }
      }
    }
  }

  for(int i=0;i<(int)calFlag_v.size();i++) {
    if(calFlag_v[i]) {
      calRun_v.push_back(run_v[i]);
      calDatime_v.push_back(datime_v[i]);
      for(int m=0;m<8;m++)
        if(temp_v[m][0].size()!=0)
          for(int s=0;s<4;s++)
            calTemp_v[m][s].push_back(temp_v[m][s][i]);
    }
  }

  cout << "Runs: " << run_v.size() << endl;
  cout << "Calibration runs: " << calRun_v.size() << endl;
  cout << "Date and Time records: " <<  datime_v.size() << endl;
  mod=0;
  for(int m=0;m<8;m++) {
    if(temp_v[m][0].size()!=0) {
      mod++;
      if(mod==1) cout << "Temperature records: " << temp_v[m][0].size() << endl;
    }
  }
  cout << "Active modules: " << mod << endl;

  if(mod>0) return 1;
  else return 0;

}

Double_t CalculatePhotopeakStats(TH1F *spectra)
{
  Int_t binMax = spectra->GetMaximumBin();
  Double_t Ephotopeak = spectra->GetXaxis()->GetBinCenter(binMax);
  return Ephotopeak;
}

Float_t GetPhiAngleDeg(const Float_t dx, const Float_t dy)
{
	Float_t ang = -9999;
  ang = TMath::ATan2(dy,dx)*TMath::RadToDeg();
	return ang;
}

Float_t GetPhiError(Float_t pitch, Float_t dx, Float_t dy)
{
	Float_t err;
	Float_t dist=TMath::Sqrt(dx*dx+dy*dy);
	return err=1./TMath::Sqrt(6)*abs(pitch/dist)*TMath::RadToDeg(); // Taken from Makek2019.
}
Float_t GetPhiError(Float_t pitch, Float_t dist)
{
	Float_t err;
	// Float_t dist=TMath::Sqrt(dx*dx+dy*dy);
	return err=1./TMath::Sqrt(6)*abs(pitch/dist)*TMath::RadToDeg(); // Taken from Makek2019.
}

Float_t GetAngleFromEnergy(const Float_t incomingEnergy, const Float_t depositedEnergy)
{
  Float_t eRestMass=511;
	Float_t angle = TMath::ACos(1.+eRestMass/incomingEnergy-eRestMass/(incomingEnergy-depositedEnergy));
  return angle;
}

TH1F* AcceptanceCorrection(TH1F* dPhi_h, TH1F* dPhiAcc_h, TH1F* dPhiCorr_h)
{
	TH1F* dPhi2_h=(TH1F*)dPhi_h->Clone();
	TH1F* dPhiAcc2_h=(TH1F*)dPhiAcc_h->Clone();
	TH1F* dPhiCorr2_h=(TH1F*)dPhiCorr_h->Clone();
	if(dPhi2_h->Integral()>0 && dPhiAcc2_h->Integral()>0) {
		dPhiAcc2_h->Scale(dPhiAcc2_h->GetNbinsX()/dPhiAcc2_h->Integral());
		dPhiCorr2_h->Divide(dPhi2_h,dPhiAcc2_h);
	} else {
		if(dPhi2_h->Integral()<=0) cout << dPhi2_h->GetTitle() << " has integral <=0" << endl;
		if(dPhiAcc2_h->Integral()<=0) cout << dPhiAcc2_h->GetTitle() << " has integral <=0" << endl;
		cout << "Acceptance correction not applied" << endl;
	}
	return dPhiCorr2_h;
}
