#include <TSystem.h> // needed for AccessPathName()
#include <TFile.h>
#include <TString.h>

using namespace std;

void fileRecov(TString fname)
{
	TString newFile = fname;
	newFile.ReplaceAll(".root","_recov.root");

	cout << "Attempting to recover file '" << fname << "'. Please be patient, this may take a while!" << endl;
	TFile *file = TFile::Open(fname);
	if(file->TestBit(TFile::kRecovered)) { // test if recovery was successful
		TTree* tree=(TTree*)file->Get("data");
		TFile* recov=new TFile(newFile,"RECREATE");
		tree->CloneTree();
		recov->Write();
		recov->Close();

		cout << "Successfully recovered file '" << fname << "' to new file '" << newFile << "'" << endl;
	} else {
		cout << "Failed to recover file ''" << fname << "'" << endl;
	}
	file->Close();
}
