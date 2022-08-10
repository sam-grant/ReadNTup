// Read ROOT trees 
// Sam Grant

// Count numbers of quality tracks
#include <iostream>
#include <vector>

#include "Plotter.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TRandom3.h"

using namespace std;

double omega_a = 1.439311;// (average for mu+ at BNL) 0.00143934*1e3;//1.439311; // rad/us 0.00143934; // kHz from gm2const, it's an angular frequency though...
double g2Period = TMath::TwoPi() / omega_a;//s * 1e-3; // us

double pLo = 1000; 
double pHi = 2500;

TTree *InitTree(string fileName, string treeName) { 

  // ++++++++++++++ Open tree and load branches ++++++++++++++
  // Get file
  TFile *fin = TFile::Open(fileName.c_str());
  cout<<"\nOpened tree:\t"<<fileName<<" "<<fin<<endl;

  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  cout<<"\nOpened tree:"<<treeName<<" "<<tree<<" from file "<<fileName<<" "<<fin<<endl;

  return tree;

}

void print(std::set<long long> const &s) {
  for (auto const& i: s) {
    std::cout << i << "\n";
  }
}

void Run(TTree *tree, bool truth) {

  std::set<long long> id;

  // Get branches (using header file)
  // Get branches (using header file)
  InitBranches br(tree);

  int64_t nEntries = tree->GetEntries();
  
  //TH1D *id = new TH1D("id", "", int(350), 0, 3500); 
  //TH1D *momentum_cuts = new TH1D("Momentum_Cuts", ";Momentum [MeV];Vertices", int(350), 0, 3500);

  int64_t nQualEntries = 0; 

  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);
    
    if(br.station != 18) continue;

    long long run = br.runNum;
    long long subRun = br.subRunNum;
    long long event = br.eventNum;

    id.insert( run*1000000+subRun*1000+event );

/*    
    int station = br.station;
    float trueTime = br.trueTime; 

    cout<<"\nrunNum: "<<runNum<<endl;
    cout<<"subRunNum: "<<subRunNum<<endl;
    cout<<"eventNum: "<<eventNum<<endl;
    cout<<"station: "<<station<<endl;
    cout<<"trueTime: "<<trueTime<<endl;

    if(station==0) station = 1;

    double uniqueID = run*1000000+subRun*1000+eventdouble(station)*double(subRunNum)*double(eventNum)*double(trueTime);

    cout<<uniqueID<<endl
*/


  }

  print(id);

 //output->cd();
  //id->Write();

  return;

}

int main(int argc, char *argv[]) {

  bool truth = true;

  string inFileName = argv[1]; 
  string outFileName = argv[2];

  string treeName = "trackerNTup/TrackerMCDecayTree";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName.c_str());
  //TFile *fout= new TFile(outFileName.c_str(), "RECREATE");

  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  //cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;
    
  // Fill histograms
  Run(tree, truth);

  // Close
  fin->Close();
  //fout->Close();

  //cout<<"\nDone."<<endl;
   
  return 0;

}
