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


void Run(TTree *tree, TFile *output, bool truth) {

  // Get branches (using header file)
  // Get branches (using header file)
  InitBranches br(tree);

  int64_t nEntries = tree->GetEntries();
  
  TH1D *momentum_raw = new TH1D("Momentum_Raw", ";Momentum [MeV];Vertices", int(350), 0, 3500); 
  TH1D *momentum_cuts = new TH1D("Momentum_Cuts", ";Momentum [MeV];Vertices", int(350), 0, 3500);

  int64_t nQualEntries = 0; 

  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    double time;
    double x; double y; double z; 
    double px; double py; double pz; 

    // Time is already in us
    if(truth) { 
      time = br.trueTime; // 
      px = br.trueVertexMomX; py = br.trueVertexMomY; pz = br.trueVertexMomZ;
    } else { 
      time = br.recoTime;
      px = br.recoVertexMomX; py = -br.recoVertexMomY; pz = br.recoVertexMomZ;
    }

    TVector3 eMom(px, py, pz);

    double p = eMom.Mag();

    bool vertexQual = br.passVertexQuality;

    if(!vertexQual) continue;

    momentum_raw->Fill(p);
   
    if (time < g2Period*7) continue; 
    if (p < pLo || p > pHi) continue;

    momentum_cuts->Fill(p);

  }

  output->cd();

  momentum_raw->Write();
  momentum_cuts->Write();

  return;

}

int main(int argc, char *argv[]) {

  bool truth = true;

  string inFileName = argv[1]; 
  string outFileName = argv[2];

  string treeName = "trackerNTup/TrackerMCDecayTree";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName.c_str());
  TFile *fout= new TFile(outFileName.c_str(), "RECREATE");

  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;
   
  // Fill histograms
  Run(tree, fout, truth);

  // Close
  fin->Close();
  fout->Close();

  cout<<"\nDone."<<endl;
   
  return 0;

}
