// Read ROOT trees 
// Sam Grant

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

double omegaAMagic = 0.00143934; // from gm2geom consts / kHz 
double g2Period = (2*TMath::Pi()/omegaAMagic) * 1e-3; // 4.3653239 us
double mMu = 105.6583715; // MeV
double aMu = 11659208.9e-10; 
double gmagic = std::sqrt( 1.+1./aMu );
double pmax = 1.01 * mMu * gmagic;
double eMass = 0.510999;
double T_c = 149.2 * 1e-3; // cyclotron period [us]

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

void Run(TTree *tree, TFile *output) { 

  // Get branches (using header file)
  InitBranches br(tree);

  TH1D *momentum_raw = new TH1D("Momentum_Raw", ";Momentum [MeV];Decays", int(350), 0, 3500); 
  TH1D *momentum_cuts = new TH1D("Momentum_Cuts", ";Momentum [MeV];Cuts", int(350), 0, 3500); 

  int64_t nEntries = tree->GetEntries();

  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    double time;
    double x; double y; double z; 
    double px; double py; double pz; 


    // Positron world momentum 
    TVector3 eMom(br.posiInitPX, br.posiInitPY, br.posiInitPZ); 

    double p = eMom.Mag();

    momentum_raw->Fill(p);

    if(time < g2Period*7) continue; 
    if(p < pLo || p > pHi) continue;

    momentum_cuts->Fill(p);
         
  }

  output->cd();

  momentum_raw->Write();
  momentum_cuts->Write();


  return;

}


int main(int argc, char *argv[]) {

  
  string inFileName = argv[1]; 
  string outFileName = argv[2];
  string null = argv[3];

  string treeName = "phaseAnalyzer/g2phase";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName.c_str());
  TFile *fout = new TFile(outFileName.c_str(), "RECREATE");
  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;

  Run(tree, fout);
  
  fin->Close();
  fout->Close();

  return 0;
}

