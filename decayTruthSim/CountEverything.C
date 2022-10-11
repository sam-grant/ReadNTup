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
#include "TGraph2D.h"

using namespace std;

double omega_a = 1.439311; // (average for mu+ at BNL) 0.00143934*1e3;//1.439311; // rad/us 0.00143934; // kHz from gm2const, it's an angular frequency though...
double g2Period = TMath::TwoPi() / omega_a; // us
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

  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  return tree;

}

void Run(TTree *tree, TFile *output) {

  int64_t nEntries = tree->GetEntries();

  // Entries is "tracks", one muon is a fill
  cout<<nEntries<<endl;

  return;

}


int main(int argc, char *argv[]) {

  string inFileName = argv[1]; 

  string treeName = "phaseAnalyzer/g2phase";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName.c_str());

  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  Run(tree, fin);

  fin->Close();

  return 0;

}

