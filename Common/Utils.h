/*

Samuel Grant

Common includes and constants.

*/

#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TRandom3.h"

using namespace std;

// Constants (updated Oct 2023)
double omega_a = 1.439337; // This is the unweighted average for FNAL Run-2/3 in rad/us, where the uncertainty is 9.521669287914112e-07 rad/us. I used 1.439311 rad/us from BNL in my thesis. 
double g2Period = TMath::TwoPi() / omega_a; // us
double mMu = 105.6583715; // MeV
double aMu = 116592059e-11; // This is the combined BNL and FNAL Run-2/3 result. Previously used 116592089e-11 from BNL. 
double gmagic = std::sqrt( 1.+1./aMu );
double pmax = 1.01 * mMu * gmagic; 
double eMass = 0.510999;
double T_c = 149.2 * 1e-3; // cyclotron period [us]
double R_magic = 7112; // mm

// Common functions

// Open ROOT tree and load branches
TTree *InitTree(string fileName, string treeName) { 

  // Get file
  TFile *fin = TFile::Open(fileName.c_str());
  cout<<"\nOpened tree:\t"<<fileName<<" "<<fin<<endl;

  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  cout<<"\nOpened tree:"<<treeName<<" "<<tree<<" from file "<<fileName<<" "<<fin<<endl;

  return tree;

}

// Modulate time over some number of g-2 periods
double ModTime(double time, int nPeriods = 1) {

  double g2fracTime = time / (nPeriods * g2Period);
  int g2fracTimeInt = int(g2fracTime);
  double g2ModTime = (g2fracTime - g2fracTimeInt) * nPeriods * g2Period;

  return g2ModTime;

}

// Randomise out the fast rotation
double RandomisedTime(TRandom3 *rand, double time) { 
  return rand->Uniform(time-T_c/2, time+T_c/2);
}