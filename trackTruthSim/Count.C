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


void Run(TTree *tree, bool quality, bool truth) {

  // Get branches (using header file)
  InitBranches br(tree);

  double targetPerc = 0;
  double muAngleMax = 0;

  int64_t nEntries = tree->GetEntries();

  int64_t nQualEntries = 0; 

  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    bool vertexQual = br.passVertexQuality;

    if(vertexQual) nQualEntries++;

    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }

  }

  cout<<"quality vertices=\t"<<nQualEntries<<endl;

  return;

}

int main(int argc, char *argv[]) {

  bool quality = true;
  bool truth = true;

  string inFileName = argv[1]; 

  string treeName = "trackerNTup/TrackerMCDecayTree";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName .c_str());
  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;
   
  // Fill histograms
  Run(tree, quality, truth);

  // Close
  fin->Close();

  cout<<"\nDone."<<endl;
   
  return 0;
}
