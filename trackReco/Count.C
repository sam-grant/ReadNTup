// Read ROOT trees 
// Count raw tracks, quality tracks, raw vertices, quality vertices
// Sam Grant
// For data

#include <iostream>
#include <vector>

#include "Plotter.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TGraphErrors.h"

using namespace std;

double mMu = 105.6583715; // MeV
double aMu = 11659208.9e-10; 
double gmagic = std::sqrt( 1.+1./aMu );
double pmax = 1.01 * mMu * gmagic; // 3127.1144

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

void Run(TTree *tree, TFile *output, string dataset="Run-1a") {

  // Get branches (using header file)
  InitBranches br(tree);

  int64_t nEntries = tree->GetEntries();
  
  TH1D *momentum_raw = new TH1D("Momentum_Raw", ";Momentum [MeV];Tracks", int(350), 0, 3500); 
  TH1D *momentum_qualityTracks = new TH1D("Momentum_QualityTracks", ";Momentum [MeV];Tracks", int(350), 0, 3500); 
  TH1D *momentum_qualityVertices = new TH1D("Momentum_QualityVertices", ";Momentum [MeV];Vertices", int(350), 0, 3500);

  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    TVector3 eMom(br.decayVertexMomX, -br.decayVertexMomY, br.decayVertexMomZ); 

    // Not sure why this would happen?
//    if(eMom.Mag() <= 0) continue;

    momentum_raw->Fill(eMom.Mag());
    if(br.passTrackQuality) momentum_qualityTracks->Fill(eMom.Mag());
    if(br.passDecayVertexQuality) momentum_qualityVertices->Fill(eMom.Mag());

  }

  output->cd();

  momentum_raw->Write();
  momentum_qualityTracks->Write();
  momentum_qualityVertices->Write();

  return;

}

int main(int argc, char *argv[]) { 
  
  string inFileName = argv[1];
  string outFileName = argv[2];
  string dataset = argv[3];
  
   string treeName = "trackAndTrackCalo/tree";

   // Open tree and load branches
   TFile *fin = TFile::Open(inFileName .c_str());
   // Get tree
   TTree *tree = (TTree*)fin->Get(treeName.c_str());

   cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;

   // Book output
   TFile *fout = new TFile(outFileName.c_str(),"RECREATE");
   
   // Fill histograms
   Run(tree, fout, dataset);

   // Close
   fout->Close();
   fin->Close();

   cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;
   
   return 0;

}
