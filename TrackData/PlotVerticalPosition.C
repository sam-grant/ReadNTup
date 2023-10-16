// Read ROOT trees 
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

double omegaAMagic = 0.00143934; // from gm2geom consts / kHz 
double g2Period = (2*TMath::Pi()/omegaAMagic) * 1e-3; // 4.3653239 us
double mMu = 105.6583715; // MeV
double aMu = 11659208.9e-10; 
double gmagic = std::sqrt( 1.+1./aMu );
double pmax = 1.01 * mMu * gmagic; // 3127.1144
double T_c = 149.2 * 1e-3; // cyclotron period [us]

double pLo = 750; 
double pHi = 2750;

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
void Run(TTree *tree, TFile *output, string dataset="Run-1a", bool quality=true) { 


  string stns[] = {"S12", "S18", "S12S18"}; 

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  TH1D *Y_[n_stn];

  for (int i_stn = 0; i_stn < n_stn; i_stn++) Y_[i_stn] = new TH1D((stns[i_stn]+"_Y").c_str(), ";Vertical decay position [mm];Vertices", 180, -60, 60); 

  // Get branches (using header file)
  InitBranches br(tree);

  double targetPerc = 0;
  int64_t nEntries = tree->GetEntries();
   
  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    // quality variables
    bool hitVol = br.hitVolume;
    double pVal = br.trackPValue;
    bool vertexQual = br.passDecayVertexQuality;

    int stn = br.station; 

    double time = br.decayTime * 1e-3; // ns -> us
    double y = br.decayVertexPosY; 

    if (quality) {

      if (time < g2Period*7) continue; 
      if (!vertexQual) continue;

    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////

    int stn_id = -1;
    if(stn==12) stn_id = 0;
    else if(stn==18) stn_id = 1;
    else cerr<<"Station "<<stn<<" not recognised";

    Y_[stn_id]->Fill(y);

    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }
 
  }

  // Combine stations
  Y_[2]->Add(Y_[0], Y_[1]);

  // Write to output
  // Set output directory
  output->mkdir("SimultaneousAnalysis");

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    output->cd("SimultaneousAnalysis");

    Y_[i_stn]->Write();

  }

  cout<<"Written plots."<<endl;

  return;

}

int main(int argc, char *argv[]) {

   bool quality = true;

   string inFileName = argv[1]; 
   string outFileName = argv[2]; 
   string dataset = argv[3];
/*
   string inFileName = "/gm2/data/g2be/Production/Trees/Run1/trackRecoTrees_15921.root";
   string outFileName = "trackRecoPlots_15921_test.root";
   string dataset = "Run-1a";*/

   string treeName = "trackAndTrackCalo/tree";

   // Open tree and load branches
   TFile *fin = TFile::Open(inFileName .c_str());
   // Get tree
   TTree *tree = (TTree*)fin->Get(treeName.c_str());

   cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;

   // Book output

   TFile *fout = new TFile(outFileName.c_str(),"RECREATE");
   
   // Fill histograms
    Run(tree, fout, dataset, quality); 

   // Close
   fout->Close();
   fin->Close();

   cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;
   
   return 0;

}
