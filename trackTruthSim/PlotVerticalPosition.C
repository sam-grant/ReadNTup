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

void Run(TTree *tree, TFile *output, bool quality, bool truth) { 

  // Set the number of periods for the longer modulo plots
  int moduloMultiple = 4; 

  double boostFactor = 5e3*(1/gmagic);

  string stns[] = {"S0S12S18", "S12S18", "S0", "S12", "S18"}; 

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  TH1D *Y_[n_stn];

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    Y_[i_stn] = new TH1D((stns[i_stn]+"_Y").c_str(), ";Vetical decay position [mm];Tracks", 180, -60, 60); 

  }

  // Get branches (using header file)
  InitBranches br(tree);

  double targetPerc = 0;
  int64_t nEntries = tree->GetEntries();

  double muAngleMax = 0;
  
  // For FR randomisation
  TRandom3 *rand = new TRandom3(12345);

  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    int stn = br.station; 

    bool hitVol = br.hitVolume;
    bool vertexQual = br.passVertexQuality;

    double time;
    double y;

    // Time is already in us
    if(truth) { 

      time = br.trueTime; // 
      y = br.trueVertexPosY;

    } else { 

      time = br.recoTime;
      y = br.recoVertexPosY;

    }

    // Time cuts
    if (quality) {

      if (time < g2Period*7) continue; 

      // Should include hit vol and pval cuts
      if (!vertexQual) continue;

    }

    // Fill stations invidually according the station array
    int stn_id = -1;
    if(stn==0) stn_id = 2;
    else if(stn==12) stn_id = 3;
    else if(stn==18) stn_id = 4;
    else cerr<<"Station "<<stn<<" not recognised";
    
    Y_[stn_id]->Fill(y);

    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }


  }

  // Combine stations S12&S18
  Y_[1]->Add(Y_[3], Y_[4]);

  // Combine stations S0&S12&S18
  Y_[0]->Add(Y_[1], Y_[2]);


  // Write to output
  // Set output directory
  output->mkdir("SimultaneousAnalysis"); output->mkdir("MomentumBinnedAnalysis");  

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    output->cd("SimultaneousAnalysis");

    Y_[i_stn]->Write();

  }

  cout<<"Written plots."<<endl;

  return;

}

int main(int argc, char *argv[]) {

  bool quality = true;
  bool truth = true;

  string inFileName = argv[1]; 
  string outFileName = argv[2];

  cout<<inFileName<<endl;
  string treeName = "trackerNTup/TrackerMCDecayTree";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName .c_str());
  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;

  // Book output

  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");
   
  // Fill histograms
  Run(tree, fout, quality, truth);

  // Close
  fout->Close();
  fin->Close();

  cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;
   
  return 0;
}
