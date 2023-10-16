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

double omegaAMagic = 1.439311; // BNL mu+ average // 0.00143934; // from gm2geom consts / kHz 
double g2Period = (2*TMath::Pi()/omegaAMagic); // * 1e-3; // 4.3653239 us
double mMu = 105.6583715; // MeV
double aMu = 11659208.9e-10; 
double gmagic = std::sqrt( 1.+1./aMu );
double pmax = 1.01 * mMu * gmagic; // 3127.1144
double T_c = 149.2 * 1e-3; // cyclotron period [us]

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

  // For direct comp. with decays
  double boostFactor = 5e3*(1/gmagic);

  TH2D *thetaY_vs_p_[n_stn];
  // TH2D *decayX_vs_decayZ_[n_stn];

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    thetaY_vs_p_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Momentum").c_str(), ";Momentum [MeV]; #theta_{y} [mrad] / 50 MeV ", int(pmax)/50, 0, pmax, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
    // decayX_vs_decayZ_[i_stn] = new TH2D((stns[i_stn]+"_DecayZ_vs_DecayX").c_str(), ";Decay vertex position Z [mm];Decay vertex position X [mm]", 800, -8000, 8000, 800, -8000, 8000);

  }

  // Get branches (using header file)
  InitBranches br(tree);

  double targetPerc = 0;
  int64_t nEntries = tree->GetEntries();

  double muAngleMax = 0;

  int64_t counter = 0;
   
  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    // quality variables
    double time = br.decayTime * 1e-3; // ns -> us

    double x = br.decayVertexPosX; double y = br.decayVertexPosY; double z = br.decayVertexPosZ; 
    double px = br.decayVertexMomX; double py = -br.decayVertexMomY; double pz = br.decayVertexMomZ; 
    bool hitVol = br.hitVolume;
    double pVal = br.trackPValue;
    bool vertexQual = br.passDecayVertexQuality;

    // Time cuts
    if (quality) {
      if(time < g2Period*7) continue;
      if (time < g2Period*12 && dataset == "Run-1d") continue; // delayed start time for EG
      if (!vertexQual) continue;
      //if (time > g2Period*15) continue;
    }

    int stn = br.station; 

    // Positron world momentum
    TVector3 eMom(px, py, pz); 
    TVector3 ePos(x, y, z);
      
    /////////////////////////////////
    //                             //
    // Define the vertical angle   //
    //                             //
    /////////////////////////////////

    double p = eMom.Mag();

    // If the you have the vertex cut switched on, py must be negative
    double theta_y = asin(py/p) * 1e3;

    // End of variable definitions.
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    int stn_id = -1;
    if(stn==12) stn_id = 0;
    else if(stn==18) stn_id = 1;
    else cerr<<"Station "<<stn<<" not recognised";

    thetaY_vs_p_[stn_id]->Fill(eMom.Mag(), theta_y);
    // decayX_vs_decayZ_[stn_id]->Fill(ePos.z(), ePos.x());

  }

  // Combine stations
  thetaY_vs_p_[2]->Add(thetaY_vs_p_[0], thetaY_vs_p_[1]);
  // decayX_vs_decayZ_[2]->Add(decayX_vs_decayZ_[0], decayX_vs_decayZ_[1]);

  // Write to output
  // Set output directory
  output->mkdir("SanityPlots"); output->cd("SanityPlots");

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    thetaY_vs_p_[i_stn]->Write();
    // decayX_vs_decayZ_[i_stn]->Write();

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
