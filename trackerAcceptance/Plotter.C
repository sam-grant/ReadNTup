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

void Run(TTree *treeAllDecays, TTree *treeTracks, TFile *output, bool quality) {

  string stns[] = {"S0S12S18", "S12S18", "S0", "S12", "S18"}; 

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  TH2D *thetaY_vs_Y_tracks_[n_stn]; 
  TH2D *thetaY_vs_Y_allDecays = new TH2D("ThetaY_vs_Y", ";Decay y-position [mm]; #theta_{y} [mrad] / 1 mm ", 180, -60, 60, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic); 

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    thetaY_vs_Y_tracks_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Y").c_str(), ";Decay y-position [mm]; #theta_{y} [mrad] / 1 mm ", 180, -60, 60, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic); 

  }

  // Get branches (using header file)
  InitBranchesAllDecays allDecays(treeAllDecays);

  for(int64_t entry = 0; entry < treeAllDecays->GetEntries(); entry++) {

    treeAllDecays->GetEntry(entry);

    double y = allDecays.posiInitPosY;
    double p = allDecays.posiInitP;
    double py = allDecays.posiInitPY;
    double theta_y = asin(py/p) * 1e3; 

    thetaY_vs_Y_allDecays->Fill(y, theta_y);

  }

  InitBranchesTracks tracks(treeTracks);


  for(int64_t entry = 0; entry < treeAllDecays->GetEntries(); entry++) {

    treeTracks->GetEntry(entry);

    if(quality && !tracks.passVertexQuality) continue;

    double y = tracks.decayVertexPosY;
    double p = tracks.decayVertexMom;
    double py = tracks.decayVertexMomY;
    double theta_y = asin(py/p) * 1e3; 

    int stn = tracks.station;

    int stn_id = -1;
    if(stn==0) stn_id = 2;
    else if(stn==12) stn_id = 3;
    else if(stn==18) stn_id = 4;
    else {
      // cerr<<"Station "<<stn<<" not recognised";
      cout<<"WARNING: Station "<<stn<<" not recognised"<<endl;
    }

    thetaY_vs_Y_tracks_[stn_id]->Fill(y, theta_y);

  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

  // Combine stations S12&S18
  thetaY_vs_Y_tracks_[1]->Add(thetaY_vs_Y_tracks_[3], thetaY_vs_Y_tracks_[4]);
  // Combine stations S0&S12&S18
  thetaY_vs_Y_tracks_[0]->Add(thetaY_vs_Y_tracks_[1], thetaY_vs_Y_tracks_[2]);

  // Write to output
  // Set output directory
  output->mkdir("AllDecays"); output->cd("AllDecays");

  thetaY_vs_Y_allDecays->Write();

  output->mkdir("Tracks"); output->cd("Tracks");

  for (int i_stn = 0; i_stn < n_stn; i_stn++) {

    thetaY_vs_Y_tracks_[i_stn]->Write();

  }


  cout<<"Written plots."<<endl;

  return;

}

int main() { // int argc, char *argv[]) {

  bool quality = true;

  string inFileName = "trackerAcceptanceTrees.test.root"; // argv[1]; 
  string outFileName = "trackerAcceptancePlots.test.root"; // argv[2];

  string treeNameAllDecays = "trackerAcceptanceTrees/AllDecays";
  string treeNameTracks = "trackerAcceptanceTrees/TrackReco";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName .c_str());
  // Get tree
  TTree *treeAllDecays = (TTree*)fin->Get(treeNameAllDecays.c_str());
  TTree *treeTracks = (TTree*)fin->Get(treeNameTracks.c_str());

  cout<<"\nOpened tree:\t"<<treeNameAllDecays<<" "<<treeAllDecays<<" from file "<<inFileName<<" "<<fin<<endl;
  cout<<"\nOpened tree:\t"<<treeNameTracks<<" "<<treeTracks<<" from file "<<inFileName<<" "<<fin<<endl;

  // Book output

  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");
   
  // Fill histograms
  Run(treeAllDecays, treeTracks, fout, quality);

  // Close
  fout->Close();
  fin->Close();

  cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;
   
  return 0;
}
