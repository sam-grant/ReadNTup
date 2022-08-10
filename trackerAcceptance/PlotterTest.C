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

void Run(TTree *treeTracks, TFile *output, bool quality) {

  // Book histograms: slice y-position; slice momentum

  string stns[] = {"S12S18", "S12", "S18"}; 

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  TH1D *Y_tracks_[n_stn]; 

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 
    
    Y_tracks_[i_stn] = new TH1D((stns[i_stn]+"_Y").c_str(), ";Decay y-position [mm];Tracks", 20, -50, 50);

  }

  InitBranchesTracks tracks(treeTracks);

  cout<<"-----> Filling track histograms"<<endl;

  for(int64_t entry = 0; entry < treeAllDecays->GetEntries(); entry++) {

    treeTracks->GetEntry(entry);

    double y = tracks.decayVertexPosY;
    int stn = tracks.station;

    if(quality) {
      if(!tracks.passVertexQuality) continue;
      if(tracks.pValue < 0.05) continue;
      if(tracks.hitVolume) continue; 
    } 

    int stn_id = -1;
    if(stn==12) stn_id = 1;
    else if(stn==18) stn_id = 2;
    else {
      // cerr<<"Station "<<stn<<" not recognised";
      cout<<"WARNING: Station "<<stn<<" not recognised"<<endl;
      continue;
    }

    if(quality && p > pLo && p < pHi) {

      Y_tracks_[stn_id]->Fill(y);

    } else if(!quality) {

      Y_tracks_[stn_id]->Fill(y);

    }

  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

  // Combine stations S12&S18
  cout<<"-----> Combining stations"<<endl;

  Y_tracks_[0]->Add(Y_tracks_[1], Y_tracks_[2]);

  cout<<"-----> Writing track histograms"<<endl;

  output->mkdir("Tracks"); 

  output->cd("Tracks"); 

  for (int i_stn = 0; i_stn < n_stn; i_stn++) {

    Y_tracks_[i_stn]->Write();

  }

  cout<<"-----> Written plots."<<endl;

  return;

}

int main(int argc, char *argv[]) {

  bool quality = true;

  string inFileName = argv[1]; // "trackerAcceptanceTrees.test.root"; // argv[1]; 
  string outFileName = argv[2]; // "trackerAcceptancePlots.test.root"; // argv[2];

  string treeNameTracks = "trackerAcceptanceTrees/TrackReco";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName .c_str());
  // Get tree
  
  TTree *treeTracks = (TTree*)fin->Get(treeNameTracks.c_str());

  
  cout<<"\nOpened tree:\t"<<treeNameTracks<<" "<<treeTracks<<" from file "<<inFileName<<" "<<fin<<endl;

  // Book output

  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");
   
  // Fill histograms
  Run(treeTracks, fout, quality);

  // Close
  fout->Close();
  fin->Close();

  cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;
   
  return 0;
}
