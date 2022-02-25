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

  string stns[] = {"S12S18", "S12", "S18"}; 

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  int step = 250; 
  int nSlices = pmax/step;

  // Simpler binning
  TH1D *Y_allDecays = new TH1D("Y", ";Decay y-position [mm]; Decays", 20, -50, 50);
  TH1D *thetaY_allDecays = new TH1D("ThetaY", "#theta_{y} [mrad]; Decays", 28, -80, 80); 
  // More representative binning
  TH2D *thetaY_vs_Y_allDecays = new TH2D("ThetaY_vs_Y", ";Decay y-position [mm]; #theta_{y} [mrad]", 48, -60, 60, 48, -TMath::Pi()*gmagic, TMath::Pi()*gmagic); 
 
  TH1D *Y_tracks_[n_stn]; 
  TH1D *thetaY_tracks_[n_stn]; 
  TH2D *thetaY_vs_Y_tracks_[n_stn]; 

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    Y_tracks_[i_stn] = new TH1D((stns[i_stn]+"_Y").c_str(), ";Decay y-position [mm];Tracks", 20, -50, 50);
    thetaY_tracks_[i_stn] = new TH1D((stns[i_stn]+"_ThetaY").c_str(), ";#theta_{y} [mrad];Tracks", 28, -80, 80); 
    thetaY_vs_Y_tracks_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Y").c_str(), ";Decay y-position [mm]; #theta_{y} [mrad]", 48, -60, 60, 48, -TMath::Pi()*gmagic, TMath::Pi()*gmagic); 

  }

  // mom-slices
  vector<TH1D*> mom_slices_;//[i_stn] = new TH1D((stns[i_stn]+"_Momentum").c_str(), ";Track momentum [MeV];Tracks", int(pmax), 0, pmax); 
  vector<TH1D*> Y_allDecays_slices_;
  vector<TH1D*> thetaY_allDecays_slices_;
  vector<TH2D*> thetaY_vs_Y_allDecays_slices_;

  vector<TH1D*> mom_slices_;
  vector<TH1D*> Y_tracks_slices_[n_stn];
  vector<TH1D*> thetaY_tracks_slices_[n_stn];
  vector<TH2D*> thetaY_vs_Y_tracks_slices_[n_stn];


  for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

    int lo = i_slice*step; 
    int hi = step + lo;

    std::string stepStr = to_string(lo)+"_"+to_string(hi);

    Y_allDecays_slices_.push_back(new TH1D(("Y_"+stepStr).c_str(), ";Decay y-position [mm];Tracks", 20, -50, 50));
    thetaY_allDecays_slices_.push_back(new TH1D(("ThetaY_"+stepStr).c_str(), ";#theta_{y} [mrad];Tracks", 28, -80, 80)); 
    thetaY_vs_Y_allDecays_slices_.push_back(new TH2D(("ThetaY_vs_Y_"+stepStr).c_str(), ";Decay y-position [mm]; #theta_{y} [mrad]", 48, -60, 60, 48, -TMath::Pi()*gmagic, TMath::Pi()*gmagic)); 

    for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

      Y_tracks_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_Y_"+stepStr).c_str(), ";Decay y-position [mm];Tracks", 20, -50, 50));
      thetaY_tracks_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_ThetaY_"+stepStr).c_str(), ";#theta_{y} [mrad];Tracks", 28, -80, 80)); 
      thetaY_vs_Y_tracks_slices_[i_stn].push_back(new TH2D((stns[i_stn]+"_ThetaY_vs_Y_"+stepStr).c_str(), ";Decay y-position [mm]; #theta_{y} [mrad]", 48, -60, 60, 48, -TMath::Pi()*gmagic, TMath::Pi()*gmagic)); 

    }

  }


  // Get branches (using header file)
  InitBranchesAllDecays allDecays(treeAllDecays);

  cout<<"-----> Filling decay histograms"<<endl;

  for(int64_t entry = 0; entry < treeAllDecays->GetEntries(); entry++) {

    treeAllDecays->GetEntry(entry);

    double y = allDecays.posiInitPosY;
    double p = allDecays.posiInitP;
    double py = allDecays.posiInitPY;
    double theta_y = asin(py/p) * 1e3; 

    //if(quality && (p < pLo || p > pHi)) continue;

    thetaY_allDecays->Fill(theta_y);
    Y_allDecays->Fill(y);
    thetaY_vs_Y_allDecays->Fill(y, theta_y);

    // Slice y
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = i_slice*step; 
      int hi = step + lo;

      if(p >= double(lo) && p < double(hi)) { 

        thetaY_allDecays_slices_.at(i_slice)->Fill(theta_y);
        Y_allDecays_slices_.at(i_slice)->Fill(y);
        thetaY_vs_Y_allDecays_slices_.at(i_slice)->Fill(y, theta_y);
      
      }

    }

  }

  InitBranchesTracks tracks(treeTracks);

  cout<<"-----> Filling track histograms"<<endl;

  for(int64_t entry = 0; entry < treeAllDecays->GetEntries(); entry++) {

    treeTracks->GetEntry(entry);

    double y = tracks.decayVertexPosY;
    double p = tracks.decayVertexMom;
    double py = -tracks.decayVertexMomY;
    double theta_y = asin(py/p) * 1e3; 
    int stn = tracks.station;

    if(quality) {
      if(!tracks.passVertexQuality) continue;
      if(tracks.pValue < 0.05) continue;
      if(tracks.hitVolume) continue; 
      if(p < pLo || p > pHi) continue;
    } 

    int stn_id = -1;
    if(stn==12) stn_id = 1;
    else if(stn==18) stn_id = 2;
    else {
      // cerr<<"Station "<<stn<<" not recognised";
      cout<<"WARNING: Station "<<stn<<" not recognised"<<endl;
    }

    thetaY_tracks_[stn_id]->Fill(theta_y);
    Y_tracks_[stn_id]->Fill(y);
    thetaY_vs_Y_tracks_[stn_id]->Fill(y, theta_y);


    // Slice y
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = -45 + i_slice*step; 
      int hi = step + lo;

      if(y >= double(lo) && y < double(hi)) { 

        thetaY_tracks_slices_[stn_id].at(i_slice)->Fill(theta_y);
        Y_tracks_slices_[stn_id].at(i_slice)->Fill(y);
        thetaY_vs_Y_tracks_slices_[stn_id].at(i_slice)->Fill(y, theta_y);
      
      }

    }

  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

  // Write to output
  // Set output directory
  output->mkdir("AllDecays"); output->cd("AllDecays");

  cout<<"-----> Writing decay histograms"<<endl;

  Y_allDecays->Write();
  thetaY_allDecays->Write();
  thetaY_vs_Y_allDecays->Write();

  output->mkdir("AllDecays/VertPosBins"); output->cd("AllDecays/VertPosBins");

  for ( int i_slice(0); i_slice < nSlices; i_slice++ ) {

    Y_allDecays_slices_.at(i_slice)->Write(); // Add(Y_allDecays_slices_.at(i_slice), Y_allDecays_slices_.at(i_slice));
    thetaY_allDecays_slices_.at(i_slice)->Write();//Add(thetaY_allDecays_slices_.at(i_slice), thetaY_allDecays_slices_.at(i_slice));
    thetaY_vs_Y_allDecays_slices_.at(i_slice)->Write();//Add(thetaY_vs_Y_allDecays_slices_.at(i_slice), thetaY_vs_Y_allDecays_slices_.at(i_slice));
      
  }



  // Combine stations S12&S18
  cout<<"-----> Combining stations"<<endl;
  Y_tracks_[0]->Add(Y_tracks_[1], Y_tracks_[2]);
  thetaY_tracks_[0]->Add(thetaY_tracks_[1], thetaY_tracks_[2]);
  thetaY_vs_Y_tracks_[0]->Add(thetaY_vs_Y_tracks_[1], thetaY_vs_Y_tracks_[2]);

   for ( int i_slice(0); i_slice < nSlices; i_slice++ ) {

    Y_tracks_slices_[0].at(i_slice)->Add(Y_tracks_slices_[1].at(i_slice), Y_tracks_slices_[2].at(i_slice));
    thetaY_tracks_slices_[0].at(i_slice)->Add(thetaY_tracks_slices_[1].at(i_slice), thetaY_tracks_slices_[2].at(i_slice));
    thetaY_vs_Y_tracks_slices_[0].at(i_slice)->Add(thetaY_vs_Y_tracks_slices_[1].at(i_slice), thetaY_vs_Y_tracks_slices_[2].at(i_slice));
      
  }

  cout<<"-----> Writing track histograms"<<endl;

  output->mkdir("Tracks"); output->cd("Tracks");
  output->mkdir("Tracks/VertPosBins"); 

  for (int i_stn = 0; i_stn < n_stn; i_stn++) {

    Y_tracks_[i_stn]->Write();
    thetaY_tracks_[i_stn]->Write();
    thetaY_vs_Y_tracks_[i_stn]->Write();

     for ( int i_slice(0); i_slice < nSlices; i_slice++ ) {

      output->cd("Tracks/VertPosBins");

      Y_tracks_slices_[i_stn].at(i_slice)->Write();//Add(Y_tracks_slices_[i_stn].at(i_slice), Y_tracks_slices_[i_stn].at(i_slice));
      thetaY_tracks_slices_[i_stn].at(i_slice)->Write();//Add(thetaY_tracks_slices_[i_stn].at(i_slice), thetaY_tracks_slices_[i_stn].at(i_slice));
      thetaY_vs_Y_tracks_slices_[i_stn].at(i_slice)->Write();//Add(thetaY_vs_Y_tracks_slices_[i_stn].at(i_slice), thetaY_vs_Y_tracks_slices_[i_stn].at(i_slice));
      
    }

  }


  cout<<"-----> Written plots."<<endl;

  return;

}

int main(int argc, char *argv[]) {

  bool quality = true;

  string inFileName = argv[1]; // "trackerAcceptanceTrees.test.root"; // argv[1]; 
  string outFileName = argv[2]; // "trackerAcceptancePlots.test.root"; // argv[2];

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
