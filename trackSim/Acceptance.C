// Read ROOT trees 
// Sam Grant

// Plot vertical momentum in slices of y-position

// I really don't understand what this does, is it even used in the final analysis or was this just when i was losing my mind trying random nonsense?

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
double Rmagic = 7112; // mm 

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

void Run(TTree *tree, TFile *output, bool quality = true, bool truth = false) {

  double boostFactor = 5e3*(1/gmagic);
  double momBoostFactor = 1.;

  string stns[] = {"S0S12S18", "S12S18", "S0", "S12", "S18"}; 

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  TH1D *momentumY_[n_stn];
  TH1D *momentumX_[n_stn];
  TH1D *decayY_[n_stn];
  TH1D *decayX_[n_stn];
  TH2D *decayX_vs_decayY_[n_stn];

  vector<TH1D *> momentumY_slices_[n_stn];
  vector<TH1D *> momentumX_slices_[n_stn];
  vector<TH1D *> decayY_slices_[n_stn];
  vector<TH1D *> decayX_slices_[n_stn];
  vector<TH2D *> decayX_vs_decayY_1_slices_[n_stn]; // Vertical
  vector<TH2D *> decayX_vs_decayY_2_slices_[n_stn]; // Vertical

  // Slice vertical position
  int step = 10; // mm 
  int nSlices = 100/step; 

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

      momentumY_[i_stn] = new TH1D((stns[i_stn]+"_MomentumY").c_str(), ";Decay intial momentum Y [MeV];Decays", 700, -70, +70); 
      momentumX_[i_stn] = new TH1D((stns[i_stn]+"_MomentumX").c_str(), ";Decay intial momentum X [MeV];Decays", 700, -70, +70); //int(pmax), -pmax*momBoostFactor, pmax*momBoostFactor); 
      decayY_[i_stn] = new TH1D((stns[i_stn]+"_DecayY").c_str(), ";Decay intial position Y [MeV];Decays", 500, -50, +50);
      decayX_[i_stn] = new TH1D((stns[i_stn]+"_DecayX").c_str(), ";Decay intial position X [MeV];Decays", 500, -50, +50);
      decayX_vs_decayY_[i_stn] = new TH2D((stns[i_stn]+"_DecayX_vs_DecayY").c_str(), ";Decay intial position X [mm];Decay intial Y [mm]", 500, -50, +50, 500, -50, +50);

      // Slice momentum
      for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

        int lo = -50 + i_slice*step; 
        int hi = -50 + step + i_slice*step;

        TH1D *momentumY_slice = new TH1D((stns[i_stn]+"_MomentumY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay intial momentum Y [MeV];Decays", 1000, -70, +70); 
        momentumY_slices_[i_stn].push_back(momentumY_slice);

        TH1D *momentumX_slice = new TH1D((stns[i_stn]+"_MomentumX_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay intial momentum X [MeV];Decays", 1000, -70, +70); // int(pmax), -pmax*momBoostFactor, pmax*momBoostFactor); 
        momentumX_slices_[i_stn].push_back(momentumX_slice);

        TH1D *decayY_slice = new TH1D((stns[i_stn]+"_DecayY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay intial position Y [MeV];Decays", 1000, -100, +100);
        decayY_slices_[i_stn].push_back(decayY_slice);

        TH1D *decayX_slice = new TH1D((stns[i_stn]+"_DecayX_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay intial position X [MeV];Decays", 1000, -100, +100);
        decayX_slices_[i_stn].push_back(decayX_slice);

        TH2D *decayX_vs_decayY_1_slice = new TH2D((stns[i_stn]+"_DecayX_vs_DecayY_1_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay intial position X [mm];Decay intial Y [mm]", 500, -50, +50, 500, -50, +50);
        decayX_vs_decayY_1_slices_[i_stn].push_back(decayX_vs_decayY_1_slice);

        TH2D *decayX_vs_decayY_2_slice = new TH2D((stns[i_stn]+"_DecayX_vs_DecayY_2_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay intial position X [mm];Decay intial Y [mm]", 500, -50, +50, 500, -50, +50);
        decayX_vs_decayY_2_slices_[i_stn].push_back(decayX_vs_decayY_2_slice);

    }

  }

  // Get branches (using header file)
  InitBranches br(tree);

  double targetPerc = 0;
  int64_t nEntries = tree->GetEntries();

  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    double time;
    double x; double y; double z; 
    double px; double py; double pz; 

    // Time is already in us
    if(truth) { 

      x = br.trueVertexPosX; y = br.trueVertexPosY; z = br.trueVertexPosZ;
      px = br.trueVertexMomX; py = br.trueVertexMomY; pz = br.trueVertexMomZ;

    } else { 

      x = br.recoVertexPosX; y = br.recoVertexPosY; z = br.recoVertexPosZ;
      px = br.recoVertexMomX; py = -br.recoVertexMomY; pz = br.recoVertexMomZ;

    }

    bool hitVol = br.hitVolume;
    double pVal = br.pValue;
    bool vertexQual = br.passVertexQuality;

    // Time cuts
    if (quality) {

      // Should include hit vol and pval cuts
      if (!vertexQual) continue;

    }

    int stn = br.station; 

    TVector3 eMom(px, py, pz); 
    TVector3 ePos(x, y, z);

    // RingAngle is angle from x axis, from 0 to 2pi
    double ringAngle = atan2(ePos.z(), ePos.x()); // br.posiInitPosZ, br.posiInitPosX);  
    if (ringAngle < 0) ringAngle += TMath::TwoPi();

    // Rotate into AAR
    eMom.RotateY(ringAngle);
    ePos.RotateY(ringAngle);

    y = ePos.y(); 
    x = ePos.x() - Rmagic; 
    py = eMom.y(); 
    px = eMom.x(); 

    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Fill stations invidually according the station array
    int stn_id = -1;
    if(stn==0) stn_id = 2;
    else if(stn==12) stn_id = 3;
    else if(stn==18) stn_id = 4;
    else cerr<<"Station "<<stn<<" not recognised";

    // Fill histograms
    decayY_[stn_id]->Fill(y);
    decayX_[stn_id]->Fill(x);
    momentumY_[stn_id]->Fill(py);
    momentumX_[stn_id]->Fill(px);
    decayX_vs_decayY_[stn_id]->Fill(x, y);

    // Slice y-position
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = -50 + i_slice*step; 
      int hi = -50 + step + i_slice*step;

      if(y >= double(lo) && y < double(hi)) { 

        decayY_slices_[stn_id].at(i_slice)->Fill(y);
        momentumY_slices_[stn_id].at(i_slice)->Fill(py);
        decayX_vs_decayY_1_slices_[stn_id].at(i_slice)->Fill(x, y);

      }

      if(x >= double(lo) && x < double(hi)) { 

        decayX_slices_[stn_id].at(i_slice)->Fill(x);
        momentumX_slices_[stn_id].at(i_slice)->Fill(px);
        decayX_vs_decayY_2_slices_[stn_id].at(i_slice)->Fill(x, y);

      }

    }

    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }


  }

  // Very ugly code below

  // Combine stations S12&S18
   
  decayY_[1]->Add(decayY_[3], decayY_[4]); 
  decayX_[1]->Add(decayX_[3], decayX_[4]); 
  momentumY_[1]->Add(momentumY_[3], momentumY_[4]); 
  momentumX_[1]->Add(momentumX_[3], momentumX_[4]);  
  decayX_vs_decayY_[1]->Add(decayX_vs_decayY_[3], decayX_vs_decayY_[4]); 

   for ( int i_slice(0); i_slice < nSlices; i_slice++ ) {

    decayY_slices_[1].at(i_slice)->Add(decayY_slices_[3].at(i_slice), decayY_slices_[4].at(i_slice)); 
    decayX_slices_[1].at(i_slice)->Add(decayX_slices_[3].at(i_slice), decayX_slices_[4].at(i_slice)); 
    momentumY_slices_[1].at(i_slice)->Add(momentumY_slices_[3].at(i_slice), momentumY_slices_[4].at(i_slice)); 
    momentumX_slices_[1].at(i_slice)->Add(momentumX_slices_[3].at(i_slice), momentumX_slices_[4].at(i_slice)); 
    decayX_vs_decayY_1_slices_[1].at(i_slice)->Add(decayX_vs_decayY_1_slices_[3].at(i_slice), decayX_vs_decayY_1_slices_[4].at(i_slice)); 
    decayX_vs_decayY_2_slices_[1].at(i_slice)->Add(decayX_vs_decayY_2_slices_[3].at(i_slice), decayX_vs_decayY_2_slices_[4].at(i_slice)); 

  }

  // Combine stations S0&S12&S18

  decayY_[0]->Add(decayY_[1], decayY_[2]); 
  decayX_[0]->Add(decayX_[1], decayX_[2]); 
  momentumY_[0]->Add(momentumY_[1], momentumY_[2]); 
  momentumX_[0]->Add(momentumX_[1], momentumX_[2]);
  decayX_vs_decayY_[0]->Add(decayX_vs_decayY_[1], decayX_vs_decayY_[2]); 

   for ( int i_slice(0); i_slice < nSlices; i_slice++ ) {

    decayY_slices_[0].at(i_slice)->Add(decayY_slices_[1].at(i_slice), decayY_slices_[2].at(i_slice)); 
    decayX_slices_[0].at(i_slice)->Add(decayX_slices_[1].at(i_slice), decayX_slices_[2].at(i_slice)); 
    momentumY_slices_[0].at(i_slice)->Add(momentumY_slices_[1].at(i_slice), momentumY_slices_[2].at(i_slice)); 
    momentumX_slices_[0].at(i_slice)->Add(momentumX_slices_[1].at(i_slice), momentumX_slices_[2].at(i_slice)); 
    decayX_vs_decayY_1_slices_[0].at(i_slice)->Add(decayX_vs_decayY_1_slices_[1].at(i_slice), decayX_vs_decayY_1_slices_[2].at(i_slice)); 
    decayX_vs_decayY_2_slices_[0].at(i_slice)->Add(decayX_vs_decayY_2_slices_[1].at(i_slice), decayX_vs_decayY_2_slices_[2].at(i_slice)); 

  }

  // Write to output
  // Set output directory
  output->mkdir("SanityPlots");
  output->mkdir("VerticalMomentumSlices");  
  output->mkdir("RadialMomentumSlices");  

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    output->cd("SanityPlots");

    decayY_[i_stn]->Write();
    decayX_[i_stn]->Write();
    momentumY_[i_stn]->Write();
    momentumX_[i_stn]->Write();
    decayX_vs_decayY_[i_stn]->Write();

    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) {

      output->cd("VerticalMomentumSlices"); 
      decayY_slices_[i_stn].at(i_slice)->Write();
      momentumY_slices_[i_stn].at(i_slice)->Write();
      decayX_vs_decayY_1_slices_[i_stn].at(i_slice)->Write();

      output->cd("RadialMomentumSlices");
      decayX_slices_[i_stn].at(i_slice)->Write();
      momentumX_slices_[i_stn].at(i_slice)->Write();
      decayX_vs_decayY_2_slices_[i_stn].at(i_slice)->Write();

    }

  }

  cout<<"Written plots."<<endl;

  return;

}

int main(int argc, char *argv[]) {

  bool quality = true;
  bool truth = false;

  string inFileName = argv[1]; 
  string outFileName = argv[2];

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
