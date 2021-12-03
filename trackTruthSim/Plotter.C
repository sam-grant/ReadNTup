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

double ModTime(double time, int nPeriods = 1) {

  double g2fracTime = time / (nPeriods * g2Period);
  int g2fracTimeInt = int(g2fracTime);
  double g2ModTime = (g2fracTime - g2fracTimeInt) * nPeriods * g2Period;

  return g2ModTime;

}

double RandomisedTime(TRandom3 *rand, double time) { 
  return rand->Uniform(time-T_c/2, time+T_c/2);
}

void Run(TTree *tree, TFile *output, bool quality, bool truth) {

  // Set the number of periods for the longer modulo plots
  int moduloMultiple = 4; 

  double boostFactor = 5e3*(1/gmagic);

  string stns[] = {"S0S12S18", "S12S18", "S0", "S12", "S18"}; 

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  TH1D *momentum_[n_stn];
  TH1D *momY_[n_stn];
  TH1D *momX_[n_stn];
  TH1D *momZ_[n_stn];
  TH2D *decayZ_vs_decayX_[n_stn];
  TH1D *wiggle_[n_stn];
  TH1D *wiggle_mod_[n_stn];
  TH1D *wiggle_mod_long_[n_stn];
  TH1D *thetaY_[n_stn];
  TH2D *thetaY_vs_time_[n_stn];
  TH2D *thetaY_vs_time_20ns_[n_stn];
  TH2D *thetaY_vs_time_50ns_[n_stn];
  TH2D *thetaY_vs_time_mod_[n_stn];
  TH2D *thetaY_vs_time_mod_20ns_[n_stn];
  TH2D *thetaY_vs_time_mod_50ns_[n_stn];
  TH2D *thetaY_vs_time_mod_long_[n_stn];
  TH2D *thetaY_vs_time_mod_long_20ns_[n_stn];
  TH2D *thetaY_vs_time_mod_long_50ns_[n_stn];

  // Other scans 
  vector<TH1D*> thetaY_mom_slices_[n_stn];
  vector<TH1D*> Y_mom_slices_[n_stn];
  vector<TH1D*> pY_mom_slices_[n_stn];
  vector<TH1D*> p_mom_slices_[n_stn];

  // Momentum scans of theta_t modulo 
  vector<TH2D*> thetaY_vs_time_mod_slices_[n_stn];
  vector<TH2D*> thetaY_vs_time_mod_20ns_slices_[n_stn];
  vector<TH2D*> thetaY_vs_time_mod_50ns_slices_[n_stn];
  vector<TH2D*> thetaY_vs_time_mod_long_slices_[n_stn];
  vector<TH2D*> thetaY_vs_time_mod_long_20ns_slices_[n_stn];
  vector<TH2D*> thetaY_vs_time_mod_long_50ns_slices_[n_stn];

  // Slice momentum
  int step = 250; 
  int nSlices = pmax/step;

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    momentum_[i_stn] = new TH1D((stns[i_stn]+"_Momentum").c_str(), ";Track momentum [MeV];Tracks", int(pmax), 0, pmax); 
    momY_[i_stn] = new TH1D((stns[i_stn]+"_MomentumY").c_str(), ";Track momentum Y [MeV];Tracks", 1000, -60, 60); 
    momX_[i_stn] = new TH1D((stns[i_stn]+"_MomentumX").c_str(), ";Track momentum X [MeV];Tracks", int(pmax), -pmax, pmax); 
    momZ_[i_stn] = new TH1D((stns[i_stn]+"_MomentumZ").c_str(), ";Track momentum Z [MeV];Tracks", int(pmax), -pmax, pmax); 
    decayZ_vs_decayX_[i_stn] = new TH2D((stns[i_stn]+"_DecayZ_vs_DecayX").c_str(), ";Decay vertex position X [mm];Decay vertex position Z [mm]", 800, -8000, 8000, 800, -8000, 8000);
    wiggle_[i_stn] = new TH1D((stns[i_stn]+"_Wiggle").c_str(), ";Decay time [#mus];Tracks", 2700, 0, 2700*T_c);
    wiggle_mod_[i_stn] = new TH1D((stns[i_stn]+"_Wiggle_Modulo").c_str(), ";t_{g#minus2}^{mod} [#mus];Tracks / 149.2 ns", 29, 0, g2Period); 
    wiggle_mod_long_[i_stn] = new TH1D((stns[i_stn]+"_Wiggle_Modulo_Long").c_str(), ";Time modulo [#mus];Tracks / 149.2 ns", 41*moduloMultiple, 0, g2Period*moduloMultiple); 
    thetaY_[i_stn] = new TH1D((stns[i_stn]+"_ThetaY").c_str(), ";#theta_{y} [mrad];Tracks", 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_time_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time").c_str(), ";Decay time [#mus]; #theta_{y} [mrad] / 149.2 ns ", 2700, 0, 2700*T_c, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_time_20ns_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_20ns").c_str(), ";Decay time [#mus]; #theta_{y} [mrad] / 20 ns ", 20000, 0, 20000*20e-3, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_time_50ns_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_50ns").c_str(), ";Decay time [#mus]; #theta_{y} [mrad] / 50 ns ", 8000, 0, 8000*50e-3, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_time_mod_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo").c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 149.2 ns", 29, 0, g2Period, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_time_mod_50ns_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_50ns").c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 50 ns", 87, 0, g2Period, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_time_mod_20ns_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_20ns").c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 20 ns", 174, 0, g2Period, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_time_mod_long_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_Long").c_str(), ";Time modulo [#mus]; #theta_{y} [mrad] / 149.2 ns", 29*moduloMultiple, 0, g2Period*moduloMultiple, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_time_mod_long_20ns_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_Long_20ns").c_str(), ";Time modulo [#mus]; #theta_{y} [mrad] / 20 ns", 174*moduloMultiple, 0, g2Period*moduloMultiple, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_time_mod_long_50ns_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_Long_50ns").c_str(), ";Time modulo [#mus]; #theta_{y} [mrad] / 50 ns", 87*moduloMultiple, 0, g2Period*moduloMultiple, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
   
    // Slice momentum
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = 0 + i_slice*step; 
      int hi = step + i_slice*step;

      TH1D *h_p_mom_slices = new TH1D((stns[i_stn]+"_Momentum_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Track momentum [MeV];Tracks",  int(pmax), 0, pmax);
      p_mom_slices_[i_stn].push_back(h_p_mom_slices);

      TH1D *h_pY_mom_slices = new TH1D((stns[i_stn]+"_MomentumY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Track momentum Y MeV];Tracks",  1000, -60, 60);
      pY_mom_slices_[i_stn].push_back(h_pY_mom_slices);

      TH1D *h_Y_mom_slice = new TH1D((stns[i_stn]+"_Y_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Vertical decay position [mm];Tracks",  180, -60, 60);
      Y_mom_slices_[i_stn].push_back(h_Y_mom_slice);

      TH2D *h_thetaY_vs_time_mod_slice = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 149.2 ns", 29, 0, g2Period, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
      thetaY_vs_time_mod_slices_[i_stn].push_back(h_thetaY_vs_time_mod_slice);

      TH2D *h_thetaY_vs_time_mod_20ns_slice = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_20ns_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 20 ns", 174, 0, g2Period, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
      thetaY_vs_time_mod_20ns_slices_[i_stn].push_back(h_thetaY_vs_time_mod_20ns_slice);

      TH2D *h_thetaY_vs_time_mod_50ns_slice = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_50ns_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 50 ns", 87, 0, g2Period, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
      thetaY_vs_time_mod_50ns_slices_[i_stn].push_back(h_thetaY_vs_time_mod_50ns_slice);

      TH2D *h_thetaY_vs_time_mod_long_slice = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_Long_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 149.2 ns", 29*moduloMultiple, 0, g2Period*moduloMultiple, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
      thetaY_vs_time_mod_long_slices_[i_stn].push_back(h_thetaY_vs_time_mod_long_slice);

      TH2D *h_thetaY_vs_time_mod_long_20ns_slice = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_Long_20ns_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 20 ns", 174*moduloMultiple, 0, g2Period*moduloMultiple, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
      thetaY_vs_time_mod_long_20ns_slices_[i_stn].push_back(h_thetaY_vs_time_mod_long_20ns_slice);

      TH2D *h_thetaY_vs_time_mod_long_50ns_slice = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_Long_50ns_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 50 ns", 87*moduloMultiple, 0, g2Period*moduloMultiple, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
      thetaY_vs_time_mod_long_50ns_slices_[i_stn].push_back(h_thetaY_vs_time_mod_long_50ns_slice);

      TH1D *h_thetaY_mom_slice = new TH1D((stns[i_stn]+"_ThetaY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";#theta_{y} [mrad];Tracks",  500, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
      thetaY_mom_slices_[i_stn].push_back(h_thetaY_mom_slice);

    }

  }

  // Get branches (using header file)
  InitBranches br(tree);

  double targetPerc = 0;
  int64_t nEntries = tree->GetEntries();

  double muAngleMax = 0;

  int64_t counter = 0;
  
  // For FR randomisation
  TRandom3 *rand = new TRandom3(12345);

  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    double time;
    double x; double y; double z; 
    double px; double py; double pz; 

    // Time is already in us
    if(truth) { 

      time = br.trueTime; // 
      x = br.trueVertexPosX; y = br.trueVertexPosY; z = br.trueVertexPosZ;
      px = br.trueVertexMomX; py = br.trueVertexMomY; pz = br.trueVertexMomZ;

    } else { 

      time = br.recoTime;
      x = br.recoVertexPosX; y = br.recoVertexPosY; z = br.recoVertexPosZ;
      px = br.recoVertexMomX; py = -br.recoVertexMomY; pz = br.recoVertexMomZ;

    }

    time = RandomisedTime(rand, time); // randomise out the FR

    bool hitVol = br.hitVolume;
    double pVal = br.pValue;
    bool vertexQual = br.passVertexQuality;

    // Time cuts
    if (quality) {

      if (time < g2Period*7 || time > g2Period*70) continue; 

      // Should include hit vol and pval cuts
      // if (!vertexQual) continue;

    }

    int stn = br.station; 

    TVector3 eMom(px, py, pz); 
    TVector3 ePos(x, y, z);

    double g2ModTime = ModTime(time);
    double longModTime = ModTime(time, moduloMultiple);

    /////////////////////////////////
    //                             //
    // Define the vertical angle   //
    //                             //
    /////////////////////////////////

    // We need to do this in a way that doesn't involve the z-component of momentum
    // It should be an entirely tranvserse quantity 

    double p = eMom.Mag();

    // If the you have the vertex cut switched on, py must be negative
    double theta_y = asin(py/p);

    theta_y = theta_y * 1e3;

    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Fill stations invidually according the station array
    int stn_id = -1;
    if(stn==0) stn_id = 2;
    else if(stn==12) stn_id = 3;
    else if(stn==18) stn_id = 4;
    else cerr<<"Station "<<stn<<" not recognised";

    decayZ_vs_decayX_[stn_id]->Fill(x, z);

    // g-2 cuts. See Fienberg thesis figure 2.10
    if(p > 1900  && p < pmax) {
      wiggle_[stn_id]->Fill(time);
      wiggle_mod_[stn_id]->Fill(g2ModTime);
      if(time>8*g2Period) wiggle_mod_long_[stn_id]->Fill(longModTime);
    } 

    momY_[stn_id]->Fill(py);
    momX_[stn_id]->Fill(px);
    momZ_[stn_id]->Fill(pz);

    // EDM cuts
    if( quality && p > pLo && p < pHi) { 

      momentum_[stn_id]->Fill(p);
      thetaY_[stn_id]->Fill(theta_y);
      thetaY_vs_time_[stn_id]->Fill(time, theta_y);
      thetaY_vs_time_20ns_[stn_id]->Fill(time, theta_y);
      thetaY_vs_time_50ns_[stn_id]->Fill(time, theta_y);
      thetaY_vs_time_mod_[stn_id]->Fill(g2ModTime, theta_y);
      thetaY_vs_time_mod_50ns_[stn_id]->Fill(g2ModTime, theta_y);

      if(time > 8*g2Period) {
        thetaY_vs_time_mod_long_[stn_id]->Fill(longModTime, theta_y);
        thetaY_vs_time_mod_long_20ns_[stn_id]->Fill(longModTime, theta_y); 
        thetaY_vs_time_mod_long_50ns_[stn_id]->Fill(longModTime, theta_y);    
      }  

    } else if(!quality) { 

      momentum_[stn_id]->Fill(p);
      thetaY_[stn_id]->Fill(theta_y);
      thetaY_vs_time_[stn_id]->Fill(time, theta_y);
      thetaY_vs_time_20ns_[stn_id]->Fill(time, theta_y);
      thetaY_vs_time_50ns_[stn_id]->Fill(time, theta_y);
      thetaY_vs_time_mod_[stn_id]->Fill(g2ModTime, theta_y);
      thetaY_vs_time_mod_20ns_[stn_id]->Fill(g2ModTime, theta_y);
      thetaY_vs_time_mod_50ns_[stn_id]->Fill(g2ModTime, theta_y);

      if(time > 8*g2Period) {
        thetaY_vs_time_mod_long_[stn_id]->Fill(longModTime, theta_y);
        thetaY_vs_time_mod_long_20ns_[stn_id]->Fill(longModTime, theta_y);
        thetaY_vs_time_mod_long_50ns_[stn_id]->Fill(longModTime, theta_y);    
      }  

    }

    // Slice momentum 
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = 0 + i_slice*step; 
      int hi = step + i_slice*step;

      if(p >= double(lo) && p < double(hi)) { 

        thetaY_vs_time_mod_slices_[stn_id].at(i_slice)->Fill(g2ModTime, theta_y);
        thetaY_vs_time_mod_20ns_slices_[stn_id].at(i_slice)->Fill(g2ModTime, theta_y);
        thetaY_vs_time_mod_50ns_slices_[stn_id].at(i_slice)->Fill(g2ModTime, theta_y);
      
        if(time > 8*g2Period) {
          thetaY_vs_time_mod_long_slices_[stn_id].at(i_slice)->Fill(longModTime, theta_y);
          thetaY_vs_time_mod_long_20ns_slices_[stn_id].at(i_slice)->Fill(longModTime, theta_y);
          thetaY_vs_time_mod_long_50ns_slices_[stn_id].at(i_slice)->Fill(longModTime, theta_y);    
        }  

        // Other scans 
        thetaY_mom_slices_[stn_id].at(i_slice)->Fill(theta_y);
        Y_mom_slices_[stn_id].at(i_slice)->Fill(y);
        pY_mom_slices_[stn_id].at(i_slice)->Fill(py);
        p_mom_slices_[stn_id].at(i_slice)->Fill(p);

      }

    }

    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }

  }

  // Very ugly code below

  // Combine stations S12&S18
  momentum_[1]->Add(momentum_[3], momentum_[4]);
  wiggle_[1]->Add(wiggle_[3], wiggle_[4]);
  wiggle_mod_[1]->Add(wiggle_mod_[3], wiggle_mod_[4]);
  wiggle_mod_long_[1]->Add(wiggle_mod_long_[3], wiggle_mod_long_[4]);
  thetaY_[1]->Add(thetaY_[3], thetaY_[4]);
  thetaY_vs_time_[1]->Add(thetaY_vs_time_[3], thetaY_vs_time_[4]);
  thetaY_vs_time_20ns_[1]->Add(thetaY_vs_time_20ns_[3], thetaY_vs_time_20ns_[4]);
  thetaY_vs_time_50ns_[1]->Add(thetaY_vs_time_50ns_[3], thetaY_vs_time_50ns_[4]);
  thetaY_vs_time_mod_[1]->Add(thetaY_vs_time_mod_[3], thetaY_vs_time_mod_[4]);
  thetaY_vs_time_mod_20ns_[1]->Add(thetaY_vs_time_mod_20ns_[3], thetaY_vs_time_mod_20ns_[4]);
  thetaY_vs_time_mod_50ns_[1]->Add(thetaY_vs_time_mod_50ns_[3], thetaY_vs_time_mod_50ns_[4]);
  thetaY_vs_time_mod_long_[1]->Add(thetaY_vs_time_mod_long_[3], thetaY_vs_time_mod_long_[4]);
  thetaY_vs_time_mod_long_20ns_[1]->Add(thetaY_vs_time_mod_long_20ns_[3], thetaY_vs_time_mod_long_20ns_[4]);
  thetaY_vs_time_mod_long_50ns_[1]->Add(thetaY_vs_time_mod_long_50ns_[3], thetaY_vs_time_mod_long_50ns_[4]);
  decayZ_vs_decayX_[1]->Add(decayZ_vs_decayX_[3], decayZ_vs_decayX_[4]);
  momX_[1]->Add(momX_[3], momX_[4]);
  momY_[1]->Add(momY_[3], momY_[4]);
  momZ_[1]->Add(momZ_[3], momZ_[4]);

   for ( int i_slice(0); i_slice < nSlices; i_slice++ ) {

      thetaY_vs_time_mod_slices_[1].at(i_slice)->Add(thetaY_vs_time_mod_slices_[3].at(i_slice), thetaY_vs_time_mod_slices_[4].at(i_slice));
      thetaY_vs_time_mod_20ns_slices_[1].at(i_slice)->Add(thetaY_vs_time_mod_20ns_slices_[3].at(i_slice), thetaY_vs_time_mod_20ns_slices_[4].at(i_slice));
      thetaY_vs_time_mod_50ns_slices_[1].at(i_slice)->Add(thetaY_vs_time_mod_50ns_slices_[3].at(i_slice), thetaY_vs_time_mod_50ns_slices_[4].at(i_slice));
      thetaY_vs_time_mod_long_slices_[1].at(i_slice)->Add(thetaY_vs_time_mod_long_slices_[3].at(i_slice), thetaY_vs_time_mod_long_slices_[4].at(i_slice));
      thetaY_vs_time_mod_long_20ns_slices_[1].at(i_slice)->Add(thetaY_vs_time_mod_long_20ns_slices_[3].at(i_slice), thetaY_vs_time_mod_long_20ns_slices_[4].at(i_slice));
      thetaY_vs_time_mod_long_50ns_slices_[1].at(i_slice)->Add(thetaY_vs_time_mod_long_50ns_slices_[3].at(i_slice), thetaY_vs_time_mod_long_50ns_slices_[4].at(i_slice));
      thetaY_mom_slices_[1].at(i_slice)->Add(thetaY_mom_slices_[3].at(i_slice), thetaY_mom_slices_[4].at(i_slice));
      Y_mom_slices_[1].at(i_slice)->Add(Y_mom_slices_[3].at(i_slice), Y_mom_slices_[4].at(i_slice));
      pY_mom_slices_[1].at(i_slice)->Add(pY_mom_slices_[3].at(i_slice), pY_mom_slices_[4].at(i_slice));
      p_mom_slices_[1].at(i_slice)->Add(p_mom_slices_[3].at(i_slice), p_mom_slices_[4].at(i_slice));
      
  }

  // Combine stations S0&S12&S18
  momentum_[0]->Add(momentum_[1], momentum_[2]);
  wiggle_[0]->Add(wiggle_[1], wiggle_[2]);
  wiggle_mod_[0]->Add(wiggle_mod_[1], wiggle_mod_[2]);
  wiggle_mod_long_[0]->Add(wiggle_mod_long_[1], wiggle_mod_long_[2]);
  thetaY_[0]->Add(thetaY_[1], thetaY_[2]);
  thetaY_vs_time_[0]->Add(thetaY_vs_time_[1], thetaY_vs_time_[2]);
  thetaY_vs_time_20ns_[0]->Add(thetaY_vs_time_20ns_[1], thetaY_vs_time_20ns_[2]);
  thetaY_vs_time_50ns_[0]->Add(thetaY_vs_time_50ns_[1], thetaY_vs_time_50ns_[2]);
  thetaY_vs_time_mod_[0]->Add(thetaY_vs_time_mod_[1], thetaY_vs_time_mod_[2]);
  thetaY_vs_time_mod_20ns_[0]->Add(thetaY_vs_time_mod_20ns_[1], thetaY_vs_time_mod_20ns_[2]);
  thetaY_vs_time_mod_50ns_[0]->Add(thetaY_vs_time_mod_50ns_[1], thetaY_vs_time_mod_50ns_[2]);
  thetaY_vs_time_mod_long_[0]->Add(thetaY_vs_time_mod_long_[1], thetaY_vs_time_mod_long_[2]);
  thetaY_vs_time_mod_long_20ns_[0]->Add(thetaY_vs_time_mod_long_20ns_[1], thetaY_vs_time_mod_long_20ns_[2]);
  thetaY_vs_time_mod_long_50ns_[0]->Add(thetaY_vs_time_mod_long_50ns_[1], thetaY_vs_time_mod_long_50ns_[2]);
  decayZ_vs_decayX_[0]->Add(decayZ_vs_decayX_[1], decayZ_vs_decayX_[2]);
  momX_[0]->Add(momX_[1], momX_[2]);
  momY_[0]->Add(momY_[1], momY_[2]);
  momZ_[0]->Add(momZ_[1], momZ_[2]);

   for ( int i_slice(0); i_slice < nSlices; i_slice++ ) {

      thetaY_vs_time_mod_slices_[0].at(i_slice)->Add(thetaY_vs_time_mod_slices_[1].at(i_slice), thetaY_vs_time_mod_slices_[2].at(i_slice));
      thetaY_vs_time_mod_20ns_slices_[0].at(i_slice)->Add(thetaY_vs_time_mod_20ns_slices_[1].at(i_slice), thetaY_vs_time_mod_20ns_slices_[2].at(i_slice));
      thetaY_vs_time_mod_50ns_slices_[0].at(i_slice)->Add(thetaY_vs_time_mod_50ns_slices_[1].at(i_slice), thetaY_vs_time_mod_50ns_slices_[2].at(i_slice));
      thetaY_vs_time_mod_long_slices_[0].at(i_slice)->Add(thetaY_vs_time_mod_long_slices_[1].at(i_slice), thetaY_vs_time_mod_long_slices_[2].at(i_slice));
      thetaY_vs_time_mod_long_20ns_slices_[0].at(i_slice)->Add(thetaY_vs_time_mod_long_20ns_slices_[1].at(i_slice), thetaY_vs_time_mod_long_20ns_slices_[2].at(i_slice));
      thetaY_vs_time_mod_long_50ns_slices_[0].at(i_slice)->Add(thetaY_vs_time_mod_long_50ns_slices_[1].at(i_slice), thetaY_vs_time_mod_long_50ns_slices_[2].at(i_slice));
      thetaY_mom_slices_[0].at(i_slice)->Add(thetaY_mom_slices_[1].at(i_slice), thetaY_mom_slices_[2].at(i_slice));
      Y_mom_slices_[0].at(i_slice)->Add(Y_mom_slices_[1].at(i_slice), Y_mom_slices_[2].at(i_slice));
      pY_mom_slices_[0].at(i_slice)->Add(pY_mom_slices_[1].at(i_slice), pY_mom_slices_[2].at(i_slice));
      p_mom_slices_[0].at(i_slice)->Add(p_mom_slices_[1].at(i_slice), p_mom_slices_[2].at(i_slice));
      
  }

  // Write to output
  // Set output directory
  output->mkdir("SimultaneousAnalysis"); output->mkdir("MomentumBinnedAnalysis");  

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    output->cd("SimultaneousAnalysis");

    momentum_[i_stn]->Write();
    wiggle_[i_stn]->Write();
    wiggle_mod_[i_stn]->Write();
    wiggle_mod_long_[i_stn]->Write();
    thetaY_[i_stn]->Write();
    thetaY_vs_time_[i_stn]->Write();
    thetaY_vs_time_20ns_[i_stn]->Write(); 
    thetaY_vs_time_50ns_[i_stn]->Write(); 
    thetaY_vs_time_mod_[i_stn]->Write();
    thetaY_vs_time_mod_20ns_[i_stn]->Write();
    thetaY_vs_time_mod_50ns_[i_stn]->Write();
    thetaY_vs_time_mod_long_[i_stn]->Write();
    thetaY_vs_time_mod_long_20ns_[i_stn]->Write();
    thetaY_vs_time_mod_long_50ns_[i_stn]->Write();
    decayZ_vs_decayX_[i_stn]->Write();
    momX_[i_stn]->Write();
    momY_[i_stn]->Write();
    momZ_[i_stn]->Write();

   for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) {

      output->cd("MomentumBinnedAnalysis"); 

      thetaY_vs_time_mod_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_time_mod_20ns_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_time_mod_50ns_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_time_mod_long_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_time_mod_long_20ns_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_time_mod_long_50ns_slices_[i_stn].at(i_slice)->Write();
      thetaY_mom_slices_[i_stn].at(i_slice)->Write();
      Y_mom_slices_[i_stn].at(i_slice)->Write();
      pY_mom_slices_[i_stn].at(i_slice)->Write();
      p_mom_slices_[i_stn].at(i_slice)->Write();
      
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
