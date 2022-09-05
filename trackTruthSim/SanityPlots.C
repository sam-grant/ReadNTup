// Produce sanity plots from MC trees
// Position histograms, vertical angle histograms, momentum histograms, vertical angle verus momentum, vertical angle vs position
// Could also add t_mod plots? 

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

double omega_a = 1.439311;// (average for mu+ at BNL) 0.00143934*1e3;//1.439311; // rad/us 0.00143934; // kHz from gm2const, it's an angular frequency though...
double g2Period = TMath::TwoPi() / omega_a;//s * 1e-3; // us
double mMu = 105.6583715; // MeV
double aMu = 11659208.9e-10; 
double gmagic = std::sqrt( 1.+1./aMu );
double pmax = 1.01 * mMu * gmagic;
double eMass = 0.510999;
double T_c = 149.2 * 1e-3; // cyclotron period [us]
double R_magic = 7112; // mm 

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

  string stns[] = {"S0S12S18", "S12S18", "S0", "S12", "S18"}; 

  cout<<"---> Booking histograms"<<endl;

  // For direct comp. with decays
  double boostFactor = 5e3*(1/gmagic);

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  TH1D *p_[n_stn];
  TH1D *pX_[n_stn];
  TH1D *pY_[n_stn];
  TH1D *pZ_[n_stn];
  TH1D *R_[n_stn];
  TH1D *Y_[n_stn];
  TH1D *Phi_[n_stn];
  TH2D *decayX_vs_decayZ_[n_stn];
  TH2D *decayY_vs_decayR_[n_stn];
  TH1D *thetaY_[n_stn];
  TH2D *thetaY_vs_t_[n_stn];
  TH2D *thetaY_vs_p_[n_stn];
  TH2D *thetaY_vs_Y_[n_stn];
  TH2D *thetaY_vs_R_[n_stn];
  TH2D *thetaY_vs_Phi_[n_stn];

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    p_[i_stn] = new TH1D((stns[i_stn]+"_Momentum").c_str(), ";Momentum [MeV];Tracks", int(pmax), 0, pmax); 
    pY_[i_stn] = new TH1D((stns[i_stn]+"_MomentumY").c_str(), ";Momentum Y [MeV];Tracks", 1000, -60, 60); 
    pX_[i_stn] = new TH1D((stns[i_stn]+"_MomentumX").c_str(), ";Momentum X [MeV];Tracks", int(pmax), -pmax, pmax); 
    pZ_[i_stn] = new TH1D((stns[i_stn]+"_MomentumZ").c_str(), ";Momentum Z [MeV];Tracks", int(pmax), -pmax, pmax); 
    R_[i_stn] = new TH1D((stns[i_stn]+"_R").c_str(), ";R [mm];Tracks", 120, -60, 60); 
    Y_[i_stn] = new TH1D((stns[i_stn]+"_Y").c_str(), ";Y [mm];Tracks", 120, -60, 60);
    Phi_[i_stn] = new TH1D((stns[i_stn]+"_Phi").c_str(), ";Ring azimuthal angle [rad];Tracks", 600, 0, TMath::TwoPi()); 
    decayX_vs_decayZ_[i_stn] = new TH2D((stns[i_stn]+"_DecayX_vs_DecayZ").c_str(), ";Decay vertex position Z [mm];Decay vertex position X [mm]", 800, -8000, 8000, 800, -8000, 8000);
    decayY_vs_decayR_[i_stn] = new TH2D((stns[i_stn]+"_DecayY_vs_DecayR").c_str(), ";Decay vertex position R [mm];Decay vertex position Y [mm]", 120, -60, 60, 120, -60, 60);
    thetaY_[i_stn] = new TH1D((stns[i_stn]+"_ThetaY").c_str(), ";#theta_{y} [mrad];Tracks", 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic); // a bit approximate 
    thetaY_vs_t_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time").c_str(), ";Decay time [#mus]; #theta_{y} [mrad] / 149.2 ns ", 2700, 0, 2700*T_c, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_p_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Momentum").c_str(), ";Momentum [MeV]; #theta_{y} [mrad] / 10 MeV ", 300, 0, 3000, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);//int(pmax), 0, pmax, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_R_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_R").c_str(), ";R [mm]; #theta_{y} [mrad] / MeV ", 120, -60, 60, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_Y_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Y").c_str(), ";Y [mm]; #theta_{y} [mrad] / MeV ", 120, -60, 60, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_Phi_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Phi").c_str(), ";Ring azimuthal angle [rad]; #theta_{y} [mrad] / MeV ", 600, 0, TMath::TwoPi(), 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    
  } // stn loop

  // Get branches (using header file)
  InitBranches br(tree);

  double targetPerc = 0;
  int64_t nEntries = tree->GetEntries();

  double muAngleMax = 0;
  
  cout<<"---> Event loop"<<endl;

  for(int64_t entry = 0; entry < nEntries; entry++) {

  	// cout<<"---> Getting entry"<<endl;

    tree->GetEntry(entry);

    // Quality variables
	  bool hitVol = br.hitVolume;
    double pVal = br.pValue;
    bool vertexQual = br.passVertexQuality;

    // Time cuts
    if (quality) {
     if(!vertexQual) continue; 
     if(hitVol) continue; // just to be sure
     if(pVal < 0.05) continue;
    }

    // cout<<"---> Getting variables"<<endl;

    // Physics variables

    double t;
    double x; double y; double z; double r; double phi;
    double px; double py; double pz; 

    if(truth) { 

      t = br.trueTime; // 
      x = br.trueVertexPosX; y = br.trueVertexPosY; z = br.trueVertexPosZ;
      px = br.trueVertexMomX; py = br.trueVertexMomY; pz = br.trueVertexMomZ;
      r = sqrt(pow(br.trueVertexPosX,2)+pow(br.trueVertexPosZ,2)) - R_magic;
      phi = atan2(br.trueVertexPosZ, br.trueVertexPosX);

    } else { 

      t = br.recoTime;
      x = br.recoVertexPosX; y = br.recoVertexPosY; z = br.recoVertexPosZ;
      px = br.recoVertexMomX; py = -br.recoVertexMomY; pz = br.recoVertexMomZ;
      r = sqrt(pow(br.recoVertexPosX,2)+pow(br.recoVertexPosZ,2)) - R_magic;
      phi = atan2(br.recoVertexPosZ, br.recoVertexPosX);

    }

    // Needed?
    TVector3 eMom(px, py, pz); 
    TVector3 ePos(x, y, z);

    if (phi  < 0) phi += TMath::TwoPi();

    double theta_y = asin(py/eMom.Mag()) * 1e3; // mrad

    // Stations
    int stn = br.station; 

    // ~1.5% of time, truth vertex just reports zero momentum!? 
    if(eMom.Mag() == 0) continue;

    // Time cut
    if(t < g2Period*7) continue;

    // Fill stations invidually according the station array
    int stn_id = -1;
    if(stn==0) stn_id = 2;
    else if(stn==12) stn_id = 3;
    else if(stn==18) stn_id = 4;
    else cerr<<"Station "<<stn<<" not recognised";

    // cout<<"---> Filling histograms"<<endl;

    // Fill histograms 
    p_[stn_id]->Fill(eMom.Mag());
    pX_[stn_id]->Fill(px);
    pY_[stn_id]->Fill(py);
    pZ_[stn_id]->Fill(pz);
    R_[stn_id]->Fill(r);
    Y_[stn_id]->Fill(y);
    Phi_[stn_id]->Fill(phi);
    decayX_vs_decayZ_[stn_id]->Fill(z, x);
    decayY_vs_decayR_[stn_id]->Fill(r, y);
    thetaY_[stn_id]->Fill(theta_y);
    thetaY_vs_t_[stn_id]->Fill(t, theta_y);
    thetaY_vs_p_[stn_id]->Fill(eMom.Mag(), theta_y);
    thetaY_vs_R_[stn_id]->Fill(r, theta_y);
    thetaY_vs_Y_[stn_id]->Fill(y, theta_y);
    thetaY_vs_Phi_[stn_id]->Fill(phi, theta_y);

/*    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }*/

  }

  cout<<"---> Merging stations"<<endl;
  // Combine stations, presumably there is a smarter way of doing this

  // S12&S18
  p_[1]->Add(p_[3], p_[4]);
  pX_[1]->Add(pX_[3], pX_[4]);
  pY_[1]->Add(pY_[3], pY_[4]);
  pZ_[1]->Add(pZ_[3], pZ_[4]);
  R_[1]->Add(R_[3], R_[4]);
  Phi_[1]->Add(Phi_[3], Phi_[4]);
  Y_[1]->Add(Y_[3], Y_[4]);
  decayX_vs_decayZ_[1]->Add(decayX_vs_decayZ_[3], decayX_vs_decayZ_[4]);
  decayY_vs_decayR_[1]->Add(decayY_vs_decayR_[3], decayY_vs_decayR_[4]);
  thetaY_[1]->Add(thetaY_[3], thetaY_[4]);
  thetaY_vs_t_[1]->Add(thetaY_vs_t_[3], thetaY_vs_t_[4]);
  thetaY_vs_p_[1]->Add(thetaY_vs_p_[3], thetaY_vs_p_[4]);
  thetaY_vs_R_[1]->Add(thetaY_vs_R_[3], thetaY_vs_R_[4]);
  thetaY_vs_Phi_[1]->Add(thetaY_vs_Phi_[3], thetaY_vs_Phi_[4]);

  // S0&S12&S18
  p_[0]->Add(p_[1], p_[2]);
  pX_[0]->Add(pX_[1], pX_[2]);
  pY_[0]->Add(pY_[1], pY_[2]);
  pZ_[0]->Add(pZ_[1], pZ_[2]);
  R_[0]->Add(R_[1], R_[2]);
  Phi_[0]->Add(Phi_[1], Phi_[2]);
  Y_[0]->Add(Y_[1], Y_[2]);
  decayX_vs_decayZ_[0]->Add(decayX_vs_decayZ_[1], decayX_vs_decayZ_[2]);
  decayY_vs_decayR_[0]->Add(decayY_vs_decayR_[1], decayY_vs_decayR_[2]);
  thetaY_[0]->Add(thetaY_[1], thetaY_[2]);
  thetaY_vs_t_[0]->Add(thetaY_vs_t_[1], thetaY_vs_t_[2]);
  thetaY_vs_p_[0]->Add(thetaY_vs_p_[1], thetaY_vs_p_[2]);
  thetaY_vs_R_[0]->Add(thetaY_vs_R_[1], thetaY_vs_R_[2]);
  thetaY_vs_Phi_[0]->Add(thetaY_vs_Phi_[1], thetaY_vs_Phi_[2]);

  // Write to output

  cout<<"---> Writing to output"<<endl;

  // Set output directory
  output->mkdir("SanityPlots"); output->cd("SanityPlots");

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

	p_[i_stn]->Write();
	pX_[i_stn]->Write();
	pY_[i_stn]->Write();
	pZ_[i_stn]->Write();
	R_[i_stn]->Write();
	Phi_[i_stn]->Write();
  Y_[i_stn]->Write();
	decayX_vs_decayZ_[i_stn]->Write();
	decayY_vs_decayR_[i_stn]->Write();
	thetaY_[i_stn]->Write();
	thetaY_vs_t_[i_stn]->Write();
	thetaY_vs_p_[i_stn]->Write();
	thetaY_vs_R_[i_stn]->Write();
  thetaY_vs_Y_[i_stn]->Write();
	thetaY_vs_Phi_[i_stn]->Write();

  }

  cout<<"---> Done"<<endl;

  return;

}

int main(int argc, char *argv[]) {

  bool quality = true;
  bool truth = true; // default is true

  string inFileName = argv[1]; 
  string outFileName = argv[2];
  string truthString = argv[3];

  if(truthString=="truth") truth = true;
  else truth = false;

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