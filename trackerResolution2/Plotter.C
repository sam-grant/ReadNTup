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
double R_magic = 7112;

// I've had to introduce some 
/*double pLo = 750; 
double pHi = 2750;*/

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

void Run(TTree *treeTracks, TFile *output, bool quality, double pLo = 1000, double pHi = 2500) {

  // Position, momentum and vertical angle resolutions over a momentum range and as a function of momentum 
  TH1D *y_res = new TH1D("y_res", ";#Delta vertical decay position [mm];Vertices",  60, -30, 30);
  TH1D *r_res = new TH1D("r_res", ";#Delta radial decay position [mm];Vertices",  60, -30, 30);
  TH1D *p_res = new TH1D("p_res", ";#Delta momentum [MeV];Vertices", 200, -100, 100);
  TH1D *theta_y_res = new TH1D("theta_y_res", ";#Delta #theta_{y} [mrad];Vertices", 120, -30, 30);

  TH1D *h_theta_y_true = new TH1D("theta_y_true", ";True #theta_{y} [mrad];Vertices", 400, -100, 100);
  TH1D *h_theta_y_reco = new TH1D("theta_y_reco", ";Reco #theta_{y} [mrad];Vertices", 400, -100, 100);

  // Position, momentum and vertical angle resolutions over a momentum range and as a function of momentum 
  TH2D *y_res_vs_p = new TH2D("y_res_vs_p", ";True momentum [MeV];#Delta vertical decay position [mm]",  30, 0, 3000, 60, -30, 30);
  TH2D *r_res_vs_p = new TH2D("r_res_vs_p", ";True momentum [MeV];#Delta radial decay position [mm]",  30, 0, 3000, 60, -30, 30);
  TH2D *p_res_vs_p = new TH2D("p_res_vs_p", ";True momentum [MeV];#Delta momentum [MeV]", 30, 0, 3000, 200, -100, 100);
  TH2D *theta_y_res_vs_p = new TH2D("theta_y_res_vs_p", ";True momentum [MeV];#Delta #theta_{y} [mrad]", 30, 0, 3000, 120, -30, 30);

  // Fill histograms
  InitBranches br(treeTracks);

  cout<<"-----> Filling histograms"<<endl;

  for(int64_t entry = 0; entry < treeTracks->GetEntries(); entry++) {

    treeTracks->GetEntry(entry);

    if(quality) {
      if(!br.passVertexQuality) continue;
      if(br.pValue < 0.05) continue;
      if(br.hitVolume) continue; 
    } 

    double y_true = br.trueVertexPosY;
    double r_true = sqrt(pow(br.trueVertexPosX,2)+pow(br.trueVertexPosZ,2)) - R_magic;
    double p_true = sqrt(pow(br.trueVertexMomX,2)+pow(br.trueVertexMomY,2)+pow(br.trueVertexMomZ,2));
    double theta_y_true = asin(br.trueVertexMomY/p_true) * 1e3; 

    double y_reco = br.recoVertexPosY;
    double r_reco = sqrt(pow(br.recoVertexPosX,2)+pow(br.recoVertexPosZ,2)) - R_magic;
    double p_reco = sqrt(pow(br.recoVertexMomX,2)+pow(br.recoVertexMomY,2)+pow(br.recoVertexMomZ,2));
    double theta_y_reco = asin(-br.recoVertexMomY/p_reco) * 1e3; 

    int stn = br.station;

    if(stn!=0) continue;

    if(br.trueTime < g2Period*7) continue;

    if(p_true > pLo && p_true < pHi) {

      y_res->Fill(y_true - y_reco);
      r_res->Fill(r_true - r_reco);
      p_res->Fill(p_true - p_reco);
      theta_y_res->Fill(theta_y_true - theta_y_reco);
      
      h_theta_y_true->Fill(theta_y_true);
      h_theta_y_reco->Fill(theta_y_reco);

    } 

    // 2D hists (no need for mom cut)
    y_res_vs_p->Fill(p_true, y_true - y_reco);
    r_res_vs_p->Fill(p_true, r_true - r_reco);
    p_res_vs_p->Fill(p_true, p_true - p_reco);
    theta_y_res_vs_p->Fill(p_true, theta_y_true - theta_y_reco);

  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

  string momSliceString = to_string(int(pLo))+"_"+to_string(int(pHi))+"_MeV";
  // Write to output

  // Set output directory

  output->mkdir(("ThetaY_"+momSliceString).c_str());
  output->mkdir(("Resolution_"+momSliceString).c_str());
  output->mkdir("ResolutionVsMomentum");

  cout<<"-----> Writing 1D hists"<<endl;

  output->cd(("ThetaY_"+momSliceString).c_str());

  h_theta_y_true->Write();
  h_theta_y_reco->Write();

  output->cd(("Resolution_"+momSliceString).c_str());

  y_res->Write();
  r_res->Write();
  p_res->Write();
  theta_y_res->Write();

  cout<<"-----> Writing 2D hists"<<endl;

  output->cd("ResolutionVsMomentum");

  y_res_vs_p->Write();
  r_res_vs_p->Write();
  p_res_vs_p->Write();
  theta_y_res_vs_p->Write();

  cout<<"-----> Written plots."<<endl;

  return;

}

int main(int argc, char *argv[]) {

  bool quality = true;
  bool timeCuts = true;

  string inFileName = argv[1]; 
  string outFileName = argv[2];

  string treeNameTracks = "trackerNTup/TrackerMCDecayTree";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName .c_str());

  // Get tree
  TTree *treeTracks = (TTree*)fin->Get(treeNameTracks.c_str());

  cout<<"\nOpened tree:\t"<<treeNameTracks<<" "<<treeTracks<<" from file "<<inFileName<<" "<<fin<<endl;

  // Book output
  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");

  Run(treeTracks, fout, quality);

  // Close
  fout->Close();
  fin->Close();

  cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;
   
  return 0;
}
