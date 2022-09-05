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

// Use the output from data to reweight 
double GetWeighting(TH2D *map, double theta_y, double p) {

  int i = map->GetXaxis()->FindBin(p);
  int j = map->GetYaxis()->FindBin(theta_y);

  double weighting = map->GetBinContent(i, j);

  if(isnan(weighting)) {
    return 0;
  }

  return weighting;

}

void Run(TTree *tree, TFile *output, bool quality, bool truth, string dataset = "Run-1a") {

  // Get weighting file
  string weightingFileName = "correctionHists/thetaYvsMomentum_"+dataset+"_BQ_noVertCorr.root";
  if(dataset=="Run-1d") weightingFileName = "correctionHists/thetaYvsMomentum_"+dataset+"_50usStartTime_BQ_noVertCorr.root";
  
  TFile *weightingFile = TFile::Open(weightingFileName.c_str());

  cout<<"---> Got weighting file "<<weightingFileName<<", "<<weightingFile<<endl;

  //string stns[] = {"S0S12S18", "S12S18", "S0", "S12", "S18"}; 
  string stns[] = {"S12S18", "S12", "S18"}; 

  cout<<"---> Booking histograms"<<endl;

  // For direct comp. with decays
  double boostFactor = 5e3*(1/gmagic);

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  TH2D *thetaY_vs_p_[n_stn];
  TH2D *weightHist_[n_stn];

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    thetaY_vs_p_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Momentum").c_str(), ";Momentum [MeV]; #theta_{y} [mrad] / 50 MeV ", int(pmax)/50, 0, pmax, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);//int(pmax), 0, pmax, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    weightHist_[i_stn] = (TH2D*)weightingFile->Get(("SanityPlots/"+stns[i_stn]+"_ThetaY_vs_Momentum").c_str());
    cout<<weightHist_[i_stn]<<endl;
    weightHist_[i_stn]->Scale(1./weightHist_[i_stn]->GetMaximum()); // Normalise 

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
    //int stn_id = -1;
    //if(stn==0) stn_id = 2;
    //else if(stn==12) stn_id = 3;
    //else if(stn==18) stn_id = 4;

    int stn_id = -1;
    //if(stn==0) stn_id = 2;
    if(stn==12) stn_id = 1;
    else if(stn==18) stn_id = 2;
    else continue;//cerr<<"Station "<<stn<<" not recognised";

    // if(weightHist_[stn_id]==0) weightHist_[stn_id]

    // For reweighting, just use S12 or S18
    //if(stn_id!=3 || stn_id!=4) continue;

    double weighting = GetWeighting(weightHist_[stn_id], theta_y, eMom.Mag());

    // cout<<"---> Filling histograms"<<endl;

    // Fill histograms 

    thetaY_vs_p_[stn_id]->Fill(eMom.Mag(), theta_y, weighting);


/*    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }*/

  }

  cout<<"---> Merging stations"<<endl;


  // Combine stations

  //thetaY_vs_p_[1]->Add(thetaY_vs_p_[3], thetaY_vs_p_[4]);
  thetaY_vs_p_[0]->Add(thetaY_vs_p_[1], thetaY_vs_p_[2]);


  // Write to output

  cout<<"---> Writing to output"<<endl;

  // Set output directory
  output->mkdir("SanityPlots"); output->cd("SanityPlots");

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

	thetaY_vs_p_[i_stn]->Write();
  }

  weightingFile->Close();

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