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
#include "TRandom3.h"

using namespace std;

double omegaAMagic = 0.00143934; // from gm2geom consts / kHz 
double g2Period = (2*TMath::Pi()/omegaAMagic) * 1e-3; // 4.3653239 us
double mMu = 105.6583715; // MeV
double aMu = 11659208.9e-10; 
double gmagic = std::sqrt( 1.+1./aMu );
double pmax = 1.01 * mMu * gmagic; // 3127.1144
double T_c = 149.2 * 1e-3; // cyclotron period [us]

// Momentum
double xmin = 750; 
double xmax = 2500; 

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

double RandomisedTime(TRandom3 *rand, double time) { 
  return rand->Uniform(time-T_c/2, time+T_c/2);
}

double VertOffsetCorr(TH1D *c_vs_mom, double theta_y, double p) { 
  double c = c_vs_mom->GetBinContent( c_vs_mom->FindBin(p) );
  return theta_y - c;
}

void Run(TTree *tree, TFile *output, bool quality) { // , string dataset = "Run-1a") {

  double boostFactor = 5e3*(1/gmagic);

  string stns[] = {"S12S18", "S12", "S18"}; 

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  // Vertical offset correction / momentum
  TString finName = "correctionPlots/verticalOffset_"+dataset+".root";
  TFile *f_c_vs_mom = TFile::Open(finName); 
  vector<TH1D*> px_c_vs_mom_; 

  TH2D *theta_y_vs_t_all_[n_stn];
  TH2D *theta_y_vs_t_750MeV_2500MeV_[n_stn];

  TH2D *theta_y_vs_p_125MeV_[n_stn];
  TH2D *theta_y_vs_p_250MeV_[n_stn];

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    // px_c_vs_mom_.push_back((TH1D*)f_c_vs_mom->Get(("VerticalOffsetPlots/"+to_string(step)+"MeV/"+stns[i_stn]+"_px_ThetaY_vs_Momentum").c_str()));

    theta_y_vs_t_all_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_All").c_str(), ";Decay time [#mus];#theta_{y} [mrad] / 4.365 #mus", 92, 0, 92*g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor); 
    theta_y_vs_t_750MeV_2500MeV_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_750MeV_2500MeV").c_str(), ";Decay time [#mus];#theta_{y} [mrad] / 4.365 #mus", 92, 0, 92*g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor); 

    theta_y_vs_p_125MeV_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Momentum_125MeV").c_str(), ";Decay vertex momentum [MeV];#theta_{y} [mrad] / 125 Mev", 24, 0, 3000, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor); 
    theta_y_vs_p_250MeV_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Momentum_250Mev").c_str(), ";Decay vertex momentum [MeV];#theta_{y} [mrad] / 250 Mev", 12, 0, 3000, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor); 

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

    double time = RandomisedTime(rand, br.decayTime * 1e-3); // randomise out the FR

    double px = br.decayVertexMomX; double py = -br.decayVertexMomY; double pz = br.decayVertexMomZ; 

    bool hitVol = br.hitVolume;
    double pVal = br.trackPValue;
    bool vertexQual = br.passDecayVertexQuality;

    // Time cuts
    if (quality) {

      // if (time < g2Period*7 || time > g2Period*70) continue; 
      // if (hitVol) continue;
      // if (pVal < 0.05) continue; 

      // Should include hit vol and pval cuts
      if (!vertexQual) continue;

    }

    // Get more analysis variables
    // double time = br.decayTime * 1e-3; // us

    int stn = br.station; 

    // Positron world momentum
    TVector3 eMom(px, py, pz); 
      
    /////////////////////////////////
    //                             //
    // Define the vertical angle   //
    //                             //
    /////////////////////////////////

    double p = eMom.Mag();

    // If the you have the vertex cut switched on, py must be negative
    double theta_y = asin(py/p);

    //double alpha = muPol.Angle(eMom);

    theta_y = theta_y * 1e3;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // All stations 
    int stn_id = 0;

    // Vertical correction
/*    if(true) {
      theta_y = VertOffsetCorr(px_c_vs_mom_.at(stn_id), theta_y, p); // step, nSlices, 
    }
*/
    theta_y_vs_t_all_[stn_id]->Fill(time, theta_y);

    if(p > xmin && p < xmax) theta_y_vs_t_750MeV_2500MeV_[stn_id]->Fill(time, theta_y);

    if (quality && (time > g2Period*7 || time < g2Period*70)) {
      theta_y_vs_p_125MeV_[stn_id]->Fill(p, theta_y);
      theta_y_vs_p_250MeV_[stn_id]->Fill(p, theta_y);
    } 

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Stn 12 or stn 18
    if(stn==12) stn_id = 1;
    else if(stn==18) stn_id = 2;
    else cerr<<"Station "<<stn<<" not recognised";

    theta_y_vs_t_all_[stn_id]->Fill(time, theta_y);
    if(p > xmin && p < xmax) theta_y_vs_t_750MeV_2500MeV_[stn_id]->Fill(time, theta_y);
    
    if (quality && (time > g2Period*7 || time < g2Period*70)) {
      theta_y_vs_p_125MeV_[stn_id]->Fill(p, theta_y);
      theta_y_vs_p_250MeV_[stn_id]->Fill(p, theta_y);
    } 


    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }

  }

  // Write to output
  // Set output directory
  output->mkdir("VerticalOffsetPlots"); output->cd("VerticalOffsetPlots");  

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    theta_y_vs_t_all_[i_stn]->Write();
    theta_y_vs_t_750MeV_2500MeV_[i_stn]->Write();
    theta_y_vs_p_125MeV_[i_stn]->Write();
    theta_y_vs_p_250MeV_[i_stn]->Write();  

  }

  cout<<"Written plots."<<endl;

  return;

}

int main(int argc, char *argv[]) {

   bool quality = true;

   string inFileName = argv[1];
   string outFileName = argv[2];

/*   string inFileName = "/gm2/data/g2be/Production/Trees/Run1/trackRecoTrees_15921.root";
   string outFileName = "trackRecoPlots_15921.root";*/

   string treeName = "trackAndTrackCalo/tree";

   // Open tree and load branches
   TFile *fin = TFile::Open(inFileName .c_str());
   // Get tree
   TTree *tree = (TTree*)fin->Get(treeName.c_str());

   cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;

   // Book output
   TFile *fout = new TFile(outFileName.c_str(),"RECREATE");
   
   // Fill histograms
   Run(tree, fout, quality);

   // Close
   fout->Close();
   fin->Close();

   cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;
   
   return 0;

}
