// Sam Grant 2021 

// We have both an offset in the vertical which is depedant on both momentum and time
// Here we produce the base histograms required for the correction: 
// TH2D of ThetaY_vs_Time in set momentum slices, binned in g-2 peroiods

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

double pLo = 750; 
double pHi = 2500;

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

void Run(TTree *tree, TFile *output, bool quality = true, string dataset = "Run-1a") {

  // ThetaY boost into lab frame
  double boostFactor = 5e3*(1/gmagic);

  // Slice momentum
  int step = 125;
  int nSlices = pmax/step;

  string stns[] = {"S12S18", "S12", "S18"}; 

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  TH2D *thetaY_vs_time_[n_stn];
  vector<TH2D *> thetaY_vs_time_slices_[n_stn];

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    thetaY_vs_time_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time").c_str(), ";Decay time [#mus];#theta_{y} [mrad] / 4.365 #mus", 92, 0, 92*g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor); 
    
    // Slice momentum
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = 0 + i_slice*step; 
      int hi = step + i_slice*step;

      TH2D *thetaY_vs_time_slice = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay time [#mus];#theta_{y} [mrad] / 4.365 #mus", 92, 0, 92*g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor); 
      thetaY_vs_time_slices_[i_stn].push_back(thetaY_vs_time_slice);

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

    double time = RandomisedTime(rand, br.decayTime * 1e-3); // randomise out the FR

    double px = br.decayVertexMomX; double py = -br.decayVertexMomY; double pz = br.decayVertexMomZ; 

    bool vertexQual = br.passDecayVertexQuality;

    // Time cuts
    if (quality && !vertexQual) continue;

    int stn = br.station; 

    // Positron world momentum
    TVector3 eMom(px, py, pz); 
    double p = eMom.Mag();
      
    /////////////////////////////////
    //                             //
    // Define the vertical angle   //
    //                             //
    /////////////////////////////////

    // If the you have the vertex cut switched on, py must be negative
    double theta_y = asin(py/p);
    theta_y = theta_y * 1e3;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // All stations 
    int stn_id = 0;

    if(p > pLo && p < pHi) thetaY_vs_time_[stn_id]->Fill(time, theta_y);

    // Slice momentum 
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = 0 + i_slice*step; 
      int hi = step + i_slice*step;

      if(p >= double(lo) && p < double(hi)) thetaY_vs_time_slices_[stn_id].at(i_slice)->Fill(time, theta_y);

    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Stn 12 or stn 18
    if(stn==12) stn_id = 1;
    else if(stn==18) stn_id = 2;
    else cerr<<"Station "<<stn<<" not recognised";

    if(p > pLo && p < pHi) thetaY_vs_time_[stn_id]->Fill(time, theta_y);

    // Slice momentum 
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = 0 + i_slice*step; 
      int hi = step + i_slice*step;

      if(p >= double(lo) && p < double(hi)) thetaY_vs_time_slices_[stn_id].at(i_slice)->Fill(time, theta_y);

    }

/*    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }*/

  }

  // Write to output
  // Set output directory
  output->mkdir("MainPlots"); output->mkdir("MomSlices");  

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    output->cd("MainPlots");

    thetaY_vs_time_[i_stn]->Write();

   for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) {

      output->cd("MomSlices"); 

      thetaY_vs_time_slices_[i_stn].at(i_slice)->Write();
      
    }

  }

  cout<<"Written plots."<<endl;

  return;

}

int main(int argc, char *argv[]) {

   bool quality = true;

   string inFileName = argv[1];
   string outFileName = argv[2];
   string dataset = argv[3];

/* string inFileName = "/gm2/data/g2be/Production/Trees/Run1/trackRecoTrees_15921.root";
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
