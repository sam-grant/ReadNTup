// Plot theta_y, extrap distance, and Y in momentum bins

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


void Run(TTree *tree, TFile *output, bool quality) { 


  string stns[] = {"S0S12S18", "S12S18", "S0", "S12", "S18"}; 

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  TH1D *p_[n_stn];
  TH1D *Y_[n_stn];
  TH1D *thetaY_[n_stn];
  TH1D *extrapDist_[n_stn];

  // momentum scans
  vector<TH1D*> p_slices_[n_stn];
  vector<TH1D*> Y_slices_[n_stn];
  vector<TH1D*> thetaY_slices_[n_stn];
  vector<TH1D*> extrapDist_slices_[n_stn];

  // Slice momentum
  int step = 250; 
  int nSlices = pmax/step;

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    p_[i_stn] = new TH1D((stns[i_stn]+"_Momentum").c_str(), ";Decay vertex momentum [MeV];Vertices", int(pmax), 0, pmax); 
    Y_[i_stn] = new TH1D((stns[i_stn]+"_Y").c_str(), ";Decay vertex position Y [mm];Vertices", 1000, -60, 60); 
 	  thetaY_[i_stn] = new TH1D((stns[i_stn]+"_ThetaY").c_str(), ";#theta_{y} [mrad];Vertices", 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
 	  extrapDist_[i_stn] = new TH1D((stns[i_stn]+"_ExtrapolatedDistance").c_str(), ";Decay extrapolated distance [mmm];Vertices", 5000, 0, 5000); // -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
     
    // Slice momentum
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = 0 + i_slice*step; 
      int hi = step + i_slice*step;

      p_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_Momentum_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay vertex momentum [MeV];Vertices", int(pmax), 0, pmax));
      Y_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_Y_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay vertex position Y [mm];Vertices", 1000, -60, 60)); 
      thetaY_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_ThetaY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";#theta_{y} [mrad];Tracks",  1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic));
      extrapDist_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_ExtrapolatedDistance_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay extrapolated distance [mmm];Vertices", 5000, 0, 5000)); 

    }

  }

  // Get branches (using header file)
  InitBranches br(tree);

  int64_t nEntries = tree->GetEntries();

  double targetPerc = 0;
  double muAngleMax = 0;

  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    bool vertexQual = br.passVertexQuality;
    int stn = br.station; 
    double time = br.decayTime * 1e-3; // us

    TVector3 ePos(br.decayVertexPosX, br.decayVertexPosY, br.decayVertexPosZ);
    TVector3 eMom(br.decayVertexMomX, -br.decayVertexMomY, br.decayVertexMomZ); 

    double p = eMom.Mag();
    double extrapDist = br.extrapolatedDistance;

    // Time cuts
    if (quality) {
      if (time < g2Period*7) continue; 
      if (!vertexQual) continue;
    }

    /////////////////////////////////
    //                             //
    // Define the vertical angle   //
    //                             //
    /////////////////////////////////

    double theta_y = asin(eMom.y()/p);

    theta_y = theta_y * 1e3;

    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Fill stations invidually according the station array
    int stn_id = -1;
    if(stn==0) stn_id = 2;
    else if(stn==12) stn_id = 3;
    else if(stn==18) stn_id = 4;
    else cerr<<"Station "<<stn<<" not recognised";

    // Fill histograms
    p_[stn_id]->Fill(p);
    Y_[stn_id]->Fill(ePos.y());
    thetaY_[stn_id]->Fill(theta_y);
    extrapDist_[stn_id]->Fill(extrapDist);

    // Slice momentum 
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

    	int lo = 0 + i_slice*step; 
    	int hi = step + i_slice*step;

    	if(p >= double(lo) && p < double(hi)) { 

    		p_slices_[stn_id].at(i_slice)->Fill(p);
    		Y_slices_[stn_id].at(i_slice)->Fill(ePos.y());
    		thetaY_slices_[stn_id].at(i_slice)->Fill(theta_y);
    		extrapDist_slices_[stn_id].at(i_slice)->Fill(extrapDist);

    	}

    }

    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }

  }

  // Combine stations S12&S18
  p_[1]->Add(p_[3],p_[4]);
  Y_[1]->Add(Y_[3], Y_[4]);
  thetaY_[1]->Add(thetaY_[3], thetaY_[4]);
  extrapDist_[1]->Add(extrapDist_[3], extrapDist_[4]);

   for ( int i_slice(0); i_slice < nSlices; i_slice++ ) {

      p_slices_[1].at(i_slice)->Add(p_slices_[3].at(i_slice), p_slices_[4].at(i_slice));
      Y_slices_[1].at(i_slice)->Add(Y_slices_[3].at(i_slice), Y_slices_[4].at(i_slice));
      thetaY_slices_[1].at(i_slice)->Add(thetaY_slices_[3].at(i_slice), thetaY_slices_[4].at(i_slice));
      extrapDist_slices_[1].at(i_slice)->Add(extrapDist_slices_[3].at(i_slice), extrapDist_slices_[4].at(i_slice));
      
  }

  // Combine stations S0&S12&S18
  p_[0]->Add(p_[1],p_[2]);
  Y_[0]->Add(Y_[1], Y_[2]);
  thetaY_[0]->Add(thetaY_[1], thetaY_[2]);
  extrapDist_[0]->Add(extrapDist_[1], extrapDist_[2]);

   for ( int i_slice(0); i_slice < nSlices; i_slice++ ) {

      p_slices_[0].at(i_slice)->Add(p_slices_[1].at(i_slice), p_slices_[2].at(i_slice));
      Y_slices_[0].at(i_slice)->Add(Y_slices_[1].at(i_slice), Y_slices_[2].at(i_slice));
      thetaY_slices_[0].at(i_slice)->Add(thetaY_slices_[1].at(i_slice), thetaY_slices_[2].at(i_slice));
      extrapDist_slices_[0].at(i_slice)->Add(extrapDist_slices_[1].at(i_slice), extrapDist_slices_[2].at(i_slice));
      
  }

  // Write to output
  // Set output directory
  output->mkdir("SimultaneousAnalysis"); output->mkdir("MomentumBinnedAnalysis");  

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    output->cd("SimultaneousAnalysis");

   	p_[i_stn]->Write();
    Y_[i_stn]->Write();
    thetaY_[i_stn]->Write();
    extrapDist_[i_stn]->Write();

   for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) {

      output->cd("MomentumBinnedAnalysis"); 

      p_slices_[i_stn].at(i_slice)->Write();
      Y_slices_[i_stn].at(i_slice)->Write();
      thetaY_slices_[i_stn].at(i_slice)->Write();
      extrapDist_slices_[i_stn].at(i_slice)->Write();
      
    }

  }

  cout<<"Written plots."<<endl;

  return;

}

int main(int argc, char *argv[]) {

  bool quality = true;

  string inFileName = argv[1]; 
  string outFileName = argv[2];

  string treeName = "trackerNTup/tracker";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName.c_str());

  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;

  // Book output

  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");
   
  // Fill histograms
  Run(tree, fout, quality); // , verticalOffsetCorrection, acceptanceCorrection);

  // Close
  fout->Close();
  fin->Close();

  cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;
   
  return 0;
}
