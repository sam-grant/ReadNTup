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

double ModTime(double time) {

  double g2fracTime = time / g2Period;
  int g2fracTimeInt = int(g2fracTime);
  double g2ModTime = (g2fracTime - g2fracTimeInt) * g2Period;

  return g2ModTime;

}

double RandomisedTime(TRandom3 *rand, double time) { 
  return rand->Uniform(time-T_c/2, time+T_c/2);
}

void Run(TTree *tree, TFile *output, bool quality) {

  // 0: WORLD; 1: AAR; 2: MRF
  int frame = 1; 

  bool boost;

  bool WORLD = false;
  bool AAR = false;
  bool MRF = false;

  // Need to clean this up
  double boostFactor;
  double momBoostFactor;

  if(frame==0) { 
    WORLD = true;
    // Factor of two is just to make sure we're not chopping the tops of the vertcial angle at low momentum 
    boostFactor = 5e3*(1/gmagic);// 2/gmagic; // 1/15;//20.;
    boost = false;
    momBoostFactor = 1;
  } 
  else if(frame==1) {
    AAR = true;
    boostFactor = 5e3*(1/gmagic); //posiMom_MRF1/15;//20.;
    boost = false; 
    momBoostFactor = 1;
  }
  else if(frame==2) {
    MRF = true;
    boostFactor = 1.0e3;///20.;//2/gmagic;
    boost = true;
    momBoostFactor = (1/(2*gmagic));
  }

  string stns[] = {"S12S18", "S12", "S18"}; 

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  TH1D *wiggle_5ns = new TH1D("S12S18_Wiggle_5ns", ";Decay time [#mus];Vertices / 5 ns", 8731, 7*g2Period, 17*g2Period);//_[n_stn];
  TH1D *wiggle_5ns_rand = new TH1D("S12S18_Wiggle_5ns_rand", ";Decay time [#mus];Vertices / 5 ns", 8731, 7*g2Period, 17*g2Period);// _[n_stn];

/*  TH2D *thetaY_vs_time_mod_50ns_[n_stn];
  TH2D *thetaY_vs_time_mod_50ns_rand_[n_stn]; 
  TH2D *thetaY_vs_time_mod_5ns_[n_stn];
  TH2D *thetaY_vs_time_mod_5ns_rand_[n_stn];
  TH2D *thetaY_vs_time_mod_149ns_[n_stn];
  TH2D *thetaY_vs_time_mod_149ns_rand_[n_stn];

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    thetaY_vs_time_mod_50ns_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_50ns").c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 50 ns", 87, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor); 
    thetaY_vs_time_mod_50ns_rand_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_50ns_Rand").c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 50 ns", 87, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor); 

    thetaY_vs_time_mod_5ns_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_5ns").c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 5 ns", 873, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor); 
    thetaY_vs_time_mod_5ns_rand_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_5ns_Rand").c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 5 ns", 873, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor); 

    thetaY_vs_time_mod_149ns_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_149ns").c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 149.2 ns", 29, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor); 
    thetaY_vs_time_mod_149ns_rand_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_149ns_Rand").c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 149.2 ns", 29, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor); 

  }*/

  // Get branches (using header file)
  InitBranches br(tree);

  double targetPerc = 0;
  int64_t nEntries = tree->GetEntries();
   
  // For FR randomisation
  TRandom3 *rand = new TRandom3(12345);

  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    // quality variables
    double time = br.decayTime * 1e-3; // ns -> us
    double randTime = RandomisedTime(rand, time); // randomise out the FR

    bool hitVol = br.hitVolume;
    double pVal = br.trackPValue;
    bool vertexQual = br.passDecayVertexQuality;

    // Time cuts
    if (quality) {

      if (time < g2Period*7 || time > g2Period*70) continue; 
      if (hitVol) continue;
      if (pVal < 0.05) continue; 
      if (!vertexQual) continue;

    }

    // Get more analysis variables
    // double time = br.decayTime * 1e-3; // us

    int stn = br.station; 

    // Positron world momentum and position
    double x = br.decayVertexPosX; double y = br.decayVertexPosY; double z = br.decayVertexPosZ; 
    double px = br.decayVertexMomX; double py = -br.decayVertexMomY; double pz = br.decayVertexMomZ; 

    TVector3 ePos(x, y, z);
    TVector3 eMom(px, py, pz); 

    double p = eMom.Mag();

    double g2ModTime = ModTime(time);

    // Not so sure about this. I think that it shouldn't make any difference.
    //double g2ModRandTime = RandomisedTime(g2ModTime);
    double g2ModRandTime = ModTime(randTime);

    // RingAngle is angle from x axis, from 0 to 2pi
    double ringAngle = atan2(z, x);    
    if (ringAngle < 0) ringAngle += TMath::TwoPi();

    // Positron angle around the ring momentum AAR
    // Z is tangential to magic mom at x and z of decay (Figure 2 of Debevec note)
    if(AAR) { 

      eMom.RotateY(ringAngle);
      ePos.RotateY(ringAngle);

    } 
      
    /////////////////////////////////
    //                             //
    // Define the vertical angle   //
    //                             //
    /////////////////////////////////

    // If the you have the vertex cut switched on, py must be negative
    double theta_y = asin(py/p);

    //double alpha = muPol.Angle(eMom);

    theta_y = theta_y * 1e3;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // All stations 
    int stn_id = 0;

    // g-2 cuts. See Fienberg thesis figure 2.10
    if(p > 1900*momBoostFactor  && p < pmax*momBoostFactor) {
      wiggle_5ns->Fill(time);
      wiggle_5ns_rand->Fill(randTime);
    } 

    // EDM cuts
    if( quality && p > 750*momBoostFactor && p < 2500*momBoostFactor) { 
/*
    thetaY_vs_time_mod_50ns_[stn_id]->Fill(g2ModTime, theta_y);
    thetaY_vs_time_mod_50ns_rand_[stn_id]->Fill(g2ModRandTime, theta_y);

    thetaY_vs_time_mod_5ns_[stn_id]->Fill(g2ModTime, theta_y);
    thetaY_vs_time_mod_5ns_rand_[stn_id]->Fill(g2ModRandTime, theta_y);

    thetaY_vs_time_mod_149ns_[stn_id]->Fill(g2ModTime, theta_y);
    thetaY_vs_time_mod_149ns_rand_[stn_id]->Fill(g2ModRandTime, theta_y);*/

    } else if(!quality) { 

    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Stn 12 or stn 18
    if(stn==12) stn_id = 1;
    else if(stn==18) stn_id = 2;
    else cerr<<"Station "<<stn<<" not recognised";

    // EDM cuts
    if( quality && p > 750*momBoostFactor && p < 2500*momBoostFactor) { 
/*
    thetaY_vs_time_mod_50ns_[stn_id]->Fill(g2ModTime, theta_y);
    thetaY_vs_time_mod_50ns_rand_[stn_id]->Fill(g2ModRandTime, theta_y);

    thetaY_vs_time_mod_5ns_[stn_id]->Fill(g2ModTime, theta_y);
    thetaY_vs_time_mod_5ns_rand_[stn_id]->Fill(g2ModRandTime, theta_y);

    thetaY_vs_time_mod_149ns_[stn_id]->Fill(g2ModTime, theta_y);
    thetaY_vs_time_mod_149ns_rand_[stn_id]->Fill(g2ModRandTime, theta_y);*/

    } else if(!quality) { 

    }


    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }

  }

  // Write to output
  // Set output directory
  output->mkdir("FastRotationPlots");
  output->cd("FastRotationPlots");

  wiggle_5ns->Write();
  wiggle_5ns_rand->Write();

/*  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    thetaY_vs_time_mod_50ns_[i_stn]->Write();
    thetaY_vs_time_mod_50ns_rand_[i_stn]->Write();

    thetaY_vs_time_mod_5ns_[i_stn]->Write();
    thetaY_vs_time_mod_5ns_rand_[i_stn]->Write();

    thetaY_vs_time_mod_149ns_[i_stn]->Write();
    thetaY_vs_time_mod_149ns_rand_[i_stn]->Write();

  }*/

  cout<<"Written plots."<<endl;

  return;

}

int main(int argc, char *argv[]) {

   bool quality = true;

   string inFileName = argv[1]; //"/gm2/data/g2be/Production/Trees/Run1/trackRecoTrees_15921.root";//"/gm2/data/g2be/Production/Trees/Run1/trackRecoTrees_15921.root";// "/pnfs/GM2/persistent/EDM/MC/dMu/TruthNTup/truthTrees.000.root";//
   string outFileName = argv[2]; //"tmp.root"; // //"tmp.root"; // argv[2]; //"tmp.root";

   // string inFileName = "/gm2/data/g2be/Production/Trees/Run1/trackRecoTrees_15921.root";//"/gm2/data/g2be/Production/Trees/Run1/trackRecoTrees_15921.root";// "/pnfs/GM2/persistent/EDM/MC/dMu/TruthNTup/truthTrees.000.root";//
   // string outFileName = "fastRotationPlots.root"; // //"tmp.root"; // argv[2]; //"tmp.root";

   string treeName = "trackAndTrackCalo/tree"; //trackerNTup/tracker";

   // Open tree and load branches
   TFile *fin = TFile::Open(inFileName .c_str());
   // Get tree
   TTree *tree = (TTree*)fin->Get(treeName.c_str());

   cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;

   // Book output

   TFile *fout = new TFile(outFileName.c_str(),"RECREATE");
   
   // Fill histograms
   Run(tree, fout, quality);
   //RunEqualStatsBins(tree, fout, quality);

   // Close
   fout->Close();
   fin->Close();

   cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;
   
   return 0;

}
