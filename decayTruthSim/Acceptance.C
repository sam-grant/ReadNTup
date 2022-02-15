// Read ROOT trees 
// Sam Grant

// Plot vertical momentum in slices of y-position

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

bool AcceptedDecay(double py, double y) { 

  if(y >= -50 && y < -40 && py >= -0.49 && py < 66.01) return true;  
  else if(y >= -40 && y < -30 && py >= -15.19 && py < 69.51) return true;  
  else if(y >= -30 && y < -20 && py >= -30.03 && py < 68.53) return true; 
  else if(y >= -20 && y < -10 && py >= -42.21 && py < 63.77) return true; 
  else if(y >= -10 && y < 0 && py >= -57.19 && py < 60.27) return true; 
  else if(y >= 0 && y < 10 && py >= -57.47 && py < 49.35) return true; 
  else if(y >= 10 && y < 20 && py >= -64.89 && py < 39.27) return true; 
  else if(y >= 20 && y < 30 && py >= -66.43 && py < 28.21) return true; 
  else if(y >= 30 && y < 40 && py >= -66.57 && py < 15.05) return true; 
  else if(y >= 40 && y < 50 && py >= -68.39 && py < -0.07) return true; 
  else return false;

}

void Run(TTree *tree, TFile *output, bool boost = false, bool acceptanceCorr = false) {

  // Somewhat arbitrary 
  double boostFactor = 5e3*(1/gmagic);
  double momBoostFactor = 1.;

  if(boost) { 
    boostFactor = 1.0e3;
    momBoostFactor = (1/(2*gmagic));
  }

  TH1D *momentumX = new TH1D("MomentumX", ";Decay intial momentum X [MeV];Decays", 700, -70, +70); //int(pmax), -pmax*momBoostFactor, pmax*momBoostFactor);  
  TH1D *momentumY = new TH1D("MomentumY", ";Decay intial momentum Y [MeV];Decays", 700, -70, +70); 

  TH1D *decayX = new TH1D("DecayX", ";Decay intial position X [MeV];Decays", 500, -50, +50);
  TH1D *decayY = new TH1D("DecayY", ";Decay intial position Y [MeV];Decays", 500, -50, +50);

  TH2D *decayX_vs_decayY = new TH2D("DecayX_vs_DecayY", ";Decay intial position X [mm];Decay intial Y [mm]", 500, -50, +50, 500, -50, +50);

  vector<TH1D*> momentumY_;
  vector<TH1D*> momentumX_;

  vector<TH1D*> decayY_;
  vector<TH1D*> decayX_;

  vector<TH2D*> decayX_vs_decayY_1_; // Vertical
  vector<TH2D*> decayX_vs_decayY_2_; // Radial

  // Slice vertical position
  int step = 10; // mm 
  int nSlices = 100/step; 

  // Slice position
  for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

    int lo = -50 + i_slice*step; 
    int hi = -50 + step + i_slice*step;

    TH1D *momentumY_slice = new TH1D(("MomentumY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay intial momentum Y [MeV];Decays", 1000, -70, +70); 
    momentumY_.push_back(momentumY_slice);

    TH1D *momentumX_slice = new TH1D(("MomentumX_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay intial momentum X [MeV];Decays", 1000, -70, +70); 
    momentumX_.push_back(momentumX_slice);

    TH1D *decayY_slice = new TH1D(("DecayY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay intial position Y [MeV];Decays", 1000, -100, +100);
    decayY_.push_back(decayY_slice);

    TH1D *decayX_slice = new TH1D(("DecayX_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay intial position X [MeV];Decays", 1000, -100, +100);
    decayX_.push_back(decayX_slice);

    TH2D *decayX_vs_decayY_1_slice = new TH2D(("DecayX_vs_DecayY_1_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay intial position X [mm];Decay intial Y [mm]", 500, -50, +50, 500, -50, +50);
    decayX_vs_decayY_1_.push_back(decayX_vs_decayY_1_slice);

    TH2D *decayX_vs_decayY_2_slice = new TH2D(("DecayX_vs_DecayY_2_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay intial position X [mm];Decay intial Y [mm]", 500, -50, +50, 500, -50, +50);
    decayX_vs_decayY_2_.push_back(decayX_vs_decayY_2_slice);

  }


  // Get branches (using header file)
  InitBranches br(tree);

  double targetPerc = 0;
  int64_t nEntries = tree->GetEntries();
   
  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    // Positron world momentum 
    TVector3 eMom(br.posiInitPX, br.posiInitPY, br.posiInitPZ); 
    TVector3 ePos(br.posiInitPosX, br.posiInitPosY, br.posiInitPosZ);
    TVector3 muPol(br.muDecayPolX, br.muDecayPolY, br.muDecayPolZ);

    // RingAngle is angle from x axis, from 0 to 2pi
    double ringAngle = atan2(ePos.z(), ePos.x()); // br.posiInitPosZ, br.posiInitPosX);  

    if (ringAngle < 0) ringAngle += TMath::TwoPi();

    // Positron angle around the ring momentum AAR
    // Z is tangential to magic mom at x and z of decay (Figure 2 of Debevec note)

    // Rotate into AAR
    eMom.RotateY(ringAngle);
    ePos.RotateY(ringAngle);
    muPol.RotateY(ringAngle);

    // Muon rest frame
    if(boost) {

      //rotation not perfect, so if original y component was 0 force that to be the case:
      if (fabs(muPol.y()) < 1E-10){
	     muPol = TVector3(muPol.x(), 0.0, muPol.z());
	     muPol = muPol.Unit();
      }

      double muP = br.muDecayP;
      double eP = br.posiInitP;
      double muE_lab = sqrt(mMu*mMu + muP*muP);
      double posiE_lab = sqrt(eMass*eMass + eP*eP);
      double gamma = sqrt(1.0 + pow( muP/mMu, 2 )); 

      // Construct muon momentum vector and rotate it
      TVector3 muMom(br.muDecayPX, br.muDecayPY, br.muDecayPZ);
      muMom.RotateY(ringAngle);

      // Construct 4-vector and boost it into the rest frame
      ROOT::Math::PxPyPzEVector muMomE_lab(muMom.x(), muMom.y(), muMom.z(), muE_lab);
      ROOT::Math::XYZVector boost = muMomE_lab.BoostToCM();

      // Get boost vector as a lorentz vector
      TVector3 boostToCM(boost.x(), boost.y(), boost.z());

      // Make positron momentum 4 vector in lab, and get boost vector as a TVector3 vector
      TLorentzVector posiMomE_lab(eMom.x(), eMom.y(), eMom.z(), posiE_lab);

      // Apply boost to positron lab vector to get it in MRF 1, using TLorentzVectors. Convert to normal TVector3 at end 
      TLorentzVector posiMomE_MRF = posiMomE_lab;
      posiMomE_MRF.Boost(boostToCM);

      TVector3 posiMom_MRF(posiMomE_MRF.Px(), posiMomE_MRF.Py(), posiMomE_MRF.Pz());

      eMom = posiMom_MRF; 


    } 

    double x = ePos.x() - Rmagic;
    double y = ePos.y();

    double px =  eMom.x();
    double py = eMom.y();

    // Crude acceptance correction
    if(acceptanceCorr && !AcceptedDecay(py, y)) continue;

    // Fill histograms
    decayY->Fill(y);
    decayX->Fill(x);
    momentumY->Fill(py);
    momentumX->Fill(px);
    decayX_vs_decayY->Fill(x, y);

    // Slice y-position
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = -50 + i_slice*step; 
      int hi = -50 + step + i_slice*step;

      if(y >= double(lo) && y < double(hi)) { 

        decayY_.at(i_slice)->Fill(y);
        momentumY_.at(i_slice)->Fill(py);
        decayX_vs_decayY_1_.at(i_slice)->Fill(x, y);

      }

      if(x >= double(lo) && x < double(hi)) { 

        decayX_.at(i_slice)->Fill(x);
        momentumX_.at(i_slice)->Fill(px);
        decayX_vs_decayY_2_.at(i_slice)->Fill(x, y);

      }

    }
 

    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }

  }

/*  cout<<"Muon max angle:\t"<<muAngleMax<<" radians"<<endl;*/
  // Write to output
  // Set output directory
  output->mkdir("SanityPlots"); 
  output->mkdir("VerticalMomentumSlices");
  output->mkdir("RadialMomentumSlices");

  output->cd("SanityPlots");

  momentumX->Write();
  momentumY->Write();
  decayX->Write();
  decayY->Write();
  decayX_vs_decayY->Write();

  for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) {

    output->cd("VerticalMomentumSlices"); 
    decayY_.at(i_slice)->Write();
    momentumY_.at(i_slice)->Write();
    decayX_vs_decayY_1_.at(i_slice)->Write();

    output->cd("RadialMomentumSlices");
    decayX_.at(i_slice)->Write();
    momentumX_.at(i_slice)->Write();
    decayX_vs_decayY_2_.at(i_slice)->Write(); 

  }

  return;

}

int main(int argc, char *argv[]) {

  bool boost = false;
  bool acceptanceCorr = false;

  string inFileName = argv[1]; 
  string outFileName = argv[2];

  string treeName = "phaseAnalyzer/g2phase";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName .c_str());
  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;

  // Book output
  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");
   
  // Fill histograms
  Run(tree, fout, boost, acceptanceCorr);

  // Close
  fout->Close();
  fin->Close();

  cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;
   
  return 0;

}
