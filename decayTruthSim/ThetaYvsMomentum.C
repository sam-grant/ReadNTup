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
double R_magic = 7112; // mm

double pLo = 1000; 
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

void Run(TTree *tree, TFile *output, bool timeCuts, bool boost) { 

  // Somewhat arbitrary... 
  double boostFactor = 5e3*(1/gmagic);
  double momBoostFactor = 1.;

  if(boost) { 
    boostFactor = 1.0e3;
    momBoostFactor = (1/(2*gmagic));
  }

  // ------ Book histograms -------

  cout<<"---> Booking histograms"<<endl;

  TH2D *thetaY_vs_p = new TH2D("ThetaY_vs_Momentum", ";Momentum [MeV]; #theta_{y} [mrad] / 50 MeV ", int(pmax)/50, 0, pmax, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);

  // Get branches (using header file)
  InitBranches br(tree);

  double targetPerc = 0;
  int64_t nEntries = tree->GetEntries();

  double muAngleMax = 0;

  int64_t counter = 0;

  cout<<"---> Event loop"<<endl;

  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    // Get physics variables
    double t = br.posiInitTime * 1e-3; // us

    // Positron world momentum 
    TVector3 eMom(br.posiInitPX, br.posiInitPY, br.posiInitPZ); 
    TVector3 ePos(br.posiInitPosX, br.posiInitPosY, br.posiInitPosZ);
    TVector3 muPol(br.muDecayPolX, br.muDecayPolY, br.muDecayPolZ);

//    double y = ePos.Y(); 
    double r = sqrt(pow(ePos.X(),2)+pow(ePos.Z(), 2)) - R_magic;

    // RingAngle is angle from x axis, from 0 to 2pi
    double ringAngle = atan2(br.posiInitPosZ, br.posiInitPosX);    
    if (ringAngle < 0) ringAngle += TMath::TwoPi();

    // Positron angle around the ring momentum AAR
    // Z is tangential to magic mom at x and z of decay (Figure 2 of Debevec note)

    // Muon rest frame
    if(boost) {

      // Rotate into AAR
      eMom.RotateY(ringAngle);
      ePos.RotateY(ringAngle);
      muPol.RotateY(ringAngle);

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

    /////////////////////////////////
    //                             //
    // Define the vertical angle   //
    //                             //
    /////////////////////////////////

    // We need to do this in a way that doesn't involve the z-component of momentum
    // It should be an entirely tranvserse quantity 

    double x = ePos.X();
    double y = ePos.Y();
    double z = ePos.Z();

    double px = eMom.X();
    double py = eMom.Y();
    double pz = eMom.Z();
    double pT = sqrt( pow(px, 2) + pow(py, 2) );
    //double p = eMom.Mag();

    double theta_y = asin(py/eMom.Mag()) * 1e3;

    double alpha = muPol.Angle(eMom);

    ////////////////////////////////////////////////////////////////////////////////////////////////

    if(timeCuts && t < g2Period*7) continue; 

    // Fill histograms 
    thetaY_vs_p->Fill(eMom.Mag(), theta_y);
 
    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }
         
  }

  // Write to output
  // Set output directory
  output->mkdir("SanityPlots"); output->cd("SanityPlots");

  cout<<"---> Writing histograms"<<endl;

  thetaY_vs_p->Write();

  return;

}

int main(int argc, char *argv[]) {

  bool boost = false;
  bool timeCuts = true; 

  string inFileName = argv[1]; // "/pnfs/GM2/persistent/EDM/MC/dMu/TruthNTup/truthTrees.000.root";//
  string outFileName = argv[2]; //"tmp.root";

  string treeName = "phaseAnalyzer/g2phase";
  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName.c_str());
  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;

  // Book output
  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");

  Run(tree, fout, timeCuts, boost);
  
  fout->Close();
  fin->Close();

  cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;

  return 0;
}

