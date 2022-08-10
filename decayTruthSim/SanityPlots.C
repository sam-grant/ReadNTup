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

double ModTime(double time, int nPeriods = 1) {

  double g2fracTime = time / (nPeriods * g2Period);
  int g2fracTimeInt = int(g2fracTime);
  double g2ModTime = (g2fracTime - g2fracTimeInt) * nPeriods * g2Period;

  return g2ModTime;

}

double RandomisedTime(TRandom3 *rand, double time) { 
  return rand->Uniform(time-T_c/2, time+T_c/2);
}

double AcceptanceWeightedAngle(TH2D *map, double theta_y, double y) { 

  int i = map->GetXaxis()->FindBin(y);
  int j = map->GetYaxis()->FindBin(theta_y);

  double weighting = map->GetBinContent(i, j);

  if(isnan(weighting)) {
    return 0;
  }

  return weighting;


}

// I don't think this helps. 
double AcceptanceWeightedAngleWithInterpolation(TH2D *map1, TH2D *map2, double theta_y, double y) { 

  // Get coordinates
  int i = map1->GetXaxis()->FindBin(y);
  int j = map1->GetYaxis()->FindBin(theta_y);

  // Get main weighting
  double weighting1 = map1->GetBinContent(i, j);
  //double err_weighting1 = map1->GetBinError(i, j);

  if(isnan(weighting1)) return weighting1;

  // For min/max momentum 
  if(map2==0) return weighting1;

  // Get secondary weighting 

  // Get coordinates
  int i2 = map2->GetXaxis()->FindBin(y);
  int j2 = map2->GetYaxis()->FindBin(theta_y);

  double weighting2 = map2->GetBinContent(i2, j2);

  // Interpolate
  double weighting = (weighting1 + weighting2)/2;

  if(isnan(weighting)) {
    //cout<<"WARNING: acceptance wieghting is nan"<<endl;
    return 0;//theta_y;
  }

  //cout<<theta_y * weighting<<endl;
  //cout<<weighting<<endl;

  return weighting;


}

void Run(TTree *tree, TFile *output, bool momCuts, bool timeCuts, bool boost, bool verticalOffsetCorrection, bool acceptanceCorr, string stn) {

  // Get vertical offset correction histograms
  TString verticalOffsetFileName = "correctionHists/verticalOffsetHists_allDecays_WORLD_250MeV_AQ.root";
  TFile *verticalOffsetFile = TFile::Open(verticalOffsetFileName);
  TH1D *verticalOffsetHist = (TH1D*)verticalOffsetFile->Get("VerticalOffsetHists/ThetaY_vs_p");

  TString acceptanceFileName = "correctionHists/acceptanceWeightingPlots.truth.root";
  TFile *acceptanceFile = TFile::Open(acceptanceFileName);

  //TH2D *acceptanceHistAllMom = (TH2D*)acceptanceFile->Get(("AcceptanceWeighting/AllMom/750_2500_MeV/"+stn+"_AcceptanceMap").c_str());
  //cout<<acceptanceHistAllMom<<endl;

  vector<TH2D*> acceptanceHists_; 
  // Set the number of periods for the longer modulo plots
  int moduloMultiple = 4; 

  // Somewhat arbitrary 
  double boostFactor = 5e3*(1/gmagic);
  double momBoostFactor = 1.;

  if(boost) { 
    boostFactor = 1.0e3;
    momBoostFactor = (1/(2*gmagic));
  }

  // ------ Book histograms -------

  TH1D *momentum = new TH1D("Momentum", ";Track momentum [MeV];Tracks", int(pmax), 0, pmax*momBoostFactor); 
  TH1D *momY = new TH1D("MomentumY", ";Track momentum Y [MeV];Tracks", 1000, -60, 60); 
  TH1D *momX = new TH1D("MomentumX", ";Track momentum X [MeV];Tracks", int(pmax), -pmax*momBoostFactor, pmax*momBoostFactor); 
  TH1D *momZ = new TH1D("MomentumZ", ";Track momentum Z [MeV];Tracks", int(pmax), -pmax*momBoostFactor, pmax*momBoostFactor); 
  TH2D *decayZ_vs_decayX = new TH2D("DecayZ_vs_DecayX", ";Decay vertex position X [mm];Decay vertex position Z [mm]", 800, -8000, 8000, 800, -8000, 8000);
  TH1D *wiggle = new TH1D("Wiggle", ";Decay time [#mus];Tracks", 2700, 0, 2700*T_c);
  TH1D *wiggle_mod = new TH1D("Wiggle_Modulo", ";t_{g#minus2}^{mod} [#mus];Tracks / 149.2 ns", 29, 0, g2Period); 
  TH1D *wiggle_mod_long = new TH1D("Wiggle_Modulo_Long", ";Time modulo [#mus];Tracks / 149.2 ns", 41*moduloMultiple, 0, g2Period*moduloMultiple); 
  TH1D *thetaY = new TH1D("ThetaY", ";#theta_{y} [mrad];Tracks", 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_time = new TH2D("ThetaY_vs_Time", ";Decay time [#mus]; #theta_{y} [mrad] / 149.2 ns ", 2700, 0, 2700*T_c, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_time_20ns = new TH2D("ThetaY_vs_Time_20ns", ";Decay time [#mus]; #theta_{y} [mrad] / 20 ns ", 20000, 0, 20000*20e-3, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_time_50ns = new TH2D("ThetaY_vs_Time_50ns", ";Decay time [#mus]; #theta_{y} [mrad] / 50 ns ", 8000, 0, 8000*50e-3, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_time_mod = new TH2D("ThetaY_vs_Time_Modulo", ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 149.2 ns", 29, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_time_mod_50ns = new TH2D("ThetaY_vs_Time_Modulo_50ns", ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 50 ns", 87, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_time_mod_20ns = new TH2D("ThetaY_vs_Time_Modulo_20ns", ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 20 ns", 174, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_time_mod_long = new TH2D("ThetaY_vs_Time_Modulo_Long", ";Time modulo [#mus]; #theta_{y} [mrad] / 149.2 ns", 29*moduloMultiple, 0, g2Period*moduloMultiple, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_time_mod_long_20ns = new TH2D("ThetaY_vs_Time_Modulo_Long_20ns", ";Time modulo [#mus]; #theta_{y} [mrad] / 20 ns", 174*moduloMultiple, 0, g2Period*moduloMultiple, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_time_mod_long_50ns = new TH2D("ThetaY_vs_Time_Modulo_Long_50ns", ";Time modulo [#mus]; #theta_{y} [mrad] / 50 ns", 87*moduloMultiple, 0, g2Period*moduloMultiple, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_momentum = new TH2D("ThetaY_vs_Momentum", ";Decay vertex momentum [MeV]; #theta_{y} [mrad] / 10 MeV ", 300, 0, 3000, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_Y = new TH2D("ThetaY_vs_Y", ";Decay y-position [mm];#theta_{y} [mrad]", 24, -60, 60, 40, -100, 100);

  // Momentum scans
  vector<TH1D*> thetaY_mom_slices_;
  vector<TH1D*> Y_mom_slices_;
  vector<TH1D*> pY_mom_slices_;
  vector<TH1D*> p_mom_slices_;
  vector<TH2D*> thetaY_vs_Y_mom_slices_;
  vector<TH2D*> thetaY_vs_time_mod_slices_;
  vector<TH2D*> thetaY_vs_time_mod_20ns_slices_;
  vector<TH2D*> thetaY_vs_time_mod_50ns_slices_;
  vector<TH2D*> thetaY_vs_time_mod_long_slices_;
  vector<TH2D*> thetaY_vs_time_mod_long_20ns_slices_;
  vector<TH2D*> thetaY_vs_time_mod_long_50ns_slices_;
  

  // Slice momentum
  int step = 250 * momBoostFactor;
  int nSlices = (pmax/step) * momBoostFactor;

  // Slice momentum
  for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

    int lo = 0 + i_slice*step; 
    int hi = step + i_slice*step;

    TH1D *h_p_mom_slice = new TH1D(("Momentum_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Track momentum [MeV];Tracks",  int(pmax), 0, pmax*momBoostFactor);
    p_mom_slices_.push_back(h_p_mom_slice);

    TH1D *h_pY_mom_slice = new TH1D(("MomentumY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Track momentum Y MeV];Tracks",  1000, -60, 60);
    pY_mom_slices_.push_back(h_pY_mom_slice);

    TH1D *h_Y_mom_slice = new TH1D(("Y_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Vertical decay position [mm];Tracks",  180, -60, 60);
    Y_mom_slices_.push_back(h_Y_mom_slice);

    TH2D *h_thetaY_vs_time_mod_slice = new TH2D(("ThetaY_vs_Time_Modulo_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 149.2 ns", 29, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
    thetaY_vs_time_mod_slices_.push_back(h_thetaY_vs_time_mod_slice);

    TH2D *h_thetaY_vs_time_mod_20ns_slice = new TH2D(("ThetaY_vs_Time_Modulo_20ns_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 20 ns", 174, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
    thetaY_vs_time_mod_20ns_slices_.push_back(h_thetaY_vs_time_mod_20ns_slice);

    TH2D *h_thetaY_vs_time_mod_50ns_slice = new TH2D(("ThetaY_vs_Time_Modulo_50ns_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 50 ns", 87, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
    thetaY_vs_time_mod_50ns_slices_.push_back(h_thetaY_vs_time_mod_50ns_slice);

    TH2D *h_thetaY_vs_time_mod_long_slice = new TH2D(("ThetaY_vs_Time_Modulo_Long_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 149.2 ns", 29*moduloMultiple, 0, g2Period*moduloMultiple, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
    thetaY_vs_time_mod_long_slices_.push_back(h_thetaY_vs_time_mod_long_slice);

    TH2D *h_thetaY_vs_time_mod_long_20ns_slice = new TH2D(("ThetaY_vs_Time_Modulo_Long_20ns_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 20 ns", 174*moduloMultiple, 0, g2Period*moduloMultiple, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
    thetaY_vs_time_mod_long_20ns_slices_.push_back(h_thetaY_vs_time_mod_long_20ns_slice);

    TH2D *h_thetaY_vs_time_mod_long_50ns_slice = new TH2D(("ThetaY_vs_Time_Modulo_Long_50ns_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 50 ns", 87*moduloMultiple, 0, g2Period*moduloMultiple, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
    thetaY_vs_time_mod_long_50ns_slices_.push_back(h_thetaY_vs_time_mod_long_50ns_slice);

    TH1D *h_thetaY_mom_slice = new TH1D(("ThetaY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";#theta_{y} [mrad];Tracks",  500, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
    thetaY_mom_slices_.push_back(h_thetaY_mom_slice);

    TH2D *thetaY_vs_Y_mom_slice = new TH2D(("ThetaY_vs_Y_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay y-position [mm];#theta_{y} [mrad]", 24, -60, 60, 40, -100, 100);
    thetaY_vs_Y_mom_slices_.push_back(thetaY_vs_Y_mom_slice);

    acceptanceHists_.push_back((TH2D*)acceptanceFile->Get(("AcceptanceWeighting/MomBins/"+stn+"_WeightMap_"+to_string(lo)+"_"+to_string(hi)).c_str()));

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

    // Get variables
    double time = br.posiInitTime * 1e-3; // us
    time = RandomisedTime(rand, time); // randomise out the FR

    // Positron world momentum 
    TVector3 eMom(br.posiInitPX, br.posiInitPY, br.posiInitPZ); 
    TVector3 ePos(br.posiInitPosX, br.posiInitPosY, br.posiInitPosZ);
    TVector3 muPol(br.muDecayPolX, br.muDecayPolY, br.muDecayPolZ);

    // Positron momentum in the lab frame
    // double p_world = br.posiInitP;

    double g2ModTime = ModTime(time);
    double longModTime = ModTime(time, moduloMultiple);

    double y = ePos.Y(); 

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

    double px = eMom.X();
    double py = eMom.Y();
    double pz = eMom.Z();
    double pT = sqrt( pow(px, 2) + pow(py, 2) );
    double p = eMom.Mag();

    double theta_y = asin(py/p);

    double alpha = muPol.Angle(eMom);

    ////////////////////////////////////////////////////////////////////////////////////////////////

    // convert into mrad (always forget to do this so I'm putting it here)
    theta_y = theta_y * 1e3;

    // Correct offset (no time dependance here)
    if(verticalOffsetCorrection) {

      double theta_y_offset = verticalOffsetHist->GetBinContent(verticalOffsetHist->FindBin(p));
      theta_y = theta_y - theta_y_offset;

    }
    
    // Apply momentum dependant acceptance weighting
     // cout<<"\ntheta_y (unweighted) = "<<theta_y<<endl;
    
    double acceptanceWeighting = 1.0;

    if(acceptanceCorr) {
      
      // Slice momentum 
      for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

        int lo = 0 + i_slice*step; 
        int hi = step + i_slice*step;

        if(p >= double(lo) && p < double(hi)) { 

          acceptanceWeighting = AcceptanceWeightedAngle(acceptanceHists_.at(i_slice), theta_y, y);

        }

      }

    }


/*
          TH2D *acceptanceHist = (TH2D*)acceptanceFile->Get(("AcceptanceWeightings/"+to_string(lo)+"_"+to_string(hi)+"_MeV/"+stn+"_AcceptanceMap").c_str());

          if(acceptanceHist==0) acceptanceHist = acceptanceHistAllMom;

          double binCentre = (lo+hi)/2;
          double new_theta_y = 0; 

          if(p<binCentre) { // Get next lowest bin 

            TH2D *acceptanceHistLo = (TH2D*)acceptanceFile->Get(("AcceptanceWeightings/"+to_string(lo-step)+"_"+to_string(hi-step)+"_MeV/"+stn+"_AcceptanceMap").c_str());
            new_theta_y = AcceptanceWeightedAngleWithInterpolation(acceptanceHist, acceptanceHistLo, theta_y, y);
          
            //cout<<"---> LO: \np = "<<p<<"\nbin centre = "<<binCentre<<"\nnew theta_y = "<<new_theta_y<<endl;


          } else if(p>=binCentre) { // Get next highest bin 

            TH2D *acceptanceHistHi = (TH2D*)acceptanceFile->Get(("AcceptanceWeightings/"+to_string(lo+step)+"_"+to_string(hi+step)+"_MeV/"+stn+"_AcceptanceMap").c_str());
            new_theta_y = AcceptanceWeightedAngleWithInterpolation(acceptanceHist, acceptanceHistHi, theta_y, y);

            //cout<<"---> HI: \np = "<<p<<"\nbin centre = "<<binCentre<<"\nnew theta_y = "<<new_theta_y<<endl;

          }

          if(!isnan(new_theta_y)) {

            theta_y = new_theta_y;

            //cout<<"---> RESULT:\nnew theta_y = "<<new_theta_y<<endl;
        
          } */

       // }


    //cout<<"theta_y (weighted) = "<<theta_y<<endl;

    // Time cuts
    if(timeCuts && time < g2Period*7) continue; 

    decayZ_vs_decayX->Fill(ePos.X(), ePos.Z());

    if(p > 1700) { //  && p < pmax*momBoostFactor) {
      wiggle->Fill(time);
      wiggle_mod->Fill(g2ModTime);
      if(time>8*g2Period) wiggle_mod_long->Fill(longModTime);
    } 

    momY->Fill(py);
    momX->Fill(px);
    momZ->Fill(pz);

    // EDM cuts
    if(momCuts && p > pLo*momBoostFactor && p < pHi*momBoostFactor) { 

      momentum->Fill(p);
      thetaY->Fill(theta_y, acceptanceWeighting);
      thetaY_vs_time->Fill(time, theta_y, acceptanceWeighting);
      thetaY_vs_time_20ns->Fill(time, theta_y, acceptanceWeighting);
      thetaY_vs_time_50ns->Fill(time, theta_y, acceptanceWeighting);
      thetaY_vs_time_mod->Fill(g2ModTime, theta_y, acceptanceWeighting);
      thetaY_vs_time_mod_50ns->Fill(g2ModTime, theta_y, acceptanceWeighting);
      thetaY_vs_momentum->Fill(p, theta_y, acceptanceWeighting);
      thetaY_vs_Y->Fill(y, theta_y, acceptanceWeighting);

      if(time > 8*g2Period) {
        thetaY_vs_time_mod_long->Fill(longModTime, theta_y, acceptanceWeighting);
        thetaY_vs_time_mod_long_20ns->Fill(longModTime, theta_y, acceptanceWeighting); 
        thetaY_vs_time_mod_long_50ns->Fill(longModTime, theta_y, acceptanceWeighting);    
      }  

    } else if(!momCuts) { 

      momentum->Fill(p);
      thetaY->Fill(theta_y, acceptanceWeighting);
      thetaY_vs_time->Fill(time, theta_y, acceptanceWeighting);
      thetaY_vs_time_20ns->Fill(time, theta_y, acceptanceWeighting);
      thetaY_vs_time_50ns->Fill(time, theta_y, acceptanceWeighting);
      thetaY_vs_time_mod->Fill(g2ModTime, theta_y, acceptanceWeighting);
      thetaY_vs_time_mod_50ns->Fill(g2ModTime, theta_y, acceptanceWeighting);
      thetaY_vs_momentum->Fill(p, theta_y, acceptanceWeighting);
      thetaY_vs_Y->Fill(y, theta_y, acceptanceWeighting);

      if(time > 8*g2Period) {
        thetaY_vs_time_mod_long->Fill(longModTime, theta_y, acceptanceWeighting);
        thetaY_vs_time_mod_long_20ns->Fill(longModTime, theta_y, acceptanceWeighting); 
        thetaY_vs_time_mod_long_50ns->Fill(longModTime, theta_y, acceptanceWeighting);    
      }  
 

    }

    // Slice momentum 
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = 0 + i_slice*step; 
      int hi = step + i_slice*step;

      if(p >= double(lo) && p < double(hi)) { 

        thetaY_vs_time_mod_slices_.at(i_slice)->Fill(g2ModTime, theta_y, acceptanceWeighting);
        thetaY_vs_time_mod_20ns_slices_.at(i_slice)->Fill(g2ModTime, theta_y, acceptanceWeighting);
        thetaY_vs_time_mod_50ns_slices_.at(i_slice)->Fill(g2ModTime, theta_y, acceptanceWeighting);
      
        if(time > 8*g2Period) {
          thetaY_vs_time_mod_long_slices_.at(i_slice)->Fill(longModTime, theta_y, acceptanceWeighting);
          thetaY_vs_time_mod_long_20ns_slices_.at(i_slice)->Fill(longModTime, theta_y, acceptanceWeighting);
          thetaY_vs_time_mod_long_50ns_slices_.at(i_slice)->Fill(longModTime, theta_y, acceptanceWeighting);    
        }  

        // Other scans 
        thetaY_mom_slices_.at(i_slice)->Fill(theta_y, acceptanceWeighting);
        Y_mom_slices_.at(i_slice)->Fill(y);
        pY_mom_slices_.at(i_slice)->Fill(py);
        p_mom_slices_.at(i_slice)->Fill(p);
        thetaY_vs_Y_mom_slices_.at(i_slice)->Fill(y, theta_y, acceptanceWeighting);

      }

    }
 
    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }

    // Get max muon angle
    if(br.muDecayPolY > muAngleMax) muAngleMax = br.muDecayPolY;
         
  }

  cout<<"Muon max angle:\t"<<muAngleMax<<" radians"<<endl;
  // Write to output
  // Set output directory
  output->mkdir("SimultaneousAnalysis"); output->mkdir("MomentumBinnedAnalysis");

  output->cd("SimultaneousAnalysis");

  momentum->Write();
  wiggle->Write();
  wiggle_mod->Write();
  wiggle_mod_long->Write();
  thetaY->Write();
  thetaY_vs_time->Write();
  thetaY_vs_momentum->Write();
  thetaY_vs_time_20ns->Write(); 
  thetaY_vs_time_50ns->Write(); 
  thetaY_vs_time_mod->Write();
  thetaY_vs_time_mod_20ns->Write();
  thetaY_vs_time_mod_50ns->Write();
  thetaY_vs_time_mod_long->Write();
  thetaY_vs_time_mod_long_20ns->Write();
  thetaY_vs_time_mod_long_50ns->Write();
  decayZ_vs_decayX->Write();
  momX->Write();
  momY->Write();
  momZ->Write();
  thetaY_vs_Y->Write();

   for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) {

    output->cd("MomentumBinnedAnalysis"); 

    thetaY_vs_time_mod_slices_.at(i_slice)->Write();
    thetaY_vs_time_mod_20ns_slices_.at(i_slice)->Write();
    thetaY_vs_time_mod_50ns_slices_.at(i_slice)->Write();
    thetaY_vs_time_mod_long_slices_.at(i_slice)->Write();
    thetaY_vs_time_mod_long_20ns_slices_.at(i_slice)->Write();
    thetaY_vs_time_mod_long_50ns_slices_.at(i_slice)->Write();
    thetaY_mom_slices_.at(i_slice)->Write();
    Y_mom_slices_.at(i_slice)->Write();
    pY_mom_slices_.at(i_slice)->Write();
    p_mom_slices_.at(i_slice)->Write();
    thetaY_vs_Y_mom_slices_.at(i_slice)->Write();
      
  }

  verticalOffsetFile->Close();
  acceptanceFile->Close();

  return;

}


int main(int argc, char *argv[]) {

  bool boost = false;
  bool momCuts = true;//false; 
  bool timeCuts = true; 

  bool verticalOffsetCorr = false;
  bool acceptanceCorr = false;
  
  string inFileName = argv[1]; // "/pnfs/GM2/persistent/EDM/MC/dMu/TruthNTup/truthTrees.000.root";//
  string outFileName = argv[2]; //"tmp.root";
  //string config = argv[3]; // "Run-1a_250MeV_BQ";
  string stn = argv[3];

  string treeName = "phaseAnalyzer/g2phase";
  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName.c_str());
  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;

  // Book output
  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");

  // Fill histograms
  // RunNormally(TTree *tree, TFile *output, bool momCuts, bool timeCuts, bool boost, bool verticalOffsetCorrection, bool acceptanceCorr, string stn) {

  Run(tree, fout, momCuts, timeCuts, boost, verticalOffsetCorr, acceptanceCorr, stn);
  //RunWithDataDrivenAcceptanceCorrection(tree, fout, config, momCuts, timeCuts, boost, verticalOffsetCorr, acceptanceCorr);
  
  fout->Close();
  fin->Close();

  cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;

  return 0;
}

