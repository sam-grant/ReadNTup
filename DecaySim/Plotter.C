/*

Samuel Grant

Produce N vs t_mod and theta_y vs t_mod histograms for "all decays". This is the first step in the EDM analysis chain. 

Output format (defined in submit.sh): edmPlots_<recoType>_<refFrame>_<binWidth>_<qualityInfo>_<corrInfo>.root

recoType: allDecays
refFrame: LAB (lab frame) or MRF (muon rest frame)
binWidth: 250MeV (settled for Run-1)
qualityInfo: noQ (no beam vertex quality cuts) or Q (beam vertex quality cuts) 
corrInfo: randCorr_vertCorr_accWeight (time randomsiation, vertical angle offset correction, acceptance weighting -- left blank for no corrections) 

Notes: 
  * For Run-1 we only perform acceptance weighting on allDecays;
  * No vertical angle correction needed for allDecays.

*/

#include "../Common/Utils.h"
#include "Plotter.h"

// Lower and upper momentum analysis cuts.
double pLo = 1000; 
double pHi = 2500;

// Acceptance weighting using Delaunay interpolation 
// https://root.cern.ch/doc/master/classTGraph2D.html#a0dfb623f2a9f55c98ebe323384cf3f0d
double AcceptanceWeighting(TGraph2D *map, double y, double theta_y) {

  double weighting = map->Interpolate(y, theta_y);

  // Return 0 for NaN bins
  if(isnan(weighting)) {
    return 0;
  }

  return weighting;

}

// Run plotting
void Run(TTree *tree, TFile *output, bool momCuts, bool timeCuts, bool boost, bool randCorr, bool vertCorr, bool accWeight, string stn) {

  // CORRECTION HISTOGRAMS NEED TO BE REGENERATED AND UPDATED 

  // Get vertical offset correction histograms 
  TString verticalOffsetFileName = "correctionHists/verticalOffsetHists_allDecays_WORLD_250MeV_AQ.root";
  TFile *verticalOffsetFile = TFile::Open(verticalOffsetFileName);
  TH1D *verticalOffsetHist = (TH1D*)verticalOffsetFile->Get("VerticalOffsetHists/ThetaY_vs_p");

  TString acceptanceFileName = "correctionHists/acceptanceMaps.thetaYvsY.truth.root"; // "correctionHists/acceptanceWeightingPlots.thetaYvsY.truth.root";
  TFile *acceptanceFile = TFile::Open(acceptanceFileName);

  cout<<"---> Opened correction files"<<endl;

  // Containers for acceptance weighting maps
  vector<TH2D*> acceptanceHists_; 
  vector<TGraph2D*> acceptanceGraphs_; 

  // Somewhat arbitrary, it just scales the maximum vertical angle
  double boostFactor = 5e3*(1/gmagic);
  double momBoostFactor = 1.;

  // Adjust if we're using the MRF
  if(boost) { 
    boostFactor = 1.0e3;
    momBoostFactor = (1/(2*gmagic));
  }

  // Book histograms

  cout<<"---> Booking histograms"<<endl;
  TH1D *momentum_noCuts = new TH1D("Momentum_NoCuts", ";Momentum [MeV];Tracks", int(pmax), 0, pmax*momBoostFactor); // Prior to any cuts
  TH1D *momentum = new TH1D("Momentum", ";Momentum [MeV];Tracks", int(pmax), 0, pmax*momBoostFactor); 
  TH1D *momY = new TH1D("MomentumY", ";Momentum Y [MeV];Tracks", 1000, -60, 60); 
  TH2D *decayZ_vs_decayX = new TH2D("DecayZ_vs_DecayX", ";Decay position X [mm];Decay position Z [mm]", 800, -8000, 8000, 800, -8000, 8000);
  TH1D *wiggle = new TH1D("Wiggle", ";Decay time [#mus];Tracks", 2700, 0, 2700*T_c); 
  TH1D *wiggle_mod_A = new TH1D("Wiggle_Modulo_A", ";t_{g#minus2}^{mod} [#mus];Decays / 149.2 ns", 29, 0, g2Period); // Using T-method momentum cut, to obtain phase
  TH1D *wiggle_mod_B = new TH1D("Wiggle_Modulo_B", ";t_{g#minus2}^{mod} [#mus];Decays / 149.2 ns", 29, 0, g2Period); // Using EDM momentum cuts, for theta_y fit demoninator 
  TH1D *thetaY = new TH1D("ThetaY", ";#theta_{y} [mrad];Decays", 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_time = new TH2D("ThetaY_vs_Time", ";Decay time [#mus]; #theta_{y} [mrad] / 149.2 ns ", 2700, 0, 2700*T_c, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_time_mod = new TH2D("ThetaY_vs_Time_Modulo", ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 149.2 ns", 29, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_Y = new TH2D("ThetaY_vs_Y", ";Decay y-position [mm];#theta_{y} [mrad]", 24, -60, 60, 40, -100, 100);
  TH2D *thetaY_vs_momentum_noCuts = new TH2D("ThetaY_vs_Momentum_NoCuts", ";Decay momentum [MeV]; #theta_{y} [mrad] / 10 MeV ", 312, 0, 3120, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
  TH2D *thetaY_vs_momentum = new TH2D("ThetaY_vs_Momentum", ";Decay momentum [MeV]; #theta_{y} [mrad] / 10 MeV ", 312, 0, 3120, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);

  // Momentum scans
  vector<TH1D*> thetaY_mom_slices_;
  vector<TH1D*> wiggle_mod_mom_slices_;
  vector<TH1D*> Y_mom_slices_;
  vector<TH1D*> pY_mom_slices_;
  vector<TH1D*> p_mom_slices_;
  vector<TH2D*> thetaY_vs_Y_mom_slices_;
  vector<TH2D*> thetaY_vs_time_mod_slices_;
  
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

    TH1D *h_wiggle_mod_mom_slice = new TH1D(("Wiggle_Modulo_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus];Tracks / 149.2 ns", 29, 0, g2Period); 
    wiggle_mod_mom_slices_.push_back(h_wiggle_mod_mom_slice);

    TH2D *h_thetaY_vs_time_mod_slice = new TH2D(("ThetaY_vs_Time_Modulo_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 149.2 ns", 29, 0, g2Period, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
    thetaY_vs_time_mod_slices_.push_back(h_thetaY_vs_time_mod_slice);

    TH1D *h_thetaY_mom_slice = new TH1D(("ThetaY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";#theta_{y} [mrad];Tracks",  500, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_mom_slices_.push_back(h_thetaY_mom_slice);

    TH2D *thetaY_vs_Y_mom_slice = new TH2D(("ThetaY_vs_Y_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay y-position [mm];#theta_{y} [mrad]", 24, -60, 60, 40, -100, 100);
    thetaY_vs_Y_mom_slices_.push_back(thetaY_vs_Y_mom_slice);

    acceptanceHists_.push_back((TH2D*)acceptanceFile->Get(("AcceptanceWeighting/MomBins/"+stn+"_WeightMapY_"+to_string(lo)+"_"+to_string(hi)).c_str()));
    acceptanceGraphs_.push_back((TGraph2D*)acceptanceFile->Get(("AcceptanceWeighting/MomBins/"+stn+"_WeightGraphY_"+to_string(lo)+"_"+to_string(hi)).c_str()));

  }

  // Get branches (using header file)
  InitBranches br(tree);

  double targetPerc = 0;
  int64_t nEntries = tree->GetEntries();

  double muAngleMax = 0;

  int64_t counter = 0;

  // For FR randomisation
  TRandom3 *rand = new TRandom3(12345);

  cout<<"---> Event loop"<<endl;
   
  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    // Get variables
    double time = br.posiInitTime * 1e-3; // us
    if(randCorr) time = RandomisedTime(rand, time); // randomise out the FR

    // Positron world momentum 
    TVector3 eMom(br.posiInitPX, br.posiInitPY, br.posiInitPZ); 
    TVector3 ePos(br.posiInitPosX, br.posiInitPosY, br.posiInitPosZ);
    TVector3 muPol(br.muDecayPolX, br.muDecayPolY, br.muDecayPolZ);

    double g2ModTime = ModTime(time);
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

    // Define the vertical angle
    double py = eMom.Y();
    double p = eMom.Mag();
    double theta_y = asin(py/p);

    // Convert from rad to mrad 
    theta_y = theta_y * 1e3;

    // Useful plots to have without cuts
    momentum_noCuts->Fill(p);
    thetaY_vs_momentum_noCuts->Fill(p, theta_y);

    // Correct offset (no time dependance here)
    if(vertCorr) {
      double theta_y_offset = verticalOffsetHist->GetBinContent(verticalOffsetHist->FindBin(p));
      theta_y = theta_y - theta_y_offset;
    }
    
    // Obtain acceptance weight for this theta_y 
    double acceptanceWeighting = 1.0;
    if(accWeight) {
      
      // Slice momentum 
      for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

        int lo = 0 + i_slice*step; 
        int hi = step + i_slice*step;

        if(p >= double(lo) && p < double(hi)) { 

          acceptanceWeighting = AcceptanceWeighting(acceptanceGraphs_.at(i_slice), y, theta_y);

        }

      }

    }

    // Early time cut
    if(timeCuts && time < g2Period*7) continue; 

    decayZ_vs_decayX->Fill(ePos.X(), ePos.Z());

    // T-method energy cut 
    if(p > 1700) { 
      wiggle->Fill(time);
      wiggle_mod_A->Fill(g2ModTime);
    } 

    // EDM cuts
    if(momCuts && p > pLo*momBoostFactor && p < pHi*momBoostFactor) { 
      momentum->Fill(p);
      momY->Fill(py);
      thetaY->Fill(theta_y, acceptanceWeighting);
      thetaY_vs_time->Fill(time, theta_y, acceptanceWeighting);
      thetaY_vs_time_mod->Fill(g2ModTime, theta_y, acceptanceWeighting);
      thetaY_vs_momentum->Fill(p, theta_y, acceptanceWeighting);
      thetaY_vs_Y->Fill(y, theta_y, acceptanceWeighting);
      wiggle_mod_B->Fill(g2ModTime);
    } else if (!momCuts) { 
      momentum->Fill(p);
      momY->Fill(py);
      thetaY->Fill(theta_y, acceptanceWeighting);
      thetaY_vs_time->Fill(time, theta_y, acceptanceWeighting);
      thetaY_vs_time_mod->Fill(g2ModTime, theta_y, acceptanceWeighting);
      thetaY_vs_momentum->Fill(p, theta_y, acceptanceWeighting);
      thetaY_vs_Y->Fill(y, theta_y, acceptanceWeighting);
      wiggle_mod_B->Fill(g2ModTime);
    }

    // Slice momentum 
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = 0 + i_slice*step; 
      int hi = step + i_slice*step;

      if(p >= double(lo) && p < double(hi)) { 

        wiggle_mod_mom_slices_.at(i_slice)->Fill(g2ModTime);
        thetaY_vs_time_mod_slices_.at(i_slice)->Fill(g2ModTime, theta_y, acceptanceWeighting);

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

    // Get max muon angle (nice sanity check)
    if(br.muDecayPolY > muAngleMax) muAngleMax = br.muDecayPolY;
         
  }

  cout<<"---> Sanity check: muon max angle:\t"<<muAngleMax<<" radians"<<endl;

  // Write to output
  cout<<"---> Writing histograms"<<endl;
  
  // Setup directories
  output->mkdir("SimultaneousAnalysis"); output->mkdir("MomentumBinnedAnalysis");
  output->cd("SimultaneousAnalysis");

  // Write histograms
  momentum_noCuts->Write();
  momentum->Write();
  wiggle->Write();
  wiggle_mod_A->Write();
  wiggle_mod_B->Write();
  thetaY->Write();
  thetaY_vs_time->Write();
  thetaY_vs_momentum->Write();
  thetaY_vs_momentum_noCuts->Write();
  thetaY_vs_time_mod->Write();
  decayZ_vs_decayX->Write();
  momY->Write();
  thetaY_vs_Y->Write();

   for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) {

    output->cd("MomentumBinnedAnalysis"); 

    wiggle_mod_mom_slices_.at(i_slice)->Write();
    thetaY_vs_time_mod_slices_.at(i_slice)->Write();
    thetaY_mom_slices_.at(i_slice)->Write();
    Y_mom_slices_.at(i_slice)->Write();
    pY_mom_slices_.at(i_slice)->Write();
    p_mom_slices_.at(i_slice)->Write();
    thetaY_vs_Y_mom_slices_.at(i_slice)->Write();
      
  }

  // Close input files
  verticalOffsetFile->Close();
  acceptanceFile->Close();

  return;

}

/*

Run test:
./Plotter.out /pnfs/GM2/persistent/EDM/MC/dMu/Trees/AllDecays/5.4e-18/decaySimTrees_000.root output/test.root None

Run full:
. submit.sh 3 None

*/


int main(int argc, char *argv[]) {

  // Get inputs from command line
  string inFileName = argv[1]; 
  string outFileName = argv[2];
  string stn = argv[3]; // Sets acceptance weighting per station

  bool boost = false;
  bool momCuts = true; 
  bool timeCuts = true; 
  bool randCorr = true; // FR
  bool vertCorr = false;
  bool accWeight = false; 

  if(stn=="S12" || stn=="S18" || stn=="S12S18") accWeight = true;

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName.c_str());

  // Get tree
  string treeName = "phaseAnalyzer/g2phase";
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;

  // Book output
  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");

  // Fill histograms
  Run(tree, fout, momCuts, timeCuts, boost, randCorr, vertCorr, accWeight, stn);
  
  // Close files
  fout->Close();
  fin->Close();

  cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;

  return 0;

}

