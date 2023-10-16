/*

Samuel Grant

Produce histograms for acceptance maps. 

Compare tracks and decays, using trees formed from Joe's jprice_mc_digit_trackAcceptance sample, 
/pnfs/GM2/persistent/EDM/MC/TrackerAcceptance/jprice_mc_digit_trackAcceptance.

Produce histograms of theta_y vs y (acceptance maps) over a specified momentum range and in momentum bins.

Also, produce histograms in vertical position bins, for illustration purposes.

This is the third step in the EDM analysis chain.

*/

#include "../Common/Utils.h"
#include "Plotter.h"

void Run(TTree *treeAllDecays, TTree *treeTracks, TFile *output, bool quality, bool momCuts, bool timeCuts, bool truth, double pLo = 750, double pHi = 2750) {

  string stns[] = {"S12S18", "S12", "S18"}; 

  int n_stn = sizeof(stns)/sizeof(stns[0]);

  // Bins and ranges

  // Vertical angle, 5 mrad bins
  double thetaY_min = -100;
  double thetaY_max = 100;
  int thetaY_nBins = (thetaY_max - thetaY_min) / 5; 

  // Vertical position, 5 mm bins
  double y_min = -60;
  double y_max = 60;
  int y_nBins = (y_max - y_min) / 5; 

  // Radial position 
  double r_min = -60;
  double r_max = 60;
  int r_nBins = (r_max - r_min) / 5;

  // Ring azimuthal angle
  double phi_min = 0;
  double phi_max = TMath::TwoPi();
  int phi_nBins = (phi_max - phi_min) / 5;

  // Book histograms 

  // Decays
  TH1D *p_decays = new TH1D("Momentum", ";Track momentum [MeV];Decays", int(pmax), 0, pmax);  
  TH1D *Y_decays = new TH1D("Y", ";Decay y-position [mm];Decays", y_nBins, y_min, y_max);
  TH1D *R_decays = new TH1D("R", ";Decay radial position [mm];Decays", r_nBins, r_min, r_max);
  TH1D *Phi_decays = new TH1D("Phi", ";Decay azimuthal angle [rad];Decays", phi_nBins, 0, TMath::TwoPi());
  TH1D *thetaY_decays = new TH1D("ThetaY", "#theta_{y} [mrad];Decays", thetaY_nBins, thetaY_min, thetaY_max); 

  TH2D *thetaY_vs_Y_decays = new TH2D("ThetaY_vs_Y", ";Decay vertical position [mm];#theta_{y} [mrad]", y_nBins, y_min, y_max, thetaY_nBins, thetaY_min, thetaY_max);
  TH2D *thetaY_vs_R_decays = new TH2D("ThetaY_vs_R", ";Decay radial position [mm];#theta_{y} [mrad]", r_nBins, r_min, r_max, thetaY_nBins, thetaY_min, thetaY_max);
  TH2D *thetaY_vs_Phi_decays = new TH2D("ThetaY_vs_Phi", ";Decay azimuthal angle [rad];#theta_{y} [mrad]", phi_nBins, phi_min, phi_max, thetaY_nBins, thetaY_min, thetaY_max);
  TH2D *thetaY_vs_p_decays = new TH2D("ThetaY_vs_p", ";Momentum [MeV];#theta_{y} [mrad]", 48, 0, pmax, thetaY_nBins, thetaY_min, thetaY_max);//-TMath::Pi()*gmagic, TMath::Pi()*gmagic);
 
  // Tracks
  TH1D *p_tracks_[n_stn];  
  TH1D *Y_tracks_[n_stn]; 
  TH1D *R_tracks_[n_stn]; 
  TH1D *Phi_tracks_[n_stn]; 
  TH1D *thetaY_tracks_[n_stn]; 

  TH2D *thetaY_vs_Y_tracks_[n_stn]; 
  TH2D *thetaY_vs_R_tracks_[n_stn];
  TH2D *thetaY_vs_Phi_tracks_[n_stn];
  TH2D *thetaY_vs_p_tracks_[n_stn]; 

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    p_tracks_[i_stn] = new TH1D((stns[i_stn]+"_Momentum").c_str(), ";Track momentum [MeV];Tracks", int(pmax), 0, pmax); 
    Y_tracks_[i_stn] = new TH1D((stns[i_stn]+"_Y").c_str(), ";Decay y-position [mm];Tracks", y_nBins, y_min, y_max);
    R_tracks_[i_stn] = new TH1D((stns[i_stn]+"_R").c_str(), ";Decay radial position [mm];Tracks", r_nBins, r_min, r_max);
    Phi_tracks_[i_stn] = new TH1D((stns[i_stn]+"_Phi").c_str(), ";Decay azimuthal angle [rad];Tracks", phi_nBins, phi_min, phi_max);
    thetaY_tracks_[i_stn] = new TH1D((stns[i_stn]+"_ThetaY").c_str(), ";#theta_{y} [mrad];Tracks", thetaY_nBins, thetaY_min, thetaY_max);
    
    thetaY_vs_Y_tracks_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Y").c_str(), ";Decay y-position [mm]; #theta_{y} [mrad]", y_nBins, y_min, y_max, thetaY_nBins, thetaY_min, thetaY_max);
    thetaY_vs_R_tracks_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_R").c_str(), ";Decay radial position [mm];#theta_{y} [mrad]", r_nBins, r_min, r_max, thetaY_nBins, thetaY_min, thetaY_max);
    thetaY_vs_Phi_tracks_[i_stn] =  new TH2D((stns[i_stn]+"_ThetaY_vs_Phi").c_str(), ";Decay azimuthal angle [rad];#theta_{y} [mrad]", phi_nBins, phi_min, phi_max, thetaY_nBins, thetaY_min, thetaY_max); 
    thetaY_vs_p_tracks_[i_stn] =  new TH2D((stns[i_stn]+"_ThetaY_vs_p").c_str(), ";Track momentum [MeV];#theta_{y} [mrad]", 48, 0, pmax, thetaY_nBins, thetaY_min, thetaY_max); 

  }

  /*  

  Slice momentum in steps of 250 MeV from 0-3000 MeV 

  */

  int p_step = 250; 
  int p_nSlices = pmax/p_step;

  // Decays 
  vector<TH1D*> p_decays_p_slices_;
  vector<TH1D*> Y_decays_p_slices_;
  vector<TH1D*> thetaY_decays_p_slices_;
  vector<TH1D*> R_decays_p_slices_; 
  vector<TH1D*> Phi_decays_p_slices_; 

  vector<TH2D*> thetaY_vs_Y_decays_p_slices_;
  vector<TH2D*> thetaY_vs_R_decays_p_slices_;
  vector<TH2D*> thetaY_vs_Phi_decays_p_slices_;
  vector<TH2D*> thetaY_vs_p_decays_p_slices_;

  // Tracks
  vector<TH1D*> p_tracks_p_slices_[n_stn];
  vector<TH1D*> Y_tracks_p_slices_[n_stn];
  vector<TH1D*> R_tracks_p_slices_[n_stn];
  vector<TH1D*> Phi_tracks_p_slices_[n_stn];
  vector<TH1D*> thetaY_tracks_p_slices_[n_stn];

  vector<TH2D*> thetaY_vs_Y_tracks_p_slices_[n_stn];
  vector<TH2D*> thetaY_vs_R_tracks_p_slices_[n_stn];
  vector<TH2D*> thetaY_vs_Phi_tracks_p_slices_[n_stn];
  vector<TH2D*> thetaY_vs_p_tracks_p_slices_[n_stn];

  for ( int i_slice = 0; i_slice < p_nSlices; i_slice++ ) { 

    int lo = i_slice*p_step; 
    int hi = p_step + i_slice*p_step;

    std::string stepStr = to_string(lo)+"_"+to_string(hi);

    // Decays 
    p_decays_p_slices_.push_back(new TH1D(("Momentum_"+stepStr).c_str(), ";Momentum [MeV];Decays", int(pmax), 0, pmax)); 
    R_decays_p_slices_.push_back(new TH1D(("R_"+stepStr).c_str(), ";Decay radial position [mm];Tracks", r_nBins, r_min, r_max));
    Phi_decays_p_slices_.push_back(new TH1D(("Phi_"+stepStr).c_str(), ";Decay azimuthal angle [rad];Tracks", phi_nBins, phi_min, phi_max));
    Y_decays_p_slices_.push_back(new TH1D(("Y_"+stepStr).c_str(), ";Decay y-position [mm];Tracks", y_nBins, y_min, y_max));
    thetaY_decays_p_slices_.push_back(new TH1D(("ThetaY_"+stepStr).c_str(), ";#theta_{y} [mrad];Tracks", thetaY_nBins, thetaY_min, thetaY_max)); 
    thetaY_vs_Y_decays_p_slices_.push_back(new TH2D(("ThetaY_vs_Y_"+stepStr).c_str(), ";Decay y-position [mm]; #theta_{y} [mrad]", y_nBins, y_min, y_max, thetaY_nBins, thetaY_min, thetaY_max)); 
    thetaY_vs_R_decays_p_slices_.push_back(new TH2D(("ThetaY_vs_R_"+stepStr).c_str(), ";Decay radial position [mm];Tracks", r_nBins, r_min, r_max, thetaY_nBins, thetaY_min, thetaY_max)); 
    thetaY_vs_Phi_decays_p_slices_.push_back(new TH2D(("ThetaY_vs_Phi_"+stepStr).c_str(), ";Decay azimuthal angle [rad];Tracks", phi_nBins, phi_min, phi_max, thetaY_nBins, thetaY_min, thetaY_max));
    thetaY_vs_p_decays_p_slices_.push_back(new TH2D(("ThetaY_vs_p_"+stepStr).c_str(), ";Momentum [MeV];#theta_{y} [mrad]", 48, 0, pmax, thetaY_nBins, thetaY_min, thetaY_max));

    for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

      // Tracks
      p_tracks_p_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_Momentum_"+stepStr).c_str(), ";Track momentum [MeV];Tracks", int(pmax), 0, pmax)); 
      Y_tracks_p_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_Y_"+stepStr).c_str(), ";Decay y-position [mm];Tracks", y_nBins, y_min, y_max));
      R_tracks_p_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_R_"+stepStr).c_str(), ";Decay radial position [mm];Tracks", r_nBins, r_min, r_max));
      Phi_tracks_p_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_Phi_"+stepStr).c_str(), ";Decay azimuthal angle [rad];Tracks", phi_nBins, phi_min, phi_max));
      thetaY_tracks_p_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_ThetaY_"+stepStr).c_str(), ";#theta_{y} [mrad];Tracks", thetaY_nBins, thetaY_min, thetaY_max)); 
      thetaY_vs_Y_tracks_p_slices_[i_stn].push_back(new TH2D((stns[i_stn]+"_ThetaY_vs_Y_"+stepStr).c_str(), ";Decay y-position [mm]; #theta_{y} [mrad]", y_nBins, y_min, y_max, thetaY_nBins, thetaY_min, thetaY_max)); 
      thetaY_vs_p_tracks_p_slices_[i_stn].push_back(new TH2D((stns[i_stn]+"_ThetaY_vs_p_"+stepStr).c_str(), ";Track momentum [MeV];#theta_{y} [mrad]", 48, 0, pmax, thetaY_nBins, thetaY_min, thetaY_max));
      thetaY_vs_R_tracks_p_slices_[i_stn].push_back(new TH2D((stns[i_stn]+"_ThetaY_vs_R_"+stepStr).c_str(), ";Decay radial position [mm];Tracks", r_nBins, r_min, r_max, thetaY_nBins, thetaY_min, thetaY_max)); 
      thetaY_vs_Phi_tracks_p_slices_[i_stn].push_back(new TH2D((stns[i_stn]+"_ThetaY_vs_Phi_"+stepStr).c_str(), ";Decay azimuthal angle [rad];Tracks", phi_nBins, phi_min, phi_max, thetaY_nBins, thetaY_min, thetaY_max)); 

    }

  }

  /*  

  Slice y in steps of 5 mm between ±45 mm 

  */

  // Vertical position slices
  int y_step = 5; 
  int y_nSlices = 18;

  // decays 
  vector<TH1D*> p_decays_y_slices_;
  vector<TH1D*> thetaY_decays_y_slices_;
  vector<TH1D*> Y_decays_y_slices_;
  vector<TH1D*> R_decays_y_slices_; 
  vector<TH1D*> Phi_decays_y_slices_; 

  vector<TH2D*> thetaY_vs_Y_decays_y_slices_;
  vector<TH2D*> thetaY_vs_R_decays_y_slices_;
  vector<TH2D*> thetaY_vs_Phi_decays_y_slices_;
  vector<TH2D*> thetaY_vs_p_decays_y_slices_;

  // tracks
  vector<TH1D*> p_tracks_y_slices_[n_stn];
  vector<TH1D*> Y_tracks_y_slices_[n_stn];
  vector<TH1D*> R_tracks_y_slices_[n_stn]; 
  vector<TH1D*> Phi_tracks_y_slices_[n_stn]; 
  vector<TH1D*> thetaY_tracks_y_slices_[n_stn];

  vector<TH2D*> thetaY_vs_Y_tracks_y_slices_[n_stn];
  vector<TH2D*> thetaY_vs_R_tracks_y_slices_[n_stn];
  vector<TH2D*> thetaY_vs_Phi_tracks_y_slices_[n_stn];
  vector<TH2D*> thetaY_vs_p_tracks_y_slices_[n_stn];

  for ( int i_slice = 0; i_slice < y_nSlices; i_slice++ ) { 

    int lo = -45 + i_slice*y_step; 
    int hi = y_step + lo;

    std::string stepStr = to_string(lo)+"_"+to_string(hi);

    p_decays_y_slices_.push_back(new TH1D(("Momentum_"+stepStr).c_str(), ";Momentum [MeV];Decays", int(pmax), 0, pmax));
    Y_decays_y_slices_.push_back(new TH1D(("Y_"+stepStr).c_str(), ";Decay y-position [mm];Tracks", y_nBins, y_min, y_max));
    R_decays_y_slices_.push_back(new TH1D(("R_"+stepStr).c_str(), ";Decay radial position [mm];Tracks", r_nBins, r_min, r_max));
    Phi_decays_y_slices_.push_back(new TH1D(("Phi_"+stepStr).c_str(), ";Decay azimuthal angle [rad];Tracks", phi_nBins, phi_min, phi_max));
    thetaY_decays_y_slices_.push_back(new TH1D(("ThetaY_"+stepStr).c_str(), ";#theta_{y} [mrad];Tracks", thetaY_nBins, thetaY_min, thetaY_max)); 
    thetaY_vs_Y_decays_y_slices_.push_back(new TH2D(("ThetaY_vs_Y_"+stepStr).c_str(), ";Decay y-position [mm]; #theta_{y} [mrad]", y_nBins, y_min, y_max, thetaY_nBins, thetaY_min, thetaY_max)); 
    thetaY_vs_R_decays_y_slices_.push_back(new TH2D(("ThetaY_vs_R_"+stepStr).c_str(), ";Decay radial position [mm];Tracks", r_nBins, r_min, r_max, thetaY_nBins, thetaY_min, thetaY_max)); 
    thetaY_vs_Phi_decays_y_slices_.push_back(new TH2D(("ThetaY_vs_Phi_"+stepStr).c_str(), ";Decay azimuthal angle [rad];Tracks", phi_nBins, phi_min, phi_max, thetaY_nBins, thetaY_min, thetaY_max));
    thetaY_vs_p_decays_y_slices_.push_back(new TH2D(("ThetaY_vs_p_"+stepStr).c_str(), ";Momentum [MeV]; #theta_{y} [mrad]", 48, 0, pmax, thetaY_nBins, thetaY_min, thetaY_max));

    for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

      p_tracks_y_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_Momentum_"+stepStr).c_str(), ";Track momentum [MeV];Tracks", int(pmax), 0, pmax)); 
      Y_tracks_y_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_Y_"+stepStr).c_str(), ";Decay y-position [mm];Tracks", y_nBins, y_min, y_max));
      R_tracks_y_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_R_"+stepStr).c_str(), ";Decay radial position [mm];Tracks", r_nBins, r_min, r_max));
      Phi_tracks_y_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_Phi_"+stepStr).c_str(), ";Decay azimuthal angle [rad];Tracks", phi_nBins, phi_min, phi_max));
      thetaY_tracks_y_slices_[i_stn].push_back(new TH1D((stns[i_stn]+"_ThetaY_"+stepStr).c_str(), ";#theta_{y} [mrad];Tracks", thetaY_nBins, thetaY_min, thetaY_max)); 
      thetaY_vs_Y_tracks_y_slices_[i_stn].push_back(new TH2D((stns[i_stn]+"_ThetaY_vs_Y_"+stepStr).c_str(), ";Decay y-position [mm]; #theta_{y} [mrad]", y_nBins, y_min, y_max, thetaY_nBins, thetaY_min, thetaY_max)); 
      thetaY_vs_R_tracks_y_slices_[i_stn].push_back(new TH2D((stns[i_stn]+"_ThetaY_vs_R_"+stepStr).c_str(), ";Decay radial position [mm];Tracks", r_nBins, r_min, r_max, thetaY_nBins, thetaY_min, thetaY_max)); 
      thetaY_vs_Phi_tracks_y_slices_[i_stn].push_back(new TH2D((stns[i_stn]+"_ThetaY_vs_Phi_"+stepStr).c_str(), ";Decay azimuthal angle [rad];Tracks", phi_nBins, phi_min, phi_max, thetaY_nBins, thetaY_min, thetaY_max)); 
      thetaY_vs_p_tracks_y_slices_[i_stn].push_back(new TH2D((stns[i_stn]+"_ThetaY_vs_p_"+stepStr).c_str(), ";Track momentum [MeV];#theta_{y} [mrad]", 48, 0, pmax, thetaY_nBins, thetaY_min, thetaY_max));

    }

  }

  // Fill histograms

  // Get decay branches
  InitBranchesAllDecays allDecays(treeAllDecays);

  cout<<"-----> Filling decay histograms"<<endl;

  for(int64_t entry = 0; entry < treeAllDecays->GetEntries(); entry++) {

    treeAllDecays->GetEntry(entry);

    double y = allDecays.posiInitPosY;
    double p = allDecays.posiInitP;
    double py = allDecays.posiInitPY;
    double theta_y = asin(py/p) * 1e3; // mrad
    double t = allDecays.posiInitTime;
    double r = sqrt(pow(allDecays.posiInitPosX,2) + pow(allDecays.posiInitPosZ,2)) - R_magic;
    // RingAngle is angle from x axis, from 0 to 2pi
    double phi = atan2(allDecays.posiInitPosZ, allDecays.posiInitPosX);    
    if (phi  < 0) phi += TMath::TwoPi();

    // Time cuts
    if(timeCuts && t < g2Period*7) continue; 

    // Cut momentum according to Run() arguments
    if(momCuts && p > pLo && p < pHi) {

      p_decays->Fill(p);
      thetaY_decays->Fill(theta_y);
      Y_decays->Fill(y);
      R_decays->Fill(r);
      Phi_decays->Fill(phi);
      thetaY_vs_Y_decays->Fill(y, theta_y);
      thetaY_vs_R_decays->Fill(r, theta_y);
      thetaY_vs_Phi_decays->Fill(phi, theta_y);
      thetaY_vs_p_decays->Fill(p, theta_y);

    } else if(!momCuts) {

      p_decays->Fill(p);
      thetaY_decays->Fill(theta_y);
      Y_decays->Fill(y);
      R_decays->Fill(r);
      Phi_decays->Fill(phi);
      thetaY_vs_Y_decays->Fill(y, theta_y);
      thetaY_vs_R_decays->Fill(r, theta_y);
      thetaY_vs_Phi_decays->Fill(phi, theta_y);
      thetaY_vs_p_decays->Fill(p, theta_y);

    }

    // Slice momentum
    for ( int i_slice = 0; i_slice < p_nSlices; i_slice++ ) { 

      int lo = i_slice*p_step; 
      int hi = p_step + lo;

      if(p >= double(lo) && p < double(hi)) { 

        p_decays_p_slices_.at(i_slice)->Fill(p);
        thetaY_decays_p_slices_.at(i_slice)->Fill(theta_y);
        Y_decays_p_slices_.at(i_slice)->Fill(y);
        R_decays_p_slices_.at(i_slice)->Fill(r);
        Phi_decays_p_slices_.at(i_slice)->Fill(phi);
        thetaY_vs_Y_decays_p_slices_.at(i_slice)->Fill(y, theta_y);
        thetaY_vs_R_decays_p_slices_.at(i_slice)->Fill(r, theta_y);
        thetaY_vs_Phi_decays_p_slices_.at(i_slice)->Fill(phi, theta_y);
        thetaY_vs_p_decays_p_slices_.at(i_slice)->Fill(p, theta_y);

      }

    }

    // Slice y
    for ( int i_slice = 0; i_slice < y_nSlices; i_slice++ ) { 

      int lo = -45 + i_slice*y_step; 
      int hi = y_step + lo;

      if(y >= double(lo) && y < double(hi)) { 

        if(momCuts && p > pLo && p < pHi) {

          p_decays_y_slices_.at(i_slice)->Fill(p);
          thetaY_decays_y_slices_.at(i_slice)->Fill(theta_y);
          Y_decays_y_slices_.at(i_slice)->Fill(y);
          R_decays_y_slices_.at(i_slice)->Fill(r);
          Phi_decays_y_slices_.at(i_slice)->Fill(phi);
          thetaY_vs_Y_decays_y_slices_.at(i_slice)->Fill(y, theta_y);
          thetaY_vs_R_decays_y_slices_.at(i_slice)->Fill(r, theta_y);
          thetaY_vs_Phi_decays_y_slices_.at(i_slice)->Fill(phi, theta_y);
          thetaY_vs_p_decays_y_slices_.at(i_slice)->Fill(p, theta_y);

        } else if(!momCuts) { 

          p_decays_y_slices_.at(i_slice)->Fill(p);
          thetaY_decays_y_slices_.at(i_slice)->Fill(theta_y);
          Y_decays_y_slices_.at(i_slice)->Fill(y);
          R_decays_y_slices_.at(i_slice)->Fill(r);
          Phi_decays_y_slices_.at(i_slice)->Fill(phi);
          thetaY_vs_Y_decays_y_slices_.at(i_slice)->Fill(y, theta_y);
          thetaY_vs_R_decays_y_slices_.at(i_slice)->Fill(r, theta_y);
          thetaY_vs_Phi_decays_y_slices_.at(i_slice)->Fill(phi, theta_y);
          thetaY_vs_p_decays_y_slices_.at(i_slice)->Fill(p, theta_y);

        }

      }

    }

  }

  // Get track branches
  InitBranchesTracks tracks(treeTracks);

  cout<<"-----> Filling track histograms"<<endl;

  for(int64_t entry = 0; entry < treeTracks->GetEntries(); entry++) {

    treeTracks->GetEntry(entry);

    double y; double p; double py; double t; double r; double phi;

    if(truth) {
      y = tracks.trueVertexPosY;
      p = tracks.trueVertexMom;
      py = tracks.trueVertexMomY;
      t = tracks.trueVertexTime;
      r = sqrt(pow(tracks.trueVertexPosX,2)+pow(tracks.trueVertexPosZ,2)) - R_magic;
      phi = atan2(tracks.trueVertexPosZ, tracks.trueVertexPosX);   

    } else { 
      y = tracks.recoVertexPosY;
      p = tracks.recoVertexMom;
      py = -tracks.recoVertexMomY;
      t = tracks.recoVertexTime;
      r = sqrt(pow(tracks.recoVertexPosX,2)+pow(tracks.recoVertexPosZ,2)) - R_magic;
      phi = atan2(tracks.recoVertexPosZ, tracks.recoVertexPosX);  
    }

    if (phi  < 0) phi += TMath::TwoPi();

    // Time cuts
    if(timeCuts && t < g2Period*7) continue; 

    double theta_y = asin(py/p) * 1e3; 
    int stn = tracks.station;

    if(quality) {
      if(!tracks.passVertexQuality) continue;
      if(tracks.pValue < 0.05) continue;
      if(tracks.hitVolume) continue; 
    } 

    int stn_id = -1;
    if(stn==12) stn_id = 1;
    else if(stn==18) stn_id = 2;
    else {
      // cerr<<"Station "<<stn<<" not recognised";
      cout<<"WARNING: Station "<<stn<<" not recognised"<<endl;
      continue;
    }

    if(momCuts && p > pLo && p < pHi) {

      p_tracks_[stn_id]->Fill(p);
      thetaY_tracks_[stn_id]->Fill(theta_y);
      Y_tracks_[stn_id]->Fill(y);
      R_tracks_[stn_id]->Fill(r);
      Phi_tracks_[stn_id]->Fill(phi);
      thetaY_vs_Y_tracks_[stn_id]->Fill(y, theta_y);
      thetaY_vs_R_tracks_[stn_id]->Fill(r, theta_y);
      thetaY_vs_Phi_tracks_[stn_id]->Fill(phi, theta_y);
      thetaY_vs_p_tracks_[stn_id]->Fill(p, theta_y);

    } else if(!momCuts) {

      p_tracks_[stn_id]->Fill(p);
      thetaY_tracks_[stn_id]->Fill(theta_y);
      Y_tracks_[stn_id]->Fill(y);
      R_tracks_[stn_id]->Fill(r);
      Phi_tracks_[stn_id]->Fill(phi);
      thetaY_vs_Y_tracks_[stn_id]->Fill(y, theta_y);
      thetaY_vs_R_tracks_[stn_id]->Fill(r, theta_y);
      thetaY_vs_Phi_tracks_[stn_id]->Fill(phi, theta_y);
      thetaY_vs_p_tracks_[stn_id]->Fill(p, theta_y);

    }

    // Slice momentum
    for ( int i_slice = 0; i_slice < p_nSlices; i_slice++ ) { 

      int lo = i_slice*p_step; 
      int hi = p_step + lo;

      if(p >= double(lo) && p < double(hi)) { 

        p_tracks_p_slices_[stn_id].at(i_slice)->Fill(p);
        thetaY_tracks_p_slices_[stn_id].at(i_slice)->Fill(theta_y);
        Y_tracks_p_slices_[stn_id].at(i_slice)->Fill(y);
        R_tracks_p_slices_[stn_id].at(i_slice)->Fill(r);
        Phi_tracks_p_slices_[stn_id].at(i_slice)->Fill(phi);
        thetaY_vs_Y_tracks_p_slices_[stn_id].at(i_slice)->Fill(y, theta_y);
        thetaY_vs_R_tracks_p_slices_[stn_id].at(i_slice)->Fill(r, theta_y);
        thetaY_vs_Phi_tracks_p_slices_[stn_id].at(i_slice)->Fill(phi, theta_y);
        thetaY_vs_p_tracks_p_slices_[stn_id].at(i_slice)->Fill(p, theta_y);

      }

    }

    // Slice y
    for ( int i_slice = 0; i_slice < y_nSlices; i_slice++ ) { 

      int lo = -45 + i_slice*y_step; 
      int hi = y_step + lo;

      if(y >= double(lo) && y < double(hi)) { 

        if(momCuts && p > pLo && p < pHi) {

          p_tracks_y_slices_[stn_id].at(i_slice)->Fill(p);
          thetaY_tracks_y_slices_[stn_id].at(i_slice)->Fill(theta_y);
          Y_tracks_y_slices_[stn_id].at(i_slice)->Fill(y);
          thetaY_vs_Y_tracks_y_slices_[stn_id].at(i_slice)->Fill(y, theta_y);
          thetaY_vs_R_tracks_y_slices_[stn_id].at(i_slice)->Fill(r, theta_y);
          thetaY_vs_Phi_tracks_y_slices_[stn_id].at(i_slice)->Fill(phi, theta_y);
          thetaY_vs_p_tracks_y_slices_[stn_id].at(i_slice)->Fill(p, theta_y);

        } else if(!momCuts) {

          p_tracks_y_slices_[stn_id].at(i_slice)->Fill(p);
          thetaY_tracks_y_slices_[stn_id].at(i_slice)->Fill(theta_y);
          Y_tracks_y_slices_[stn_id].at(i_slice)->Fill(y);
          thetaY_vs_Y_tracks_y_slices_[stn_id].at(i_slice)->Fill(y, theta_y);
          thetaY_vs_R_tracks_y_slices_[stn_id].at(i_slice)->Fill(r, theta_y);
          thetaY_vs_Phi_tracks_y_slices_[stn_id].at(i_slice)->Fill(phi, theta_y);
          thetaY_vs_p_tracks_y_slices_[stn_id].at(i_slice)->Fill(p, theta_y);

        }

      }

    }

  }

  // Write to output

  // Momentum interval string
  string momSliceString = to_string(int(pLo))+"_"+to_string(int(pHi))+"_MeV";

  // Set output directory
  output->mkdir(momSliceString.c_str()); output->mkdir((momSliceString+"/AllDecays").c_str()); 
  output->cd((momSliceString+"/AllDecays").c_str());

  cout<<"-----> Writing decay histograms"<<endl;

  p_decays->Write();
  Y_decays->Write();
  R_decays->Write();
  Phi_decays->Write();
  thetaY_decays->Write();
  thetaY_vs_Y_decays->Write();
  thetaY_vs_R_decays->Write();
  thetaY_vs_Phi_decays->Write();
  thetaY_vs_p_decays->Write();

  output->mkdir((momSliceString+"/AllDecays/MomBins").c_str()); 
  output->cd((momSliceString+"/AllDecays/MomBins").c_str());

  for ( int i_slice(0); i_slice < p_nSlices; i_slice++ ) {

    p_decays_p_slices_.at(i_slice)->Write(); 
    Y_decays_p_slices_.at(i_slice)->Write(); 
    R_decays_p_slices_.at(i_slice)->Write(); 
    Phi_decays_p_slices_.at(i_slice)->Write(); 
    thetaY_decays_p_slices_.at(i_slice)->Write();
    thetaY_vs_Y_decays_p_slices_.at(i_slice)->Write();
    thetaY_vs_R_decays_p_slices_.at(i_slice)->Write();
    thetaY_vs_Phi_decays_p_slices_.at(i_slice)->Write();
    thetaY_vs_p_decays_p_slices_.at(i_slice)->Write();  

  }

  output->mkdir((momSliceString+"/AllDecays/VertPosBins").c_str());
  output->cd((momSliceString+"/AllDecays/VertPosBins").c_str());

  for ( int i_slice(0); i_slice < y_nSlices; i_slice++ ) {

    p_decays_y_slices_.at(i_slice)->Write(); 
    Y_decays_y_slices_.at(i_slice)->Write(); 
    R_decays_y_slices_.at(i_slice)->Write(); 
    Phi_decays_y_slices_.at(i_slice)->Write(); 
    thetaY_decays_y_slices_.at(i_slice)->Write();
    thetaY_vs_Y_decays_y_slices_.at(i_slice)->Write();
    thetaY_vs_R_decays_y_slices_.at(i_slice)->Write();
    thetaY_vs_Phi_decays_y_slices_.at(i_slice)->Write();
    thetaY_vs_p_decays_y_slices_.at(i_slice)->Write();

  }

  // Combine stations S12&S18
  cout<<"-----> Combining stations"<<endl;

  p_tracks_[0]->Add(p_tracks_[1], p_tracks_[2]);
  Y_tracks_[0]->Add(Y_tracks_[1], Y_tracks_[2]);
  R_tracks_[0]->Add(R_tracks_[1], R_tracks_[2]);
  Phi_tracks_[0]->Add(Phi_tracks_[1], Phi_tracks_[2]);
  thetaY_tracks_[0]->Add(thetaY_tracks_[1], thetaY_tracks_[2]);
  thetaY_vs_Y_tracks_[0]->Add(thetaY_vs_Y_tracks_[1], thetaY_vs_Y_tracks_[2]);
  thetaY_vs_R_tracks_[0]->Add(thetaY_vs_R_tracks_[1], thetaY_vs_R_tracks_[2]);
  thetaY_vs_Phi_tracks_[0]->Add(thetaY_vs_Phi_tracks_[1], thetaY_vs_Phi_tracks_[2]);
  thetaY_vs_p_tracks_[0]->Add(thetaY_vs_p_tracks_[1], thetaY_vs_p_tracks_[2]);

  for ( int i_slice(0); i_slice < p_nSlices; i_slice++ ) {
    
    p_tracks_p_slices_[0].at(i_slice)->Add(p_tracks_p_slices_[1].at(i_slice), p_tracks_p_slices_[2].at(i_slice));
    Y_tracks_p_slices_[0].at(i_slice)->Add(Y_tracks_p_slices_[1].at(i_slice), Y_tracks_p_slices_[2].at(i_slice));
    R_tracks_p_slices_[0].at(i_slice)->Add(R_tracks_p_slices_[1].at(i_slice), R_tracks_p_slices_[2].at(i_slice));
    Phi_tracks_p_slices_[0].at(i_slice)->Add(Phi_tracks_p_slices_[1].at(i_slice), Phi_tracks_p_slices_[2].at(i_slice));
    thetaY_tracks_p_slices_[0].at(i_slice)->Add(thetaY_tracks_p_slices_[1].at(i_slice), thetaY_tracks_p_slices_[2].at(i_slice));
    thetaY_vs_Y_tracks_p_slices_[0].at(i_slice)->Add(thetaY_vs_Y_tracks_p_slices_[1].at(i_slice), thetaY_vs_Y_tracks_p_slices_[2].at(i_slice));
    thetaY_vs_p_tracks_p_slices_[0].at(i_slice)->Add(thetaY_vs_p_tracks_p_slices_[1].at(i_slice), thetaY_vs_p_tracks_p_slices_[2].at(i_slice));
  }

  for ( int i_slice(0); i_slice < y_nSlices; i_slice++ ) {
    
    p_tracks_y_slices_[0].at(i_slice)->Add(p_tracks_y_slices_[1].at(i_slice), p_tracks_y_slices_[2].at(i_slice));
    Y_tracks_y_slices_[0].at(i_slice)->Add(Y_tracks_y_slices_[1].at(i_slice), Y_tracks_y_slices_[2].at(i_slice));
    R_tracks_y_slices_[0].at(i_slice)->Add(R_tracks_y_slices_[1].at(i_slice), R_tracks_y_slices_[2].at(i_slice));
    Phi_tracks_y_slices_[0].at(i_slice)->Add(Phi_tracks_y_slices_[1].at(i_slice), Phi_tracks_y_slices_[2].at(i_slice));
    thetaY_tracks_y_slices_[0].at(i_slice)->Add(thetaY_tracks_y_slices_[1].at(i_slice), thetaY_tracks_y_slices_[2].at(i_slice));
    thetaY_vs_Y_tracks_y_slices_[0].at(i_slice)->Add(thetaY_vs_Y_tracks_y_slices_[1].at(i_slice), thetaY_vs_Y_tracks_y_slices_[2].at(i_slice));
    thetaY_vs_R_tracks_y_slices_[0].at(i_slice)->Add(thetaY_vs_R_tracks_y_slices_[1].at(i_slice), thetaY_vs_R_tracks_y_slices_[2].at(i_slice));
    thetaY_vs_Phi_tracks_y_slices_[0].at(i_slice)->Add(thetaY_vs_Phi_tracks_y_slices_[1].at(i_slice), thetaY_vs_Phi_tracks_y_slices_[2].at(i_slice));
    thetaY_vs_p_tracks_y_slices_[0].at(i_slice)->Add(thetaY_vs_p_tracks_y_slices_[1].at(i_slice), thetaY_vs_p_tracks_y_slices_[2].at(i_slice));
     
  }

  cout<<"-----> Writing track histograms"<<endl;

  output->mkdir((momSliceString+"/Tracks").c_str()); 
  output->mkdir((momSliceString+"/Tracks/MomBins").c_str());
  output->mkdir((momSliceString+"/Tracks/VertPosBins").c_str()); 

  for (int i_stn = 0; i_stn < n_stn; i_stn++) {

    output->cd((momSliceString+"/Tracks").c_str());

    p_tracks_[i_stn]->Write();
    Y_tracks_[i_stn]->Write();
    R_tracks_[i_stn]->Write();
    Phi_tracks_[i_stn]->Write();
    thetaY_tracks_[i_stn]->Write();
    thetaY_vs_Y_tracks_[i_stn]->Write();
    thetaY_vs_R_tracks_[i_stn]->Write();
    thetaY_vs_Phi_tracks_[i_stn]->Write();
    thetaY_vs_p_tracks_[i_stn]->Write();

    output->cd((momSliceString+"/Tracks/MomBins").c_str());

    for ( int i_slice(0); i_slice < p_nSlices; i_slice++ ) {
      
      p_tracks_p_slices_[i_stn].at(i_slice)->Write();
      Y_tracks_p_slices_[i_stn].at(i_slice)->Write();
      R_tracks_p_slices_[i_stn].at(i_slice)->Write();
      Phi_tracks_p_slices_[i_stn].at(i_slice)->Write();
      thetaY_tracks_p_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_Y_tracks_p_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_R_tracks_p_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_Phi_tracks_p_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_p_tracks_p_slices_[i_stn].at(i_slice)->Write();
           
    }

    output->cd((momSliceString+"/Tracks/VertPosBins").c_str());

    for ( int i_slice(0); i_slice < y_nSlices; i_slice++ ) {

      p_tracks_y_slices_[i_stn].at(i_slice)->Write();
      Y_tracks_y_slices_[i_stn].at(i_slice)->Write();
      R_tracks_y_slices_[i_stn].at(i_slice)->Write();
      Phi_tracks_y_slices_[i_stn].at(i_slice)->Write();
      thetaY_tracks_y_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_Y_tracks_y_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_R_tracks_y_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_Phi_tracks_y_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_p_tracks_y_slices_[i_stn].at(i_slice)->Write();
           
    }

  }

  cout<<"-----> Written plots."<<endl;

  return;

}

/*

  After compilation, run either like: 

  ./Plotter.out /pnfs/GM2/persistent/EDM/MC/TrackerAcceptance/jprice_mc_digit_trackAcceptance/trackerAcceptanceTrees.1567206067.root output/plots.root truth

  or over the full dataset 

  . submit.sh 3 truth

*/

int main(int argc, char *argv[]) {

  // Config variables
  bool quality = true;
  bool momCuts = true;
  bool timeCuts = true;
  bool truth = true;

  // Arguments
  string inFileName = argv[1]; 
  string outFileName = argv[2]; 
  string truthStr = argv[3];

  if(truthStr!="truth") truth=false;

  string treeNameAllDecays = "trackerAcceptanceTrees/AllDecays";
  string treeNameTracks = "trackerAcceptanceTrees/TrackReco";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName .c_str());

  // Get trees
  TTree *treeAllDecays = (TTree*)fin->Get(treeNameAllDecays.c_str());
  TTree *treeTracks = (TTree*)fin->Get(treeNameTracks.c_str());

  cout<<"\nOpened tree:\t"<<treeNameAllDecays<<" "<<treeAllDecays<<" from file "<<inFileName<<" "<<fin<<endl;
  cout<<"\nOpened tree:\t"<<treeNameTracks<<" "<<treeTracks<<" from file "<<inFileName<<" "<<fin<<endl;

  // Book output
  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");

  // Run once over the full momentum range and once over the EDM analysis momentum range 

  cout<<"-----> Running momentum slice 0 < p [MeV] < "<<pmax<<endl;

  Run(treeAllDecays, treeTracks, fout, quality, momCuts, timeCuts, truth, 0, pmax);

  cout<<"-----> Running momentum slice 1000 < p [MeV] < 2500"<<endl;

  Run(treeAllDecays, treeTracks, fout, quality, momCuts, timeCuts, truth, 1000, 2500);

  // Close
  fout->Close();
  fin->Close();

  cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;
   
  return 0;

}