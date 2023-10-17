/*

Samuel Grant

Primary plotting code for producing EDM base histograms. 

This is the first step in the data-based section of the EDM analysis chain. 

*/

#include "../Common/Utils.h"
#include "Plotter.h"

double pLo = 1000; 
double pHi = 2500;

double f_c = 6.71; // MHz
// double T_y = 1/(f_c*sqrt(n)); 
// double T_x = 1/(f_c - (f_c*sqrt(1-n)));

double RandomisedTimeVBO(TRandom3 *rand, double time, string ds) {

  double n;
  if(ds=="Run-1a" || ds=="Run-1d") n = 0.108; // Run-1a / Run-1d
  else if(ds=="Run-1b" || ds=="Run-1c") n = 0.120; // Run-1b / Run-1c
  else cerr<<"RandomisedTimeVB(): Dataset not valid!";

  double T_y = 1/(f_c*sqrt(n));

  return rand->Uniform(time-T_y/2, time+T_y/2);

}

void Run(TTree *tree, TFile *output, string dataset="Run-1a", bool quality=true, bool timeCuts=true, bool momCuts=true, bool randCorr=true, bool vertCorr=false) {

  // Tracker station array
  string stns[] = {"S12S18", "S12", "S18"}; 
  int n_stn = sizeof(stns)/sizeof(stns[0]);

  // Counting 
  TH1D *vertices_[n_stn];
  TH1D *quality_vertices_[n_stn];

  // Simultaneous analysis
  TH1D *momentum_noCuts_[n_stn];
  TH1D *momentum_[n_stn];
  TH1D *momY_[n_stn];
  TH2D *decayX_vs_decayZ_[n_stn];
  TH1D *wiggle_[n_stn];
  TH1D *wiggle_mod_A_[n_stn]; // Using T-method momentum cut, to obtain phase
  TH1D *wiggle_mod_B_[n_stn]; // Using EDM momentum cuts, for theta_y fit demoninator 
  TH1D *thetaY_[n_stn];
  TH2D *thetaY_vs_time_[n_stn];
  TH2D *thetaY_vs_time_mod_[n_stn];
  TH2D *thetaY_vs_momentum_noCuts_[n_stn];
  TH2D *thetaY_vs_momentum_vertexCuts_[n_stn];
  TH2D *thetaY_vs_momentum_allCuts_[n_stn];

  // Momentum binned analysis
  vector<TH1D*> thetaY_mom_slices_[n_stn];
  vector<TH1D*> Y_mom_slices_[n_stn];
  vector<TH1D*> pY_mom_slices_[n_stn];
  vector<TH1D*> p_mom_slices_[n_stn];
  vector<TH1D*> wiggle_mod_mom_slices_[n_stn];
  vector<TH2D*> thetaY_vs_time_mod_slices_[n_stn];

  // Slice momentum
  int step = 250;
  int nSlices = pmax/step;

  // Vertical offset correction 
  // Needs to be updated for run-1_paper branch
  TFile *vertCorrFile = TFile::Open(("correctionPlots/verticalOffsetFits_"+dataset+"_"+to_string(step)+"MeV_BQ.root").c_str()); 
  vector<vector<TF1*>> verticalOffsetFits_;

  // This is just for consitency with sim when we make theta_y reweighting ratios, it's fairly abitrary 
  double boostFactor = 5e3*(1/gmagic);

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    vector<TF1*> verticalOffsetFitsPerStn_;

    vertices_[i_stn] = new TH1D((stns[i_stn]+"_Vertices").c_str(), ";;Decay vertices", 1, 0, 1); 
    quality_vertices_[i_stn] = new TH1D((stns[i_stn]+"_Quality_Vertices").c_str(), ";;Quality Decay vertices", 1, 0, 1); 
    momentum_[i_stn] = new TH1D((stns[i_stn]+"_Momentum").c_str(), ";Decay vertex  momentum [MeV];Decay vertices", int(pmax), 0, pmax); 
    momentum_noCuts_[i_stn] = new TH1D((stns[i_stn]+"_Momentum_NoCuts").c_str(), ";Decay vertex momentum [MeV];Decay vertices", int(pmax), 0, pmax); 
    momY_[i_stn] = new TH1D((stns[i_stn]+"_MomentumY").c_str(), ";Decay vertex momentum Y [MeV];Decay vertices", 1000, -60, 60); 
    decayX_vs_decayZ_[i_stn] = new TH2D((stns[i_stn]+"_DecayX_vs_DecayZ").c_str(), ";Decay vertex position Z [mm];Decay vertex position X [mm]", 800, -8000, 8000, 800, -8000, 8000);
    wiggle_[i_stn] = new TH1D((stns[i_stn]+"_Wiggle").c_str(), ";Decay vertex time [#mus];Decay vertices", 2700, 0, 2700*T_c);
    wiggle_mod_A_[i_stn] = new TH1D((stns[i_stn]+"_Wiggle_Modulo_A").c_str(), ";t_{g#minus2}^{mod} [#mus];Decay vertices / 149.2 ns", 29, 0, g2Period); 
    wiggle_mod_B_[i_stn] = new TH1D((stns[i_stn]+"_Wiggle_Modulo_B").c_str(), ";t_{g#minus2}^{mod} [#mus];Decay vertices / 149.2 ns", 29, 0, g2Period);
    thetaY_[i_stn] = new TH1D((stns[i_stn]+"_ThetaY").c_str(), ";#theta_{y} [mrad];Decay vertices", 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_time_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time").c_str(), ";Decay time [#mus]; #theta_{y} [mrad] / 149.2 ns ", 2700, 0, 2700*T_c, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_time_mod_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo").c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 149.2 ns", 29, 0, g2Period, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
    thetaY_vs_momentum_noCuts_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Momentum_NoCuts").c_str(), ";Decay vertex momentum [MeV]; #theta_{y} [mrad] / 10 MeV ", 312, 0, 3120, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
    thetaY_vs_momentum_vertexCuts_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Momentum_VertexCuts").c_str(), ";Decay vertex momentum [MeV]; #theta_{y} [mrad] / 10 MeV ", 312, 0, 3120, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);
    thetaY_vs_momentum_allCuts_[i_stn] = new TH2D((stns[i_stn]+"_ThetaY_vs_Momentum_AllCuts").c_str(), ";Decay vertex momentum [MeV]; #theta_{y} [mrad] / 10 MeV ", 312, 0, 3120, 1000, -TMath::Pi()*boostFactor, TMath::Pi()*boostFactor);

    // Slice momentum
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = 0 + i_slice*step; 
      int hi = step + i_slice*step;

      // Store vertical offset correction fit
      TGraphErrors* gr_verticalOffsetFit = (TGraphErrors*)vertCorrFile->Get(("MomBinnedAna/"+stns[i_stn]+"_ThetaY_vs_Time_Fit_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str());
      if(gr_verticalOffsetFit==0) verticalOffsetFitsPerStn_.push_back(0);
      else verticalOffsetFitsPerStn_.push_back(gr_verticalOffsetFit->GetFunction("DoubleExponentialFunc"));

      TH1D *h_p_mom_slices = new TH1D((stns[i_stn]+"_Momentum_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay vertex momentum [MeV];Decay vertices",  int(pmax), 0, pmax);
      p_mom_slices_[i_stn].push_back(h_p_mom_slices);

      TH1D *h_pY_mom_slices = new TH1D((stns[i_stn]+"_MomentumY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Decay vertex momentum Y [MeV];Decay vertices",  1000, -60, 60);
      pY_mom_slices_[i_stn].push_back(h_pY_mom_slices);

      TH1D *h_Y_mom_slice = new TH1D((stns[i_stn]+"_Y_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Vertical decay position [mm];Decay vertices",  180, -60, 60);
      Y_mom_slices_[i_stn].push_back(h_Y_mom_slice);

      TH1D *h_thetaY_mom_slice = new TH1D((stns[i_stn]+"_ThetaY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";#theta_{y} [mrad];Decay vertices",  500, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
      thetaY_mom_slices_[i_stn].push_back(h_thetaY_mom_slice);

      TH1D *h_wiggle_mod_mom_slice = new TH1D((stns[i_stn]+"_Wiggle_Modulo_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus];Decay vertices / 149.2 ns", 29, 0, g2Period); // in ana mom range
      wiggle_mod_mom_slices_[i_stn].push_back(h_wiggle_mod_mom_slice);

      TH2D *h_thetaY_vs_time_mod_slice = new TH2D((stns[i_stn]+"_ThetaY_vs_Time_Modulo_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 149.2 ns", 29, 0, g2Period, 1000, -TMath::Pi()*gmagic, TMath::Pi()*gmagic);
      thetaY_vs_time_mod_slices_[i_stn].push_back(h_thetaY_vs_time_mod_slice);

    }

    verticalOffsetFits_.push_back(verticalOffsetFitsPerStn_);

  }

  // Get branches 
  InitBranches br(tree);

  double targetPerc = 0;
  int64_t nEntries = tree->GetEntries();

  double muAngleMax = 0;
  
  // For FR randomisation
  TRandom3 *rand = new TRandom3(12345);

  cout<<"---> Looping through entries"<<endl;
   
  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    double time = br.decayTime * 1e-3; // ns -> us
    double x = br.decayVertexPosX; double y = br.decayVertexPosY; double z = br.decayVertexPosZ; 
    double px = br.decayVertexMomX; double py = -br.decayVertexMomY; double pz = br.decayVertexMomZ; 

    if(randCorr) {
      time = RandomisedTime(rand, time); // randomise out the FR
      time = RandomisedTimeVBO(rand, time, dataset); // randomise out the VBO
    }

    // Modulated g-2 time
    double g2ModTime = ModTime(time);

    int stn = br.station; 

    // Positron world momentum
    TVector3 eMom(px, py, pz); 
    TVector3 ePos(x, y, z);

    // Define the vertical angle
    double p = eMom.Mag(); 
    double theta_y = asin(py/p);
    theta_y = theta_y * 1e3; // rad -> mrad

    // Quality variable
    bool vertexQual = br.passDecayVertexQuality;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Station IDs for wrangling arrays
    int stn_id = -1;
    if(stn==12) stn_id = 1;
    else if(stn==18) stn_id = 2;
    else cerr<<"Station "<<stn<<" not recognised";

    // Useful plots before cuts
    vertices_[stn_id]->Fill(1.0);
    momentum_noCuts_[stn_id]->Fill(p);
    thetaY_vs_momentum_noCuts_[stn_id]->Fill(p, theta_y);

    // Cuts 
    if (timeCuts) { // Early time cut
      if (time < 7*g2Period) continue; 
      if (dataset=="Run-1d" && time < 12*g2Period) continue; // Late start for Run-1d
    }
    if (quality) { // Quality cuts
      // Includes hit vol and pval cuts 
      // Cuts defined here: https://cdcvs.fnal.gov/redmine/projects/gm2tracker/repository/revisions/develop/entry/quality/TrackQuality_service.cc
      // Also see: https://gm2-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=16444
      if (!vertexQual) continue;
    }

    quality_vertices_[stn_id]->Fill(1.0);

    // Vertical offset correction (time only)
    if(vertCorr) {

      // Vertical offset correction in momentum and time
      TF1 *verticalOffsetFunc; 

      // Slice momentum 
      for ( int i_slice = 0; i_slice < nSlices; i_slice++) { 

        int lo = 0 + i_slice*step; 
        int hi = step + i_slice*step;

        if(p >= double(lo) && p < double(hi)) {

          verticalOffsetFunc = (verticalOffsetFits_.at(stn_id)).at(i_slice);

        }

      }

      theta_y = theta_y - verticalOffsetFunc->Eval(time);

    } 

    thetaY_vs_momentum_vertexCuts_[stn_id]->Fill(p, theta_y);
    decayX_vs_decayZ_[stn_id]->Fill(z, x);
 
    // T-method energy cut. See Fienberg thesis figure 2.10. 
    if(p > 1700  && p < pmax) {
      wiggle_[stn_id]->Fill(time);
      wiggle_mod_A_[stn_id]->Fill(g2ModTime);
    } 

    // EDM cuts
    if( momCuts && p > pLo && p < pHi) { 

      momentum_[stn_id]->Fill(p);
      momY_[stn_id]->Fill(py);
      thetaY_[stn_id]->Fill(theta_y);
      thetaY_vs_time_[stn_id]->Fill(time, theta_y);
      thetaY_vs_time_mod_[stn_id]->Fill(g2ModTime, theta_y);
      thetaY_vs_momentum_allCuts_[stn_id]->Fill(p, theta_y);
      wiggle_mod_B_[stn_id]->Fill(g2ModTime);

    } else if(!momCuts) { 

      momentum_[stn_id]->Fill(p);
      momY_[stn_id]->Fill(py);
      thetaY_[stn_id]->Fill(theta_y);
      thetaY_vs_time_[stn_id]->Fill(time, theta_y);
      thetaY_vs_time_mod_[stn_id]->Fill(g2ModTime, theta_y);
      thetaY_vs_momentum_allCuts_[stn_id]->Fill(p, theta_y);
      wiggle_mod_B_[stn_id]->Fill(g2ModTime);

    }


    // Slice momentum 
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = 0 + i_slice*step; 
      int hi = step + i_slice*step;

      if(p >= double(lo) && p < double(hi)) { 

        wiggle_mod_mom_slices_[stn_id].at(i_slice)->Fill(g2ModTime);
        thetaY_vs_time_mod_slices_[stn_id].at(i_slice)->Fill(g2ModTime, theta_y);

        // Other scans 
        thetaY_mom_slices_[stn_id].at(i_slice)->Fill(theta_y);
        Y_mom_slices_[stn_id].at(i_slice)->Fill(y);
        pY_mom_slices_[stn_id].at(i_slice)->Fill(py);
        p_mom_slices_[stn_id].at(i_slice)->Fill(p);

      }

    }

    if(100*float(entry) / nEntries > targetPerc) {
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }


  }
  // Combine stations
  vertices_[0]->Add(vertices_[1], vertices_[2]);
  quality_vertices_[0]->Add(quality_vertices_[1], quality_vertices_[2]);
  momentum_[0]->Add(momentum_[1], momentum_[2]);
  momentum_noCuts_[0]->Add(momentum_noCuts_[1], momentum_noCuts_[2]);
  wiggle_[0]->Add(wiggle_[1], wiggle_[2]);
  wiggle_mod_A_[0]->Add(wiggle_mod_A_[1], wiggle_mod_A_[2]);
  wiggle_mod_B_[0]->Add(wiggle_mod_B_[1], wiggle_mod_B_[2]);
  thetaY_[0]->Add(thetaY_[1], thetaY_[2]);
  thetaY_vs_time_[0]->Add(thetaY_vs_time_[1], thetaY_vs_time_[2]);
  thetaY_vs_time_mod_[0]->Add(thetaY_vs_time_mod_[1], thetaY_vs_time_mod_[2]);
  decayX_vs_decayZ_[0]->Add(decayX_vs_decayZ_[1], decayX_vs_decayZ_[2]);
  momY_[0]->Add(momY_[1], momY_[2]);
  thetaY_vs_momentum_noCuts_[0]->Add(thetaY_vs_momentum_noCuts_[1], thetaY_vs_momentum_noCuts_[2]);
  thetaY_vs_momentum_vertexCuts_[0]->Add(thetaY_vs_momentum_vertexCuts_[1], thetaY_vs_momentum_vertexCuts_[2]);
  thetaY_vs_momentum_allCuts_[0]->Add(thetaY_vs_momentum_allCuts_[1], thetaY_vs_momentum_allCuts_[2]);

   for ( int i_slice(0); i_slice < nSlices; i_slice++ ) {

      wiggle_mod_mom_slices_[0].at(i_slice)->Add(wiggle_mod_mom_slices_[1].at(i_slice), wiggle_mod_mom_slices_[2].at(i_slice));
      thetaY_vs_time_mod_slices_[0].at(i_slice)->Add(thetaY_vs_time_mod_slices_[1].at(i_slice), thetaY_vs_time_mod_slices_[2].at(i_slice));
      thetaY_mom_slices_[0].at(i_slice)->Add(thetaY_mom_slices_[1].at(i_slice), thetaY_mom_slices_[2].at(i_slice));
      Y_mom_slices_[0].at(i_slice)->Add(Y_mom_slices_[1].at(i_slice), Y_mom_slices_[2].at(i_slice));
      pY_mom_slices_[0].at(i_slice)->Add(pY_mom_slices_[1].at(i_slice), pY_mom_slices_[2].at(i_slice));
      p_mom_slices_[0].at(i_slice)->Add(p_mom_slices_[1].at(i_slice), p_mom_slices_[2].at(i_slice));
      
  }

  // Write to output

  output->mkdir("SimultaneousAnalysis"); output->mkdir("MomentumBinnedAnalysis");  

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    output->cd("SimultaneousAnalysis");

    vertices_[i_stn]->Write();
    quality_vertices_[i_stn]->Write();
    momentum_[i_stn]->Write();
    momentum_noCuts_[i_stn]->Write();
    wiggle_[i_stn]->Write();
    wiggle_mod_A_[i_stn]->Write();
    wiggle_mod_B_[i_stn]->Write();
    thetaY_[i_stn]->Write();
    thetaY_vs_time_[i_stn]->Write();
    thetaY_vs_time_mod_[i_stn]->Write();
    decayX_vs_decayZ_[i_stn]->Write();
    momY_[i_stn]->Write();
    thetaY_vs_momentum_noCuts_[i_stn]->Write();
    thetaY_vs_momentum_vertexCuts_[i_stn]->Write();
    thetaY_vs_momentum_allCuts_[i_stn]->Write();

   for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) {

      output->cd("MomentumBinnedAnalysis"); 

      thetaY_mom_slices_[i_stn].at(i_slice)->Write();
      Y_mom_slices_[i_stn].at(i_slice)->Write();
      pY_mom_slices_[i_stn].at(i_slice)->Write();
      p_mom_slices_[i_stn].at(i_slice)->Write();
      wiggle_mod_mom_slices_[i_stn].at(i_slice)->Write();
      thetaY_vs_time_mod_slices_[i_stn].at(i_slice)->Write();
      
    }

  }

  cout<<"---> Written plots"<<endl;

  vertCorrFile->Close();

  return;

}

/*

Run test:

./Plotter.out /gm2/data/g2be/Production/Trees/Run1/trackRecoTrees_15921.root output/trackRecoPlots_15921_test.root Run-1a

Run full:

. submit.sh 3 Truth 

OR 

. submit.sh 3 Run-1a

(modify submit.sh accordingly)

*/

int main(int argc, char *argv[]) {

  bool quality = true;
  bool momCuts = true;
  bool timeCuts = true;
  bool randCorr = true;
  bool vertCorr = false;

  string inFileName = argv[1]; 
  string outFileName = argv[2]; 
  string dataset = argv[3];

  string treeName = "trackAndTrackCalo/tree";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName .c_str());
  
  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;

  // Book output
  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");
  
  // Fill histograms
  Run(tree, fout, dataset, quality, timeCuts, momCuts, randCorr, vertCorr);

  // Close
  fout->Close();
  fin->Close();

  cout<<"\nDone. Histogram written to:\t"<<outFileName<<" "<<fout<<endl;
  
  return 0;

}
