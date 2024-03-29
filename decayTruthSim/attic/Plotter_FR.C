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

using namespace std;

double omegaAMagic = 0.00143934; // from gm2geom consts / kHz 
double g2Period = (2*TMath::Pi()/omegaAMagic) * 1e-3; // 4.3653239 us
double mMu = 105.6583715; // MeV
double aMu = 11659208.9e-10; 
double gmagic = std::sqrt( 1.+1./aMu );
double pmax = 1.01 * mMu * gmagic;

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

void Run(TTree *tree, TFile *output, bool quality) {

   // Book histograms
   // Bin out FR
   double binWidth = 0.148936;
   TH1D *momentum = new TH1D("Momentum", ";Track momentum [MeV];Tracks", int(pmax), 0, pmax); 
   TH1D *wiggle = new TH1D("Wiggle", ";Decay time [#mus];Tracks", 2700, 0, 2700*0.148936);
   TH1D *wiggle_mod = new TH1D("Wiggle_Modulo", ";t_{g#minus2}^{mod} [#mus];Tracks / 50 ns", 87, 0, g2Period); // int(g2Period/0.148936), 0, g2Period); // 87 bins for 50 ns
   TH1D *thetaY = new TH1D("ThetaY", ";#theta_{y} [mrad];Tracks", 180, -60, 60);
   TH2D *thetaY_vs_time = new TH2D("ThetaY_vs_Time", ";Decay time [#mus]; #theta_{y} [mrad] / 149 ns ", 2700, 0, 2700*0.148936, 180, -60, 60);
   TH2D *thetaY_vs_time_mod = new TH2D("ThetaY_vs_Time_Modulo", ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 50 ns", 87, 0,g2Period, 180, -60, 60);
   
   // Momentum scans of mod 
   vector<TH2D*> mom_slices_;
   vector<TH2D*> mom_sym_cuts_;
   vector<TH2D*> mom_min_scan_;
   vector<TH2D*> mom_max_scan_;

   // Other scans 
   vector<TH1D*> thetaY_mom_slices_;
   vector<TH1D*> Y_mom_slices_;
   vector<TH1D*> pY_mom_slices_;
   //vector<TH2D*> thetaA_mom_slices_;

   // Slice momentum in 200 MeV up to 3100 MeV
   int step = 200;
   int nSlices = pmax / step;

   for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = 0 + i_slice*step; 
      int hi = step + i_slice*step;

      // Mod 
      TH2D *h_mom_slice = new TH2D(("ThetaY_vs_Time_Modulo_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad]", 87, 0, g2Period, 180, -60, 60);
      mom_slices_.push_back(h_mom_slice);

      TH1D *h_thetaY_mom_slice = new TH1D(("ThetaY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";#theta_{y} [mrad];Tracks", 180, -60, 60);
      thetaY_mom_slices_.push_back(h_thetaY_mom_slice);

      TH1D *h_Y_mom_slice = new TH1D(("Y_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Vertical decay position [mm];Tracks", 180, -60, 60);
      Y_mom_slices_.push_back(h_Y_mom_slice);

      TH1D *h_pY_mom_slices = new TH1D(("MomentumY_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";Track momentum Y MeV];Tracks", 140, -70, 70);
      pY_mom_slices_.push_back(h_pY_mom_slices);

   }

   // Symmetric cuts
   for ( int i_slice = 0; i_slice < nSlices/2; i_slice++ ) { 

      int lo = 400 + i_slice*step; 
      int hi = 3000 - i_slice*step;

      TH2D *h_tmp = new TH2D(("ThetaY_vs_Time_Modulo_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad]", 87, 0, g2Period, 180, -60, 60);
      mom_sym_cuts_.push_back(h_tmp);

   }

    // Pmin scan
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = 0 + i_slice*step;
      int hi = pmax;

      TH2D *h_tmp = new TH2D(("ThetaY_vs_Time_Modulo_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad]", 87, 0, g2Period, 180, -60, 60);
      mom_min_scan_.push_back(h_tmp);

    } 

    // pmax scan
    for ( int i_slice = 0; i_slice < (nSlices-(600/step)); i_slice++ ) { 

      int lo = 600;// 0 + i_slice*step;
      int hi = 3000 - i_slice*step;

      TH2D *h_tmp = new TH2D(("ThetaY_vs_Time_Modulo_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str(), ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad]", 87, 0, g2Period, 180, -60, 60);
      mom_max_scan_.push_back(h_tmp);

    } 
   // ==================================================
//   int slice_count[nSlices];// = 0;
//   int sym_count[nSlices/2];// = 0;
//   int pmin_count[nSlices];// = 0;
//
//   for (int i_slice = 0; i_slice < nSlices; i_slice++) { 
//      slice_count[i_slice] = 0;
//      pmin_count[i_slice] = 0;
//      if(i_slice < nSlices/2) sym_count[i_slice] = 0;
//   }

   //int maxBinEntries = 10204; // 500e3 total / 49 files;


   // Get branches (using header file)
   InitBranches br(tree);

   int64_t nEntries = tree->GetEntries();

   for(int64_t entry = 0; entry < nEntries; entry++) {

      tree->GetEntry(entry);

      // Get variables
      double time = br.posiInitTime * 1e-3; // us
      double px = br.posiInitPX; // MeV
      double py = br.posiInitPY; // MeV
      double pz = br.posiInitPZ; // MeV
      double p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));//br.posiInitP; // MeV
      double pr = sqrt(pow(px,2)+pow(pz,2));
      double theta_y = atan2(py,pr) * 1e3; // mrad
      double g2ModTime = ModTime(time);
      double y = br.posiInitPosY;

      // Time cuts
      if(quality && (time < 30 || time > 300)) continue; 
 
      // g-2 cuts. See Fienberg thesis figure 2.10
      if(quality && (p > 1900 && p < pmax)) {
         wiggle->Fill(time);
         wiggle_mod->Fill(g2ModTime);
      }

      // EDM cuts
      if(quality && (p > 700 && p < 2400)) { 

         momentum->Fill(p);
         thetaY->Fill(theta_y);
         thetaY_vs_time->Fill(time, theta_y);
         thetaY_vs_time_mod->Fill(g2ModTime, theta_y);

      }

      // Slice momentum 
      for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

         int lo = 0 + i_slice*step; 
         int hi = step + i_slice*step;

         if(p >= double(lo) && p < double(hi)) {//} && slice_count[i_slice]<maxBinEntries) {

            mom_slices_.at(i_slice)->Fill(g2ModTime, theta_y);

            // Other scans 
            thetaY_mom_slices_.at(i_slice)->Fill(theta_y);
            Y_mom_slices_.at(i_slice)->Fill(y);
            pY_mom_slices_.at(i_slice)->Fill(py);

            //slice_count[i_slice]++;  
         }

         lo = 400 + i_slice*step; 
         hi = 3000 - i_slice*step;

         if(p >= double(lo) && p < double(hi) && i_slice < nSlices/2) {//} && sym_count[i_slice]<maxBinEntries) {
            mom_sym_cuts_.at(i_slice)->Fill(g2ModTime, theta_y);
            //sym_count[i_slice]++;
         }
         

         lo = 0 + i_slice*step;
         hi = pmax;

         if(p >= double(lo) && p < double(hi)) {//&& pmin_count[i_slice]<maxBinEntries) {
            mom_min_scan_.at(i_slice)->Fill(g2ModTime, theta_y);
            //pmin_count[i_slice]++;
         }

         lo = 600;
         hi = 3000 - i_slice*step;

         if(p >= double(lo) && p < double(hi) && i_slice < (nSlices - (600/step))) {//&& pmin_count[i_slice]<maxBinEntries) {
            mom_max_scan_.at(i_slice)->Fill(g2ModTime, theta_y);
            //pmin_count[i_slice]++;
         }

      }
         

   }

   // Write to output
   // Set output directory
   output->mkdir("MainPlots"); output->cd("MainPlots");
   momentum->Write();
   wiggle->Write();
   wiggle_mod->Write();
   thetaY->Write();
   thetaY_vs_time->Write();
   thetaY_vs_time_mod->Write();

   output->mkdir("MomSlices"); output->mkdir("MomSymCuts"); output->mkdir("MomMinScan"); output->mkdir("MomMaxScan");

   for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) {
      output->cd("MomSlices"); 
      mom_slices_.at(i_slice)->Write();
      thetaY_mom_slices_.at(i_slice)->Write();
      Y_mom_slices_.at(i_slice)->Write();
      pY_mom_slices_.at(i_slice)->Write();
      output->cd("MomSymCuts"); if(i_slice < nSlices/2) mom_sym_cuts_.at(i_slice)->Write();
      output->cd("MomMinScan"); mom_min_scan_.at(i_slice)->Write();
      output->cd("MomMaxScan"); if(i_slice < (nSlices - (600/step))) mom_max_scan_.at(i_slice)->Write();
   }

   return;

}

int main(int argc, char *argv[]) {

   bool quality = true;

   string inFileName = argv[1]; // "/pnfs/GM2/persistent/EDM/MC/dMu/TruthNTup/truthTrees.000.root";//
   string outFileName = argv[2]; //"tmp.root";

   string treeName = "phaseAnalyzer/g2phase";
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
