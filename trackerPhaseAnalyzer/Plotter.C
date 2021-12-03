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
   //TH1D *momentum = new TH1D("Momentum", ";Track momentum [MeV];Tracks", int(pmax), 0, pmax); 
   //TH1D *wiggle = new TH1D("Wiggle", ";Decay time [#mus];Tracks", 2700, 0, 2700*0.148936);
   //TH1D *wiggle_mod = new TH1D("Wiggle_Modulo", ";t_{g#minus2}^{mod} [#mus];Tracks / 50 ns", 87, 0, g2Period); // int(g2Period/0.148936), 0, g2Period); // 87 bins for 50 ns
   //TH1D *thetaY = new TH1D("ThetaY", ";#theta_{y} [mrad];Tracks", 180, -60, 60);
   //TH2D *thetaY_vs_time = new TH2D("ThetaY_vs_Time", ";Decay time [#mus]; #theta_{y} [mrad] / 149 ns ", 2700, 0, 2700*0.148936, 180, -60, 60);
   //TH2D *thetaY_vs_time_mod = new TH2D("ThetaY_vs_Time_Modulo", ";t_{g#minus2}^{mod} [#mus]; #theta_{y} [mrad] / 50 ns", 87, 0,g2Period, 180, -60, 60);
   TH2D *decayX_vs_decayZ = new TH2D("decayX_vs_decayZ", ";Decay vertex position X [mm];Decay vertex position Z [mm]", 800, -8000, 8000, 800, -8000, 8000);
   
   // Get branches (using header file)
   InitBranches br(tree);

   int64_t nEntries = tree->GetEntries();

   for(int64_t entry = 0; entry < nEntries; entry++) {

      tree->GetEntry(entry);

      // Get variables
      //double time = br.posiInitTime * 1e-3; // us
      //double px = br.posiInitPX; // MeV
      //double py = br.posiInitPY; // MeV
      //double pz = br.posiInitPZ; // MeV
      //double p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));//br.posiInitP; // MeV
      //double theta_y = atan2(py,p) * 1e3; // mrad
      //double g2ModTime = ModTime(time);
      double x = br.posiInitPosX;
      //double y = br.posiInitPosY;
      double z = br.posiInitPosZ;

      decayX_vs_decayZ->Fill(x,z);

      // Time cuts
/*      if(quality) { 

         if(time < 30 || time > 300) continue;
 
         // g-2 cuts. See Fienberg thesis figure 2.10
         if(p > 1900 && p < pmax) {
            wiggle->Fill(time);
            wiggle_mod->Fill(g2ModTime);
         }

         // EDM cuts
         if(p > 700 && p < 2400) { 

         momentum->Fill(p);
         thetaY->Fill(theta_y);
         thetaY_vs_time->Fill(time, theta_y);
         thetaY_vs_time_mod->Fill(g2ModTime, theta_y);
         
         }
         
      } else { 


      }*/

   }

   // Write to output
   // Set output directory
   output->mkdir("MainPlots"); output->cd("MainPlots");
   //momentum->Write();
   //wiggle->Write();
   //wiggle_mod->Write();
   //thetaY->Write();
   //thetaY_vs_time->Write();
   //thetaY_vs_time_mod->Write();
   decayX_vs_decayZ->Write();

   //output->mkdir("MomSlices"); output->mkdir("MomSymCuts"); output->mkdir("MomMinScan"); output->mkdir("MomMaxScan");

   return;

}

int main(int argc, char *argv[]) {

   bool quality = true;

   string inFileName = argv[1]; // "phaseAnalyzerTree.root";//argv[1]; // "/pnfs/GM2/persistent/EDM/MC/dMu/TruthNTup/truthTrees.000.root";//
   string outFileName = argv[2];//"phaseAnalyzerPlots.root";//argv[2]; //"tmp.root";

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
