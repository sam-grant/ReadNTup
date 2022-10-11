// Read ROOT trees 
// Count raw tracks, quality tracks, raw vertices, quality vertices
// Sam Grant
// For data

#include <iostream>
#include <vector>

#include "Plotter.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

double mMu = 105.6583715; // MeV
double aMu = 11659208.9e-10; 
double gmagic = std::sqrt( 1.+1./aMu );
double pmax = 1.01 * mMu * gmagic; // 3127.1144

TTree *InitTree(string fileName, string treeName) { 

   // ++++++++++++++ Open tree and load branches ++++++++++++++
   // Get file
   TFile *fin = TFile::Open(fileName.c_str());
   //cout<<"\nOpened tree:\t"<<fileName<<" "<<fin<<endl;

   // Get tree
   TTree *tree = (TTree*)fin->Get(treeName.c_str());

  //cout<<"\nOpened tree:"<<treeName<<" "<<tree<<" from file "<<fileName<<" "<<fin<<endl;

   return tree;

}

void Run(TTree *tree) {

  // Get branches (using header file)
  InitBranches br(tree);

  // Entries is number of tracks
  int64_t nEntries = tree->GetEntries();

  int64_t nSubRuns = 0;
  int64_t lastID = 0;
  int64_t nFills = 0;
  int64_t nLaserFills = 0;
  int64_t nTracks = 0;
  int64_t nQTracks = 0;
  int64_t nVertices = 0;
  int64_t nQVertices = 0;

  // Define outside loop, should be one run number per tree
  unsigned int runNum = 0;
  unsigned int lastSubRunNum = 0;

  for(int64_t entry = 0; entry < nEntries; entry++) {

    tree->GetEntry(entry);

    // Exclude laser fills, is this needed?
    //if(br.inFillLaserCount>0) continue;

    // Make sure there are stations being hit
    if(br.station!=12 && br.station!=18) continue;

    // Count fills
    runNum = br.runNum;
    unsigned int subRunNum = br.subRunNum;
    unsigned int eventNum = br.eventNum;
    
    // Define a unique ID for the fill (probably overkill)
    int64_t fillID = 1e8*br.eventNum + 1e5*br.subRunNum + br.runNum;

    // Count fills
    if(fillID != lastID) { // New fill
      lastID = fillID; 
      nFills++; // all fills
      if(br.inFillLaserCount>0) nLaserFills++; // laser fills

    } 

    // Count subruns 
    if(subRunNum != lastSubRunNum) { // New subrun
      lastSubRunNum = subRunNum; 
      nSubRuns++;
    } 

    // Count tracks and vertices
    nTracks++;
    if(br.passTrackQuality) nQTracks++;
    if(br.hasDecayVertex) nVertices++;
    if(br.passDecayVertexQuality) nQVertices++;

  
  }
  
  cout<<runNum<<", "<<nSubRuns<<", "<<nFills<<", "<<nLaserFills<<", "<<nTracks<<", "<<nQTracks<<", "<<nVertices<<", "<<nQVertices<<endl;

  return;

}

int main(int argc, char *argv[]) { 
  
  string inFileName = argv[1];
  
  string treeName = "trackAndTrackCalo/tree";

  // Open tree and load branches
  TFile *fin = TFile::Open(inFileName .c_str());
  // Get tree
  TTree *tree = (TTree*)fin->Get(treeName.c_str());

  //cout<<"\nOpened tree:\t"<<treeName<<" "<<tree<<" from file "<<inFileName<<" "<<fin<<endl;
   
  // Get fills
  Run(tree);

  // Close
  fin->Close();
   
  return 0;

}
