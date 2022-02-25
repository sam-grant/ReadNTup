#ifndef PlotterAllDecays_h
#define PlotterAllDecays_h

// ROOT includes
#include "TTree.h"

using namespace std;

class InitBranches { 

public: 

   int posiTrackID_;
   float posiInitP_;
   float posiInitPX_;
   float posiInitPY_;
   float posiInitPZ_;
   float posiInitE_;
   float posiInitTime_;
   float posiInitPosX_;
   float posiInitPosY_;
   float posiInitPosZ_;

   float decayVertexPosX_;
   float decayVertexPosY_;
   float decayVertexPosZ_;
   float decayVertexMom_;
   float decayVertexMomX_;
   float decayVertexMomY_;
   float decayVertexMomZ_;
   float decayVertexTime_;
   bool hitVolume_;
   float pValue_;
   int station_;
   int nHits_;
   bool passVertexQuality_;

   // Declare constructer
   InitBranchesAllDecays(TTree *tree);
   InitBranchesTracks(TTree *tree);

};

// Constructer 1
InitBranches::InitBranchesAllDecays(TTree* tree) {

   tree->SetBranchAddress("posiTrackID", &posiTrackID);
   tree->SetBranchAddress("posiInitP", &posiInitP);
   tree->SetBranchAddress("posiInitPX", &posiInitPX);
   tree->SetBranchAddress("posiInitPY", &posiInitPY);
   tree->SetBranchAddress("posiInitPZ", &posiInitPZ);
   tree->SetBranchAddress("posiInitE", &posiInitE);
   tree->SetBranchAddress("posiInitTime", &posiInitTime);
   tree->SetBranchAddress("posiInitPosX", &posiInitPosX);
   tree->SetBranchAddress("posiInitPosY", &posiInitPosY);
   tree->SetBranchAddress("posiInitPosZ", &posiInitPosZ);

}

// Constructer 2
InitBranches::InitBranchesTracks(TTree* tree) {

   tree->SetBranchAddress("decayVertexPosX", &decayVertexPosX);
   tree->SetBranchAddress("decayVertexPosY", &decayVertexPosY);
   tree->SetBranchAddress("decayVertexPosZ", &decayVertexPosZ);
   tree->SetBranchAddress("decayVertexMom", &decayVertexMom);
   tree->SetBranchAddress("decayVertexMomX", &decayVertexMomX);
   tree->SetBranchAddress("decayVertexMomY", &decayVertexMomY);
   tree->SetBranchAddress("decayVertexMomZ", &decayVertexMomZ);
   tree->SetBranchAddress("decayVertexTime", &decayVertexTime);
   tree->SetBranchAddress("hitVolume", &hitVolume);
   tree->SetBranchAddress("pValue", &pValue);
   tree->SetBranchAddress("nHits", &nHits);
   tree->SetBranchAddress("passVertexQuality", &passVertexQuality);

}

#endif