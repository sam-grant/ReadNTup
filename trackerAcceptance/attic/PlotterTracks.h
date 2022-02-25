#ifndef PlotterTracks_h
#define PlotterTracks_h

// ROOT includes
#include "TTree.h"

using namespace std;

class InitBranchesTracks { 

public: 


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
   InitBranchesTracks(TTree *tree);

};

// Constructer
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