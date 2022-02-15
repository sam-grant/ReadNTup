#ifndef Plotter_h
#define Plotter_h

// ROOT includes
#include "TTree.h"

using namespace std;

class InitBranchesAllDecays { 

public: 

   unsigned int posiTrackID;
   float posiInitP;
   float posiInitPX;
   float posiInitPY;
   float posiInitPZ;
   float posiInitE;
   float posiInitTime;
   float posiInitPosX;
   float posiInitPosY;
   float posiInitPosZ;

   // Declare constructer
   InitBranchesAllDecays(TTree *tree);

};

// Constructer 1
InitBranchesAllDecays::InitBranchesAllDecays(TTree* tree) {

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

class InitBranchesTracks { 

public: 

   float decayVertexPosX;
   float decayVertexPosY;
   float decayVertexPosZ;
   float decayVertexMom;
   float decayVertexMomX;
   float decayVertexMomY;
   float decayVertexMomZ;
   float decayVertexTime;
   bool hitVolume;
   float pValue;
   int station;
   int nHits;
   bool passVertexQuality;

   // Declare constructer
   InitBranchesTracks(TTree *tree);

};

// Constructer 2
InitBranchesTracks::InitBranchesTracks(TTree* tree) {

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
   tree->SetBranchAddress("station", &station);
   tree->SetBranchAddress("nHits", &nHits);
   tree->SetBranchAddress("passVertexQuality", &passVertexQuality);

}

#endif