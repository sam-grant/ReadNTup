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

   float recoVertexPosX;
   float recoVertexPosY;
   float recoVertexPosZ;
   float recoVertexMom;
   float recoVertexMomX;
   float recoVertexMomY;
   float recoVertexMomZ;
   float recoVertexTime;
   float trueVertexPosX;
   float trueVertexPosY;
   float trueVertexPosZ;
   float trueVertexMom;
   float trueVertexMomX;
   float trueVertexMomY;
   float trueVertexMomZ;
   float trueVertexTime;
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

   tree->SetBranchAddress("recoVertexPosX", &recoVertexPosX);
   tree->SetBranchAddress("recoVertexPosY", &recoVertexPosY);
   tree->SetBranchAddress("recoVertexPosZ", &recoVertexPosZ);
   tree->SetBranchAddress("recoVertexMom",  &recoVertexMom);
   tree->SetBranchAddress("recoVertexMomX", &recoVertexMomX);
   tree->SetBranchAddress("recoVertexMomY", &recoVertexMomY);
   tree->SetBranchAddress("recoVertexMomZ", &recoVertexMomZ);
   tree->SetBranchAddress("recoVertexTime", &recoVertexTime);
   tree->SetBranchAddress("trueVertexPosX", &trueVertexPosX);
   tree->SetBranchAddress("trueVertexPosY", &trueVertexPosY);
   tree->SetBranchAddress("trueVertexPosZ", &trueVertexPosZ);
   tree->SetBranchAddress("trueVertexMom",  &trueVertexMom);
   tree->SetBranchAddress("trueVertexMomX", &trueVertexMomX);
   tree->SetBranchAddress("trueVertexMomY", &trueVertexMomY);
   tree->SetBranchAddress("trueVertexMomZ", &trueVertexMomZ);
   tree->SetBranchAddress("trueVertexTime", &trueVertexTime);
   tree->SetBranchAddress("hitVolume", &hitVolume);
   tree->SetBranchAddress("pValue", &pValue);
   tree->SetBranchAddress("station", &station);
   tree->SetBranchAddress("nHits", &nHits);
   tree->SetBranchAddress("passVertexQuality", &passVertexQuality);

}

class InitBranchesTracks2 { 

public: 

   float recoVertexPosX;
   float recoVertexPosY;
   float recoVertexPosZ;
   float recoVertexMomX;
   float recoVertexMomY;
   float recoVertexMomZ;
   float recoTime;
   float trueVertexPosX;
   float trueVertexPosY;
   float trueVertexPosZ;
   float trueVertexMomX;
   float trueVertexMomY;
   float trueVertexMomZ;
   float trueTime;
   bool hitVolume;
   float pValue;
   int station;
   bool passVertexQuality;

   // Declare constructer
   InitBranchesTracks2(TTree *tree);

};

// Constructer 2
InitBranchesTracks2::InitBranchesTracks2(TTree* tree) {

   tree->SetBranchAddress("recoVertexPosX", &recoVertexPosX);
   tree->SetBranchAddress("recoVertexPosY", &recoVertexPosY);
   tree->SetBranchAddress("recoVertexPosZ", &recoVertexPosZ);
   tree->SetBranchAddress("recoVertexMomX", &recoVertexMomX);
   tree->SetBranchAddress("recoVertexMomY", &recoVertexMomY);
   tree->SetBranchAddress("recoVertexMomZ", &recoVertexMomZ);
   tree->SetBranchAddress("recoTime", &recoTime);
   tree->SetBranchAddress("trueVertexPosX", &trueVertexPosX);
   tree->SetBranchAddress("trueVertexPosY", &trueVertexPosY);
   tree->SetBranchAddress("trueVertexPosZ", &trueVertexPosZ);
   tree->SetBranchAddress("trueVertexMomX", &trueVertexMomX);
   tree->SetBranchAddress("trueVertexMomY", &trueVertexMomY);
   tree->SetBranchAddress("trueVertexMomZ", &trueVertexMomZ);
   tree->SetBranchAddress("trueTime", &trueTime);
   tree->SetBranchAddress("hitVolume", &hitVolume);
   tree->SetBranchAddress("pValue", &pValue);
   tree->SetBranchAddress("station", &station);
   tree->SetBranchAddress("passVertexQuality", &passVertexQuality);

}
#endif