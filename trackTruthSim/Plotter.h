#ifndef Plotter_h
#define Plotter_h

// ROOT includes
#include "TTree.h"

using namespace std;

class InitBranches { 

public: 

   // Declaration of leaf types
   unsigned int runNum;
   unsigned int subRunNum;
   unsigned int eventNum;
   float recoVertexPosX;
   float recoVertexPosY;
   float recoVertexPosZ;
   float trueVertexPosX;
   float trueVertexPosY;
   float trueVertexPosZ;
   float recoVertexMomX;
   float recoVertexMomY;
   float recoVertexMomZ;
   float trueVertexMomX;
   float trueVertexMomY;
   float trueVertexMomZ;
   float recoTime;
   float trueTime;
   bool hitVolume;
   float pValue;
   int station;
   bool passVertexQuality;

   // Declare constructer
   InitBranches(TTree *tree);

};

// Constructer
InitBranches::InitBranches(TTree* tree) {

   tree->SetBranchAddress("runNum", &runNum);
   tree->SetBranchAddress("subRunNum", &subRunNum);
   tree->SetBranchAddress("eventNum", &eventNum);
   tree->SetBranchAddress("recoVertexPosX", &recoVertexPosX);
   tree->SetBranchAddress("recoVertexPosY", &recoVertexPosY);
   tree->SetBranchAddress("recoVertexPosZ", &recoVertexPosZ);
   tree->SetBranchAddress("trueVertexPosX", &trueVertexPosX);
   tree->SetBranchAddress("trueVertexPosY", &trueVertexPosY);
   tree->SetBranchAddress("trueVertexPosZ", &trueVertexPosZ);
   tree->SetBranchAddress("recoVertexMomX", &recoVertexMomX);
   tree->SetBranchAddress("recoVertexMomY", &recoVertexMomY);
   tree->SetBranchAddress("recoVertexMomZ", &recoVertexMomZ);
   tree->SetBranchAddress("trueVertexMomX", &trueVertexMomX);
   tree->SetBranchAddress("trueVertexMomY", &trueVertexMomY);
   tree->SetBranchAddress("trueVertexMomZ", &trueVertexMomZ);
   tree->SetBranchAddress("recoTime", &recoTime);
   tree->SetBranchAddress("trueTime", &trueTime);
   tree->SetBranchAddress("hitVolume", &hitVolume);
   tree->SetBranchAddress("pValue", &pValue);
   tree->SetBranchAddress("station", &station);
   tree->SetBranchAddress("passVertexQuality", &passVertexQuality);

}

#endif