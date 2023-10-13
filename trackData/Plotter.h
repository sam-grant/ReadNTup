#ifndef Plotter_h
#define Plotter_h

// ROOT includes
#include "TTree.h"

using namespace std;

class InitBranches { 

public: 

   // Declaration of leaf types
   unsigned int midasSerialNum;
   unsigned int runNum;
   unsigned int subRunNum;
   unsigned int eventNum;
   unsigned int islandNum;
   int inFillLaserCount;
   bool passTrackQuality;
   bool passCandidateQuality;
   bool passDecayVertexQuality; 
   bool hitVolume;
   bool hasDecayVertex;
   bool hasCaloVertex;
   int station;
   int nHits;
   int nUHits;
   int nVHits;
   float missedLayersFrac;
   float trackPValue;
   float trackT0;
   float trackMomentum;
   float trackMomentumX;
   float trackMomentumY;
   float trackMomentumZ;
   float decayTime;
   float decayVertexPosX;
   float decayVertexPosY;
   float decayVertexPosZ;
   float decayVertexPosR;
   float decayVertexMom;
   float decayVertexMomX;
   float decayVertexMomY;
   float decayVertexMomZ;
   float decayVertexMomUnc;
   float decayVertexUncR;
   float decayVertexUncY;
   float decayVertexUncPR;
   float decayVertexUncPY;
   float minDriftTime;
   float maxDriftTime;
   float maxResidual;
   float decayExtrapolatedDistance;
   float caloExtrapolatedDistance;
   // Calo branches not included 

   // Declare constructer
   InitBranches(TTree *tree);

};

// Constructer
InitBranches::InitBranches(TTree* tree) {

   tree->SetBranchAddress("midasSerialNum", &midasSerialNum);
   tree->SetBranchAddress("runNum", &runNum);
   tree->SetBranchAddress("subRunNum", &subRunNum);
   tree->SetBranchAddress("eventNum", &eventNum);
   tree->SetBranchAddress("islandNum", &islandNum);
   tree->SetBranchAddress("inFillLaserCount", &inFillLaserCount);
   tree->SetBranchAddress("passTrackQuality", &passTrackQuality);
   tree->SetBranchAddress("passCandidateQuality", &passCandidateQuality);
   tree->SetBranchAddress("passDecayVertexQuality", &passDecayVertexQuality);
   tree->SetBranchAddress("hitVolume", &hitVolume);
   tree->SetBranchAddress("hasDecayVertex", &hasDecayVertex);
   tree->SetBranchAddress("hasCaloVertex", &hasCaloVertex);
   tree->SetBranchAddress("station", &station);
   tree->SetBranchAddress("nHits", &nHits);
   tree->SetBranchAddress("nUHits", &nUHits);
   tree->SetBranchAddress("nVHits", &nVHits);
   tree->SetBranchAddress("missedLayersFrac", &missedLayersFrac);
   tree->SetBranchAddress("trackPValue", &trackPValue);
   tree->SetBranchAddress("trackT0", &trackT0);
   tree->SetBranchAddress("trackMomentum", &trackMomentum);
   tree->SetBranchAddress("trackMomentumX", &trackMomentumX);
   tree->SetBranchAddress("trackMomentumY", &trackMomentumY);
   tree->SetBranchAddress("trackMomentumZ", &trackMomentumZ);
   tree->SetBranchAddress("decayTime", &decayTime);
   tree->SetBranchAddress("decayVertexPosX", &decayVertexPosX);
   tree->SetBranchAddress("decayVertexPosY", &decayVertexPosY);
   tree->SetBranchAddress("decayVertexPosZ", &decayVertexPosZ);
   tree->SetBranchAddress("decayVertexPosR", &decayVertexPosR);
   tree->SetBranchAddress("decayVertexMom", &decayVertexMom);
   tree->SetBranchAddress("decayVertexMomX", &decayVertexMomX);
   tree->SetBranchAddress("decayVertexMomY", &decayVertexMomY);
   tree->SetBranchAddress("decayVertexMomZ", &decayVertexMomZ);
   tree->SetBranchAddress("decayVertexMomUnc", &decayVertexMomUnc);
   tree->SetBranchAddress("decayVertexUncR", &decayVertexUncR);
   tree->SetBranchAddress("decayVertexUncY", &decayVertexUncY);
   tree->SetBranchAddress("decayVertexUncPR", &decayVertexUncPR);
   tree->SetBranchAddress("decayVertexUncPY", &decayVertexUncPY);
   tree->SetBranchAddress("minDriftTime", &minDriftTime);
   tree->SetBranchAddress("maxDriftTime", &maxDriftTime);
   tree->SetBranchAddress("maxResidual", &maxResidual);
   tree->SetBranchAddress("decayExtrapolatedDistance", &decayExtrapolatedDistance);
   tree->SetBranchAddress("caloExtrapolatedDistance", &caloExtrapolatedDistance);

}

#endif


