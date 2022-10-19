# ReadNTup

**The purpose of this repository is the read in data/simulation in the form of ROOT Ntuples and output analysis level histograms.** 

- The core of each directory is a Plotter.C, a Plotter.h, a Makefile, and submit.sh/run.sh scripts. 

- Data is read in by **trackReco**. 

- EDM simulation is primarily read by **decayTruthSim** (*all decays*) and **trackTruthSim** (*reco/truth vertices*), older types of Ntuples are read by trackRecoSim.  

- Acceptance simulation is read by **trackerAcceptance**.

- Tracker resolution simulation, which is the same dataset as is used for acceptance studies, is read by **trackerResolution2**.
