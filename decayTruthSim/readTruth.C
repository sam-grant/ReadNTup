void readTruth(TTree* in);

void readTruth(){
  TChain* chain = new TChain("trackerPhaseAnalyzer/g2phase");
  
  vector<int> IDs = {
    1605115845,1605119474,1605120404,1605121195,1605121936,1605163551,1605163652,1605175087,1605175090,1605175281,1605183372,1605213661,1605214193,1605214526,1605214527,1605253108,1605254292,1605273200,1605273201,1605273210,1605273221,1605273361,1605289885,1605289892,1605289896,1605289899,1605290026,1605355172,1605358723,1605359386,1605359387,1605377164,1605377790,1605452425,1605453447,1605477504,1605477505,1605477506,1605477691,1605537341,1605583477,1605583743,1605605585,1605605586,1605626943,1605626979,1605639348,1605639460,1605645669,1605645992,1605665438,1605666878,1605670973,1605671194,1605673964,1605674177,1605676962,1605677375,1605681672,1605682642,1605684045,1605685936,1605686123,1605687615,1605687616,1605768429,1605768430,1605812225,1605812417,1605836303,1605836538,1605844930,1605845460,1605853783,1605909359,1605936909,1606050012,1606052269,1606059674,1606073213,1606075298,1606100078,1606102629,1606198565,1606226499,1606296559,1606317405,1606403666,1606422342,1606478184,1606484866,1606517047,1606522670,1606563991,1606564442,1606566556,1606566572,1606568092,1606586696,1606589186,1606662932,1606663138,1606663271,1606663634,1606665149,1606665165,1606678837,1606679465,1606680017,1606680298,1606680377,1606688217,1606688317,1606695317,1606695624,1606695913,1606696218,1606699172,1606706168,1606710620,1606710914,1606711818,1606715079,1606725921,1606726835,1606727563,1606727578,1606729379,1606729778,1606732777,1606741571,1606742847,1606743485,1606744203,1606749526,1606777484,1606777487,1606777680,1606778127,1606792830,1606793238,1606793443,1606793718,1606793963,1606808416,1606809019,1606809313,1606809632,1606809635,1606824318,1606824623,1606824921,1606825222,1606825514,1606941377,1606941687,1606941872,1607226338,1607226592,1607226593,1607241916,1607242220,1607242221,1607242516,1607242522,1607257601,1607257685,1607257834,1607258115,1607264746,1607273584,1607274975,1607275034,1607275187,1607282029,1607289522,1607291077,1607291119,1607292044,1607297673,1607305642,1607307018,1607307400,1607308194,1607488066,1607488067,1607503280,1607503781,1607529891,1607529929,1607529930,1607530602
  };

  //for 1 file check
  IDs.clear();
  IDs.push_back(1607530602);

  int maxFiles = 10;
  int nFiles = 0;
  for (auto& ID: IDs){
    if (nFiles == maxFiles) break;
    chain->Add(Form("/pnfs/GM2/persistent/Tracking/mc_gasgun_run1_2001/TrackerPhaseAnalyzerTree_2001_%i.root", ID));
    nFiles++;
  }
  readTruth((TTree*) chain);
}

void readTruth(TTree* tree){

  TCanvas* c1 = new TCanvas("c1", "", 800, 600);
  tree->Print();

  int64_t nEntries = tree->GetEntries();
  nEntries = 5E6;
  cout << "\n\tEntries = " << nEntries << ";\n" << endl;

  float vertexTime, trackMom, py, posy, trackY;
  int station;
  bool passCandQuality, passTrackQuality, passVertexQuality;
  tree->SetBranchAddress("decayTime",&vertexTime);
  tree->SetBranchAddress("trackMomentum",&trackMom);
  tree->SetBranchAddress("trackMomentumY",&py);
  tree->SetBranchAddress("trackPosY",&trackY);
  tree->SetBranchAddress("decayVertexPosY",&posy);
  tree->SetBranchAddress("station",&station);
  tree->SetBranchAddress("passCandidateQuality",&passCandQuality);
  tree->SetBranchAddress("passTrackQuality",&passTrackQuality);
  tree->SetBranchAddress("passVertexQuality",&passVertexQuality);

  //truth
  float decayTime;
  float posiPX, posiPY, posiPZ, posiP;
  float posiX, posiY, posiZ;
  tree->SetBranchAddress("muDecayTime",&decayTime);
  tree->SetBranchAddress("posiInitPosX",&posiX);
  tree->SetBranchAddress("posiInitPosY",&posiY);
  tree->SetBranchAddress("posiInitPosZ",&posiZ);
  tree->SetBranchAddress("posiInitP",&posiP);
  tree->SetBranchAddress("posiInitPX",&posiPX);
  tree->SetBranchAddress("posiInitPY",&posiPY);
  tree->SetBranchAddress("posiInitPZ",&posiPZ);

  //omega a period
  double period = 4.365411; //us

  TFile* fout = new TFile("truthOut.root", "RECREATE");
  TDirectory* noTracksDir = fout->mkdir("noTracks"); 
  TDirectory* stat12Dir = fout->mkdir("Station12"); 
  TDirectory* stat18Dir = fout->mkdir("Station18"); 
  TDirectory* stat12RecoDir = fout->mkdir("RecoStation12"); //reconstructed track passes vertex quality
  TDirectory* stat18RecoDir = fout->mkdir("RecoStation18"); //reconstructed track passes vertex quality

  for (auto& dir: {stat12Dir, stat18Dir, stat12RecoDir, stat18RecoDir}){
    //make a directory for each different y bin
    for (int yBin(0); yBin < nYBins; yBin++){
      TDirectory* yDir = dir->mkdir(Form("ybin_%i", yBin));
      yDir->cd();
      TProfile* tp = new TProfile("mom_vs_vertAngle", ";momentum [MeV]; average vertical angle [mrad]",29, 300, 3200);
    }
  }

  for (auto &statDir : {noTracksDir, stat12Dir, stat18Dir, stat12RecoDir, stat18RecoDir}){
    statDir->cd();
    //book hists

    //vs momentum plots
    TH2F* mom_vs_y = new TH2F("mom_vs_y", ";momentum [MeV]; y position [mm]", 29, 300, 3200, 44, -55, 55);
    TH2F* mom_vs_vertAngle = new TH2F("mom_vs_vertAngle", ";momentum [MeV]; vertical angle [mrad]", 29, 300, 3200, 50, -100, 100);
  }

  double targetPerc = 0;
  for(int64_t entry = 0; entry < nEntries; entry++){
    if(100*float(entry) / nEntries > targetPerc){
      cout << Form("Processed %.1f%%", 100*float(entry)/nEntries) << endl;
      targetPerc += 10;
    }
    tree->GetEntry(entry);
    
    TDirectory* statDir = (station ==12)? stat12Dir : (station ==18)? stat18Dir: noTracksDir;
    
    //get y bin
    //posiY+50 / ; 
    

    double angle = posiPY / posiP;
    angle *= 1E6; //microRad
    if (decayTime > 30.6){
      ((TH2F*) statDir->Get("mom_vs_y"))->Fill(posiP, posiY);
      ((TH2F*) statDir->Get("mom_vs_vertAngle"))->Fill(posiP, angle*1E-3); //mrad

      if (passVertexQuality) {
	TDirectory* recoDir = (station ==12)? stat12RecoDir : (station ==18)? stat18RecoDir : NULL;
	((TH2F*) recoDir->Get("mom_vs_y"))->Fill(posiP, posiY);
	((TH2F*) recoDir->Get("mom_vs_vertAngle"))->Fill(posiP, angle*1E-3); //mrad
      }
    }
  }

  fout->Write();
  fout->Close();

}
