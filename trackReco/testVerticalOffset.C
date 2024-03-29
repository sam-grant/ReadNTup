double omegaAMagic = 0.00143934; // from gm2geom consts / kHz 
double g2Period = (2*TMath::Pi()/omegaAMagic) * 1e-3; // 4.3653239 us
double mMu = 105.6583715; // MeV
double aMu = 11659208.9e-10; 
double gmagic = std::sqrt( 1.+1./aMu );
double pmax = 1.01 * mMu * gmagic; // 3127.1144
double T_c = 149.2 * 1e-3; // cyclotron period [us]

void testVerticalOffset() { 


  string stns[] = {"S12", "S18", "S12S18"}; 
  int n_stn = sizeof(stns)/sizeof(stns[0]);
  int step = 125;
  int nSlices = pmax/step;

  // Vertical offset correction 
  TFile *verticalOffsetCorrectionFile = TFile::Open("correctionPlots/verticalOffsetFits_Run-1a_125MeV_BQ.root"); 
  cout<<"FILE "<<verticalOffsetCorrectionFile<<endl;

  vector<vector<TF1*>> verticalOffsetFits_; 
  
  int count = 0;

  for (int i_stn = 0; i_stn < n_stn; i_stn++) { 

    vector<TF1*> verticalOffsetFitsPerStn_;

    // Slice momentum
    for ( int i_slice = 0; i_slice < nSlices; i_slice++ ) { 

      int lo = 0 + i_slice*step; 
      int hi = step + i_slice*step;

      // Store vertical offset correction fit
      TGraphErrors* gr_verticalOffsetFit = (TGraphErrors*)verticalOffsetCorrectionFile->Get(("MomBinnedAna/"+stns[i_stn]+"_ThetaY_vs_Time_Fit_"+std::to_string(lo)+"_"+std::to_string(hi)).c_str());

      cout<<"GRAPH "<<gr_verticalOffsetFit<<endl;
 
      if(gr_verticalOffsetFit==0) verticalOffsetFitsPerStn_.push_back(0);
      else verticalOffsetFitsPerStn_.push_back(gr_verticalOffsetFit->GetFunction("DoubleExponentialFunc"));
      
      cout<<"FUNCTION "<<gr_verticalOffsetFit->GetFunction("DoubleExponentialFunc")<<endl;

      if(gr_verticalOffsetFit->GetFunction("DoubleExponentialFunc")!=0) cout<<"FUNCTION EVAL AT 1000 MEV "<<gr_verticalOffsetFit->GetFunction("DoubleExponentialFunc")->Eval(1000)<<endl;

      if(gr_verticalOffsetFit->GetFunction("DoubleExponentialFunc")!=0) count++;

      if(count==1) {
	
	TCanvas *c = new TCanvas();
	gr_verticalOffsetFit->GetFunction("DoubleExponentialFunc")->Draw();
	c->SaveAs("test.png");
	delete c;
	break;
      }
      
    }

    verticalOffsetFits_.push_back(verticalOffsetFitsPerStn_);

  }


  return;

}
