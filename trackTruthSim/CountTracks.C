{ 
  TFile *f = TFile::Open("count_truth.root");

  TH1D *h1 = (TH1D*)f->Get("Momentum_Raw");//SimultaneousAnalysis/ThetaY_vs_Time_Modulo");
  TH1D *h2 = (TH1D*)f->Get("Momentum_Cuts");//SimultaneousAnalysis/ThetaY_vs_Time_Modulo");

  int N1 = h1->GetEntries();
  int N2 = h2->GetEntries();

  cout<<"raw = "<<N1<<endl;
  cout<<"cuts = "<<N2<<endl;

  f->Close();

}
