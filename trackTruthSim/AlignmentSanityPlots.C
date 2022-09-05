void DrawTGraphErrors(TGraphErrors *graph, TString title, TString fname) {

	TCanvas *c = new TCanvas("c","c",800,600);

	graph->SetTitle(title);
	graph->GetXaxis()->SetTitleSize(.04);
	graph->GetYaxis()->SetTitleSize(.04);
	graph->GetXaxis()->SetTitleOffset(1.1);
	graph->GetYaxis()->SetTitleOffset(1.2);
	graph->GetXaxis()->CenterTitle(true);
	graph->GetYaxis()->CenterTitle(true);
	graph->GetYaxis()->SetMaxDigits(4);
	graph->SetMarkerStyle(20); //  Full circle
	graph->Draw("ALP");

	c->SaveAs(fname+".pdf");
	c->SaveAs(fname+".png");
	c->SaveAs(fname+".C");

	delete c;

	return;

}

void Run(TString config, TString finName) {

	// /pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/1mm/trackTruthTrees.1.root
	// /pnfs/GM2/persistent/EDM/MC/dMu/Trees/Summer2022Campaign/sgrant/MainNtuple
	// /pnfs/GM2/persistent/EDM/MC/dMu/Trees/Summer2022Campaign/sgrant/Plus1mmNtuple
	// /pnfs/GM2/persistent/EDM/MC/dMu/Trees/Summer2022Campaign/sgrant/Minus1mmNtuple
	// /pnfs/GM2/persistent/EDM/MC/dMu/Trees/Summer2022Campaign/sgrant/Plus0.1degNtuple
	// /pnfs/GM2/persistent/EDM/MC/dMu/Trees/Summer2022Campaign/sgrant/Minus0.1degNtuple

	finName = "/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Summer2022Campaign/sgrant/"+config+"/"+finName;
	if(config=="1mm") finName = "/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/1mm/trackTruthTrees.1.root";

	TFile *fin = TFile::Open(finName);

	cout<<"---> Opened "<<finName<<", "<<fin<<endl;

	vector<TString> stn_ = {"0", "12", "18"};

	for(auto& stn : stn_) { 

		TGraphErrors *gr1 = (TGraphErrors*)fin->Get("Extrapolation/strawGeometry/strawGeometry_worldCoordsYn_station"+stn);

		// Original graph is massive with many duplicate points, loop through the first lot of points to speed things up
		TGraphErrors *gr2 = new TGraphErrors();

		for (int i(1); i<9; i++) {

			double y = gr1->Eval(i);

			// cout<<y<<endl;

			gr2->SetPoint(i-1, i, y);
			gr2->SetPointError(i-1, 0, 0);
		}

		DrawTGraphErrors(gr2, "Station "+stn+";Module number;World y-position [mm]", "Images/S"+stn+"_worldYn_"+config);

	}
	

	fin->Close();

	return;

}

void AlignmentSanityPlots() { 

	Run("1mm", "trackTruthTrees.1.root");
	Run("Plus1mmNtuple", "gm2tracker_hadd_58466221_1.root");
	Run("Minus1mmNtuple", "gm2tracker_hadd_58466222_1.root");
	Run("Plus0.1degNtuple", "gm2tracker_hadd_36355304_1.root");
	Run("Minus0.1degNtuple", "gm2tracker_hadd_58466326_1.root");

	return; 
}
