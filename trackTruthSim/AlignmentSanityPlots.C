void DrawTGraphErrors(TGraphErrors *graph, std::string title, std::string fname) {

	TCanvas *c = new TCanvas("c","c",800,600);

	graph->SetTitle(title.c_str());
	graph->GetXaxis()->SetTitleSize(.04);
	graph->GetYaxis()->SetTitleSize(.04);
	graph->GetXaxis()->SetTitleOffset(1.1);
	graph->GetYaxis()->SetTitleOffset(1.2);
	graph->GetXaxis()->CenterTitle(true);
	graph->GetYaxis()->CenterTitle(true);
	graph->GetYaxis()->SetMaxDigits(4);
	graph->SetMarkerStyle(20); //  Full circle
	graph->Draw("ALP");

	c->SaveAs((fname+".pdf").c_str());
	c->SaveAs((fname+".png").c_str());
	c->SaveAs((fname+".C").c_str());

	delete c;

	return;

}

void Run(string config) {

	TFile *fin = TFile::Open(("/pnfs/GM2/persistent/EDM/MC/dMu/Trees/Alignment/"+config+"/trackTruthTrees.1.root").c_str());

	vector<string> stn_ = {"0", "12", "18"};

	for(auto& stn : stn_) { 



		TGraphErrors *gr1 = (TGraphErrors*)fin->Get(("Extrapolation/strawGeometry/strawGeometry_worldCoordsYn_station"+stn).c_str());

		// Original graph is massive with many duplicate points.
		// Doing this speeds things up by a factor of 10 at least.
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

	Run("1mm");
	Run("0mm");

	return; 
}
