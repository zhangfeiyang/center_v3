{

	TFile *f = new TFile("n_C_pdf.root","read");
	TH1F *h1 = (TH1F*)f->Get("h1");
	
	for(int i=2400;i<2600;i++){
		cout<< h1->GetBinContent(i)<<"\n";
	}
	
}
