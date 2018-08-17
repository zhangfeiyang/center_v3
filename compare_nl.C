
{

	const int n = 8;
	double energy[n],PE[n],ePE[n],sigma[n],esigma[n];
	double energy0[n],PE0[n],ePE0[n],sigma0[n],esigma0[n];
	double bias[n],ebias[n];
	ifstream fin("data");
	ifstream fin0("data0");

	for(int i=0;i<n;i++){
			fin>>energy[i]>>PE[i]>>ePE[i]>>sigma[i]>>esigma[i];
			fin0>>energy0[i]>>PE0[i]>>ePE0[i]>>sigma0[i]>>esigma0[i];
			if(i==0 || i==3){
				sigma[i]/=sqrt(2);
				sigma0[i]/=sqrt(2);
				PE[i] /= 2;
				PE0[i] /= 2;
			}
			sigma[i] = sigma[i]/PE[i];
			sigma0[i] = sigma0[i]/PE0[i];
			esigma[i] /= PE[i];
			esigma0[i] /= PE0[i];
			sigma[i] *=100;
			sigma0[i] *=100;
			esigma[i] *=100;
			esigma0[i] *=100;

			PE0[i]/=1.3056438e3;
			ePE0[i]/=1.3056438e3;
			PE[i]/=1.3056438e3;
			ePE[i]/=1.3056438e3;

			/*
			PE[i]/=1.3039277466e3;
			ePE0[i]/=1.3056438e3;
			ePE[i]/=1.3039277466e3;
			*/

			PE0[i]/=energy[i];
			PE[i]/=energy[i];
			ePE0[i]/=energy[i];
			ePE[i]/=energy[i];
			
			bias[i] = 100*(PE0[i] - PE[i]);
			ebias[i] = 100*sqrt(ePE0[i]**2 + ePE[i]**2);
			//cout << TMath::Abs(bias[i]) <<"\n";
			cout << (bias[i]) <<"\n";
			
	}
	//TCanvas *c1 = new TCanvas("","",1600,600);
	//TCanvas *c1 = new TCanvas();

	//c1->Divide(2,1);
	//c1->cd(1);
	TGraphErrors *T0 = new TGraphErrors(n,energy,PE0,0,ePE0);
	TGraphErrors *T1 = new TGraphErrors(n,energy,PE,0,ePE);
	T1->SetMarkerStyle(20);
    T1->SetMarkerSize(1);
	//T1->Draw();
	T0->SetMarkerColor(kRed);
	//T0->Draw("same");

	TLegend *l = new TLegend(0.6,0.6,0.8,0.8);
	l->AddEntry(T0,"ideal case","pl");
	l->AddEntry(T1,"realistic case","pl");
	l->Draw();

	//c1->cd(2);
	TGraphErrors *t = new TGraphErrors(n,energy,bias,0,ebias);
	t->Draw();
}
