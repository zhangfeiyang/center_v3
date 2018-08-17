
{

	const int n = 8;
	double energy[n],PE[n],ePE[n],sigma[n],esigma[n];
	double energy0[n],PE0[n],ePE0[n],sigma0[n],esigma0[n];
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
	}


	TGraphErrors *T0 = new TGraphErrors(n,energy,sigma0,0,esigma0);
	TGraphErrors *T1 = new TGraphErrors(n,energy,sigma,0,esigma);
	T1->Draw("Ap");
	T1->SetMarkerStyle(20);
    T1->SetMarkerSize(1);
	T1->Draw();
	T0->SetMarkerColor(kRed);
	T0->Draw("psame");

}
