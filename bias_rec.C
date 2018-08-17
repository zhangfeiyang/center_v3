{

	const int n = 70;
	double pe[n],epe[n],tmp,e[n],ratio[n],sigma[n],esigma[n];
	double res[n],eres[n];
	double res2[n],bias[n],ebias[n];
//	TF1 *fun0 = new TF1("fun0","sqrt(2.86*2.86/x+0.71963*0.71963)",0.01,8);
//	TF1 *fun0 = new TF1("fun0","sqrt(2.993*2.993/x+0.7276*0.7276)",0.01,8);
    TF1 *fun0 = new TF1("fun0","sqrt([0]*[0]/x+[1]*[1])",0.01,8);
    ifstream fin0("result");
    double a,b;
    fin0>>tmp>>a>>b>>tmp>>tmp;
    fun0->SetParameters(a*100,b*100);
	ifstream fin("ele_data");
	for(int i=0;i<n;i++){
		fin>>pe[i]>>epe[i]>>tmp>>sigma[i]>>esigma[i];
		e[i] = (1+i)/10.0;
		res[i] = 100*sigma[i]/pe[i];
		eres[i] = 100*esigma[i]/pe[i];
		pe[i]/=1.3056438e3;
		bias[i] = (res[i]-fun0->Eval(pe[i]));
		ebias[i] = eres[i];
	}
	
	TGraphErrors *T = new TGraphErrors(n,pe,bias,0,ebias);
//	TGraphErrors *T = new TGraphErrors(n,pe,res,0,eres);
		
	T->GetYaxis()->SetRangeUser(-0.1,0.05);
	T->GetXaxis()->SetRangeUser(1,8);
    	T->SetMarkerStyle(20);
    	T->SetMarkerSize(1);
	T->GetXaxis()->SetTitle("E_{rec} [MeV]");
	T->GetYaxis()->SetTitle("bias of resolution [%]");
    	T->Draw();
//	fun0->Draw("same");
//
//    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
//    l->AddEntry(T,"data","p");
//    l->AddEntry(fun0,"prediction","l");
//    l->Draw();

}
