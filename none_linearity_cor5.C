#include "rootheader.h"
TH1D *h_K40;
TH1D *h_n_H;
TH1D *h_n_Fe;
TH1D *h_n_C;
TH1D *h_Am_C;
TH1D *h_Ge68;
TH1D *h_Cs137;
TH1D *h_Co60;
TH1D *h_Mn54;
TH1D *h;
TH1F* h_scale_K40;  
TH1F* h_scale_n_H;  
TH1F* h_scale_n_Fe;
TH1F* h_scale_n_C;  
TH1F* h_scale_Am_C;
TH1F* h_scale_Ge68;
TH1F* h_scale_Cs137;
TH1F* h_scale_Co60;
TH1F* h_scale_Mn54; 
TH1F *h_scale;
double get_ratio(double *X,double *p){
	
		double prob[5000],energy[5000];
		double xx = X[0];
        double scale;
        if(xx<(0.6617+0.511)/2)                        {h = h_Ge68;h_scale =  h_scale_K40;  }
        if(xx>(0.6617+0.511)/2 && xx<(0.6617+0.835)/2) {h = h_Cs137;h_scale = h_scale_n_H;  }
        if(xx>(0.6617+0.835)/2 && xx<(1.2500+0.835)/2) {h = h_Mn54;h_scale =  h_scale_n_Fe; }
        if(xx>(1.2500+0.835)/2 && xx<(1.4608+1.250)/2) {h = h_Co60;h_scale =  h_scale_n_C;  }
        if(xx>(1.4608+1.250)/2 && xx<(1.4608+2.223)/2) {h = h_K40;h_scale =   h_scale_Am_C; }
        if(xx>(1.4608+2.223)/2 && xx<(4.9450+2.223)/2) {h = h_n_H;h_scale =   h_scale_Ge68; }
        if(xx>(4.9450+2.223)/2 && xx<(4.9450+6.130)/2) {h = h_n_C;h_scale =   h_scale_Cs137;}
        if(xx>(4.9450+6.130)/2 && xx<(6.1300+7.500)/2) {h = h_Am_C;h_scale =  h_scale_Co60; }
        if(xx>(6.1300+7.500)/2)                        {h = h_n_Fe;h_scale =  h_scale_Mn54; }

		scale = h_scale->GetMean();
		double abs_e;
		double s1=0,s2=0;
		for(int i=0;i<5000;i++) {
				energy[i]=h->GetBinCenter(i+1);
				prob[i]=h->GetBinContent(i+1);
				s2 + = energy[i]*prob[i]; 
				abs_e = TMath::Abs(energy[i]);
				s1 + = energy[i]*prob[i]*(p[0]+p[1]*abs_e+p[2]/abs_e)/(1+p[3]*exp(-p[4]*abs_e)); 
		}
	    return s1/s2;
};


void none_linearity_cor5(){

    TFile *file_K40   = new TFile("K40_pdf.root","read");
    TFile *file_n_H   = new TFile("n_H_pdf.root","read");
    TFile *file_n_Fe  = new TFile("n_Fe_pdf.root","read");
    TFile *file_n_C   = new TFile("n_C_pdf.root","read");
    TFile *file_Am_C  = new TFile("Am_C_pdf.root","read");
    TFile *file_Ge68  = new TFile("Ge68_pdf.root","read");
    TFile *file_Cs137 = new TFile("Cs137_pdf.root","read");
    TFile *file_Co60  = new TFile("Co60_pdf.root","read");
    TFile *file_Mn54  = new TFile("Mn54_pdf.root","read");

    h_K40  = (TH1D*)file_K40->Get("h1");
    h_n_H  = (TH1D*)file_n_H->Get("h1");
    h_n_Fe = (TH1D*)file_n_Fe->Get("h1");
    h_n_C  = (TH1D*)file_n_C->Get("h1");
    h_Am_C = (TH1D*)file_Am_C->Get("h1");
    h_Ge68 = (TH1D*)file_Ge68->Get("h1");
    h_Cs137= (TH1D*)file_Cs137->Get("h1");
    h_Co60 = (TH1D*)file_Co60->Get("h1");
    h_Mn54 = (TH1D*)file_Mn54->Get("h1");

    h_scale_K40  = (TH1F*)file_K40->Get("h0");
    h_scale_n_H  = (TH1F*)file_n_H->Get("h0");
    h_scale_n_Fe = (TH1F*)file_n_Fe->Get("h0");
    h_scale_n_C  = (TH1F*)file_n_C->Get("h0");
    h_scale_Am_C = (TH1F*)file_Am_C->Get("h0");
    h_scale_Ge68 = (TH1F*)file_Ge68->Get("h0");
    h_scale_Cs137= (TH1F*)file_Cs137->Get("h0");
    h_scale_Co60 = (TH1F*)file_Co60->Get("h0");
    h_scale_Mn54 = (TH1F*)file_Mn54->Get("h0");


	const int n = 8;
	double energy[n],r[n],er[n],Erec[n],eErec[n];
	double pe[n],epe[n],sigma[n],esigma[n];
	ifstream fin("data");
	double scale = 1.3056438e3;
	//double scale = 1.3039277466e3;
	//ifstream fin("data2");
	//double er[n] = {0.00214208,0.000226145,0.000516776,6.5722e-09,0.000210549,0.000689992,0.000484337,0.00445041};
	double er[n] = {0.00337988,0.00148852,0.00179825,0.00131434,0.00153622,0.000664956,0.000898696,0.00304613};

	//for(int i=0;i<n;i++)
	//	er[i] *= 2;
	for(int i=0;i<n;i++){
			fin>>energy[i]>>pe[i]>>epe[i]>>sigma[i]>>esigma[i];
			Erec[i] = pe[i]/scale;	
			eErec[i] = epe[i]/scale;	
			if(i==0 || i==3) {
				Erec[i] /= 2;
			}
			r[i] = Erec[i]/energy[i];
			if(i<n-1) r[i]+=er[i];
			else r[i]-=er[i];
			//er[i] = eErec[i]/energy[i];
			cout<<r[i]<<"\t"<<er[i]<<endl;
	}
	TGraphErrors *T1 = new TGraphErrors(n,energy,r,0,er);
	//TGraphErrors *T1 = new TGraphErrors(n,Erec,r,0,er);
	T1->Draw("Ap");
	T1->SetMarkerStyle(20);
    T1->SetMarkerSize(0.8);

	TF1 *f1 = new TF1("f1",get_ratio,0,7,5);
	//f1->SetParameters(1.05,0.002,0.11,1.5);
	//f1->SetParameters(1.043762,0.005738488,0.123178,2.46286);
	//f1->SetParameters(1.05376,0.00277669,0.166229,3.00274);
	f1->SetParameters(1.05412,0.0069017,0.00938879,0.170967,2.95026);
	//T1->Fit("f1","","",0,1.2);
	T1->Fit("f1","","",0,8);
	ofstream fout("result");
	double *a = f1->GetParameters();
	fout<<a[0]<<"\t"<<a[1]<<"\t"<<a[2]<<"\t"<<a[3]<<"\t"<<a[4]<<"\n";
	
	double fit[n];
	for(int i=0;i<n;i++) {
		fit[i] = f1->Eval(energy[i]);
	}
	TGraph *T2 = new TGraph(n,energy,fit);
	T2->SetMarkerColor(kYellow);
	T2->Draw("psame");
	//f->Draw("same");
	
	gStyle->SetStatY(0.7);
	gStyle->SetStatX(0.9);
	gStyle->SetStatH(0.18);
	gStyle->SetStatW(0.18);

	//double th_r[n]={0.961133,0.968143,0.97784,0.994842,0.997859,1.02436,1.05276};
	////TGraph *T2 = new TGraph(n,energy,th_r);
	////T2->Draw("same");
	//TLegend *l = new TLegend(0.1,0.2,0.3,0.4);
	//l->AddEntry(T,"gamma energy non-linearity","l");
	//l->AddEntry(f,"electron energy non-linearity","l");
	//l.Draw();
}
