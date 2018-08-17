#include "rootheader.h"
TTree *t_K40;
TTree *t_n_H;
TTree *t_n_Fe;
TTree *t_n_C;
TTree *t_Am_C;
TTree *t_Ge68;
TTree *t_Cs137;
TTree *t_Co60;
TTree *t_Mn54;
TTree *t;

const int N = 10000;
double totalPE[N];

//double res_gamma2e[9] = {0.0196429,0.0161949,0.0145086,0.01355292,0.0134099,0.0128792,0.00841269,0.00729792,0.00649169};
double get_ratio(double *X,double *p){

		double prob[5000],energy[5000];	
		double xx = X[0];
		double scale;
		int id;
	
		double p0 = p[3];	
		double p1 = p[4];	
		double p2 = p[5];	
		double p3 = p[6];	
		double p4 = p[7];	

		id = 0;
        if(xx<(0.6617+0.511)/2)                        {id = 2;t = t_Ge68;}
        if(xx>(0.6617+0.511)/2 && xx<(0.6617+0.835)/2) {t = t_Cs137;}
        if(xx>(0.6617+0.835)/2 && xx<(1.2500+0.835)/2) {t = t_Mn54;}
        if(xx>(1.2500+0.835)/2 && xx<(1.4608+1.250)/2) {id = 2;t = t_Co60;}
        if(xx>(1.4608+1.250)/2 && xx<(1.4608+2.223)/2) {t = t_K40; }
        if(xx>(1.4608+2.223)/2 && xx<(4.9450+2.223)/2) {t = t_n_H; }
        if(xx>(4.9450+2.223)/2 && xx<(4.9450+6.130)/2) {t = t_n_C; }
        if(xx>(4.9450+6.130)/2 && xx<(6.1300+7.500)/2) {t = t_Am_C;}
        if(xx>(6.1300+7.500)/2)                        {t = t_n_Fe;}

		int nEle;
		double E_e;
		double abs_e;
		double sign;
		double Erec;
		double sum_sigma_2 = 0;

		for(int i=0;i<N;i++){
			t->GetEntry(i+2*N);
			totalPE[i] = 0;
			nEle = t->GetLeaf("nEle")->GetValue(0);
			for(int j=0;j<nEle;j++){
				E_e = t->GetLeaf("E_e")->GetValue(j);
				abs_e = TMath::Abs(E_e);
				if(abs_e==0) {
					Erec = (abs_e*p0+p1*abs_e*abs_e+p2)/(1+p3*exp(-p4*abs_e));
				}else{
					sign = E_e/abs_e;
					Erec = sign*(abs_e*p0+p1*abs_e*abs_e+p2)/(1+p3*exp(-p4*abs_e));
				}

				totalPE[i] += Erec;
				if(Erec>0.2){
					sum_sigma_2 += p[1]*p[1]*Erec + p[2]*p[2]*Erec*Erec;
				}else{
					sum_sigma_2 += p[0]*p[0]*Erec;
				}

			}
		}

		double sigma_2 = sum_sigma_2/N;
		double rms = TMath::RMS(N,totalPE);
		double mean= TMath::Mean(N,totalPE);

		if(id==2){
			sigma_2/=2;
			rms/=sqrt(2);
			mean/= 2;
		}

		double sigma = sqrt(sigma_2+rms*rms);
		//cout << sigma/mean <<"\n";
		return sigma/mean;			
}

void fit3_cor2(){

    TFile *file2_K40   = new TFile("K40.root","read");
    TFile *file2_n_H   = new TFile("n_H.root","read");
    TFile *file2_n_Fe  = new TFile("n_Fe.root","read");
    TFile *file2_n_C   = new TFile("n_C.root","read");
    TFile *file2_Am_C  = new TFile("Am_C.root","read");
    TFile *file2_Ge68  = new TFile("Ge68.root","read");
    TFile *file2_Cs137 = new TFile("Cs137.root","read");
    TFile *file2_Co60  = new TFile("Co60.root","read");
    TFile *file2_Mn54  = new TFile("Mn54.root","read");

    t_K40  = (TTree*)file2_K40->Get("Gamma2e");
    t_n_H  = (TTree*)file2_n_H->Get("Gamma2e");
    t_n_Fe = (TTree*)file2_n_Fe->Get("Gamma2e");
    t_n_C  = (TTree*)file2_n_C->Get("Gamma2e");
    t_Am_C = (TTree*)file2_Am_C->Get("Gamma2e");
    t_Ge68 = (TTree*)file2_Ge68->Get("Gamma2e");
    t_Cs137= (TTree*)file2_Cs137->Get("Gamma2e");
    t_Co60 = (TTree*)file2_Co60->Get("Gamma2e");
    t_Mn54 = (TTree*)file2_Mn54->Get("Gamma2e");

	const int n = 8;
	double energy[n],PE[n],ePE[n],sigma[n],tmp;
	ifstream fin("data");
	double esigma[n]={0.0594764,
0.00499388,
-0.031551,
-0.0110508,
-0.0300592,
-0.00513058,
-0.00897476,
-0.0617164};

	//for(int i=0;i<n;i++) esigma[i] = TMath::Abs(esigma[i])/100.0;

	for(int i=0;i<n;i++) esigma[i] = (esigma[i])/100.0;
	for(int i=0;i<n;i++){
			fin>>energy[i]>>PE[i]>>ePE[i]>>sigma[i]>>tmp;
			//fin>>energy[i]>>PE[i]>>ePE[i]>>sigma[i]>>esigma[i];
			if(i==0 || i==3){
				sigma[i]/=sqrt(2);
				PE[i] /= 2;
			}
			sigma[i] = sigma[i]/PE[i];
			sigma[i] += esigma[i];
			cout << sigma[i] <<"\t"<< esigma[i]<<"\n";
			//esigma[i] /= PE[i];
	}
	for(int i=0;i<n;i++) esigma[i] = TMath::Abs(esigma[i]);

	TGraphErrors *T1 = new TGraphErrors(n,energy,sigma,0,esigma);
	T1->Draw("Ap");
	T1->SetMarkerStyle(20);
    T1->SetMarkerSize(1);

	TF1 *f1 = new TF1("f1",get_ratio,0,8,8);
	f1->SetParameters(0.027,0.02803,0.007886,1.07476,0.00148426,-0.001801,0.090052,1.30358);

	f1->SetParNames("a_{0}","a","b","P_{0}","P_{1}","P_{2}","P_{3}","P_{4}");

    f1->FixParameter(3,1.07476);
    f1->FixParameter(4,0.00148426);
    f1->FixParameter(5,-0.001801);
    f1->FixParameter(6,0.090052);
    f1->FixParameter(7,1.30358);

	T1->Fit("f1","0","",0,8);
	ofstream fout("result");
	double *a = f1->GetParameters();
	fout<<a[0]<<"\t"<<a[1]<<"\t"<<a[2]<<"\t"<<a[3]<<"\t"<<a[4]<<"\n";

    double fit[n];
    for(int i=0;i<n;i++) {
        fit[i] = f1->Eval(energy[i]);
    }
    TGraph *T2 = new TGraph(n,energy,fit);
    T2->SetMarkerColor(kRed);
	T1->SetMarkerStyle(20);
    T1->SetMarkerSize(1);
	T2->SetMarkerStyle(20);
    T2->SetMarkerSize(1);
    T2->Draw("psame");

	T1->GetYaxis()->SetTitle("resolution");
	T1->GetXaxis()->SetTitle("E_{true} [MeV]");

	TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
	l->AddEntry(T1,"data","p");
	l->AddEntry(T2,"prediction","p");
	l->Draw();

    gStyle->SetStatY(0.7);
    gStyle->SetStatX(0.9);
    gStyle->SetStatH(0.18);
    gStyle->SetStatW(0.18);
}
