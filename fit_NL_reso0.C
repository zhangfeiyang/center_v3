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

double init_PE,init_sigma;
double get_ratio1(double *X,double *p){
	
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
				s2 += energy[i]*prob[i]; 
				abs_e = TMath::Abs(energy[i]);
				s1 += energy[i]*prob[i]*(p[0]+p[1]*abs_e+p[2]/abs_e)/(1+p[3]*exp(-p[4]*abs_e)); 
		}
	    return s1/s2;
};
double get_ratio2(double *X,double *p){

        double prob[5000],energy[5000];
        double xx = X[0];
        double scale;
        int id;

        double p0 = p[2];
        double p1 = p[3];
        double p2 = p[4];
        double p3 = p[5];
        double p4 = p[6];

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
                Erec += init_PE/1.3056438e3;
                sum_sigma_2 += p[0]*p[0]*Erec + p[1]*p[1]*Erec*Erec - init_sigma*init_sigma/1.3056438e3/1.3056438e3;
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

        return sigma/mean;

}
        

void fit_NL_reso0(){

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
	double energy[n],r[n],Erec[n],eErec[n];
	double pe[n],epe[n],sigma[n],esigma[n];
	ifstream fin("data0");
	double scale = 1.3056438e3;
	double er[n] = {0.00337988,0.00148852,0.00179825,0.00131434,0.00153622,0.000664956,0.000898696,0.00304613};

	for(int i=0;i<n;i++){
			fin>>energy[i]>>pe[i]>>epe[i]>>sigma[i]>>esigma[i];
			Erec[i] = pe[i]/scale;	
			eErec[i] = epe[i]/scale;	
	
			if(i==0){
				init_PE = pe[0];
				init_sigma = sigma[0];
			}
			if(i==0 || i==3) {
				Erec[i] /= 2;
				sigma[i]/=sqrt(2);
                pe[i] /= 2;
			}
			sigma[i] = sigma[i]/pe[i];
			esigma[i] = esigma[i]/pe[i];
			r[i] = Erec[i]/energy[i];
			er[i] = eErec[i]/energy[i];
			cout<<r[i]<<"\t"<<er[i]<<endl;
	}

	TGraphErrors *T1 = new TGraphErrors(n,energy,r,0,er);
	T1->Draw("Ap");
	T1->SetMarkerStyle(20);
    T1->SetMarkerSize(0.8);

	TF1 *f1 = new TF1("f1",get_ratio1,0,7,5);
	f1->SetParameters(1.05412,0.0069017,0.00938879,0.170967,2.95026);
	T1->Fit("f1","","",0,8);
	double *pars = f1->GetParameters();
	ofstream fout1("result1");
	fout1<<pars[0]<<"\t"<<pars[1]<<"\t"<<pars[2]<<"\t"<<pars[3]<<"\t"<<pars[4]<<"\n";
	
	TGraphErrors *T2 = new TGraphErrors(n,energy,sigma,0,esigma);
    T2->Draw("Ap");
    T2->SetMarkerStyle(20);
    T2->SetMarkerSize(1);

    TF1 *f2 = new TF1("f2",get_ratio2,0,8,7);
    f2->SetParameters(0.02803,0.007886,pars[0],pars[1],pars[2],pars[3],pars[4]);

    f2->SetParNames("a","b","P_{0}","P_{1}","P_{2}","P_{3}","P_{4}");
    f2->FixParameter(2,pars[0]);
    f2->FixParameter(3,pars[1]);
    f2->FixParameter(4,pars[2]);
    f2->FixParameter(5,pars[3]);
    f2->FixParameter(6,pars[4]);

	T2->Fit("f2","","",0,8);
	pars = f2->GetParameters();
	ofstream fout2("result2");
	fout2<<pars[0]<<"\t"<<pars[1]<<"\t"<<pars[2]<<"\t"<<pars[3]<<"\t"<<pars[4]<<"\n";

}
