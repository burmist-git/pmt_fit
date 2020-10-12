Int_t pmtTypicalGainR762_NA0138(){

  const Int_t n = 5;
  Double_t hv[n];
  Double_t hvErr[n];
  Double_t gain[n];
  Double_t gainErr[n];
  Double_t hvErr[n];
  Double_t gainErr[n];
  
  const Int_t nn = 100;
  Double_t hv_s_min = 300;
  Double_t hv_s_max = 1100;
  Double_t hv_s[nn];
  Double_t hv_tt[nn];
  Double_t gain_s[nn];
  Double_t gain_ss[nn];
  Double_t deltaGain[nn];
  const Int_t npar = 2;
  Double_t f2paramsFit[npar];
  f2paramsFit[0] =6.90259e-01;//k
  f2paramsFit[1] =1.92359e-01;//a
  Int_t i = 0;
  for(i = 0;i<nn;i++){
    hv_s[i] = hv_s_min + (hv_s_max - hv_s_min)*i/(nn-1);
    gain_s[i] = 2.45e-16*TMath::Power(hv_s[i],7.43);
    gain_ss[i] = fun2(&hv_s[i], f2paramsFit);
    deltaGain[i] = (gain_s[i] - gain_ss[i])/gain_ss[i];
  }

  hv[0] = 500;
  hv[1] = 600;
  hv[2] = 700;
  hv[3] = 800;
  hv[4] = 900;

  for(i = 0;i<n;i++)
    hvErr[i] = hv[i]*0.001/2;

  gain[0] = 1.8E+04; //  500 V
  gain[1] = 7.0E+04; //  600 V
  gain[2] = 2.1E+05; //  700 V
  gain[3] = 5.4E+05; //  800 V
  gain[4] = 1.0E+06; // 1000 V

  for(i = 0;i<n;i++)
    gainErr[i] = gain[i]/10.0/10;

  TGraphErrors *gr1 = new TGraphErrors(n,hv,gain,hvErr,gainErr);
  //TGraphErrors *gr1 = new TGraphErrors(n,hv,gain,0,gainErr);
  //TGraphErrors *gr1 = new TGraphErrors(n,hv,gain,0,0);
  gr1->SetMarkerStyle(21);
  gr1->SetMarkerColor(kBlack);
  gr1->SetTitle(""); 
  gr1->SetMinimum(0.0);

  TGraphErrors *gr1_s = new TGraphErrors(nn,hv_s,gain_s,0,0);
  gr1_s->SetMarkerStyle(21);
  gr1_s->SetMarkerColor(kRed);
  gr1_s->SetTitle(""); 
  gr1_s->SetMinimum(0.0);
  gr1_s->SetLineWidth(2);
  gr1_s->SetLineColor(kRed);

  TGraphErrors *gr1_ss = new TGraphErrors(nn,hv_s,gain_ss,0,0);
  gr1_ss->SetMarkerStyle(21);
  gr1_ss->SetMarkerColor(kRed);
  gr1_ss->SetTitle(""); 
  gr1_ss->SetMinimum(0.0);
  gr1_ss->SetLineWidth(2);
  gr1_ss->SetLineColor(kRed);

  TGraphErrors *gr1_deltaGain = new TGraphErrors(nn,hv_s,deltaGain,0,0);
  gr1_deltaGain->SetMarkerStyle(21);
  gr1_deltaGain->SetMarkerColor(kBlack);
  gr1_deltaGain->SetTitle(""); 
  gr1_deltaGain->SetMinimum(0.0);
  gr1_deltaGain->SetLineWidth(2);
  gr1_deltaGain->SetLineColor(kBlack);

  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  
  TCanvas *c1 = new TCanvas("c1","plots",10,10,800,600);
  //TCanvas *c1 = new TCanvas("c1","plots",10,10,800,800*3/4.0);
  c1->Clear();
  c1->SetFillColor(kWhite);

  //////////////////// FIT ///////////////
  ////const Int_t npar = 2;
  ////Double_t f2params[npar];
  ////f2params[0] = 0.73;//k
  ////f2params[1] = 0.3;//a
  //TF1 *gainF = new TF1("gainF",fun2,500,1200,npar);  
  //gr1->Fit("gainF");  
  //gr1->Draw("AP");
  //gr1->GetXaxis()->SetTitle("Absolute value of high voltage on PMT, V");
  //gr1->GetYaxis()->SetTitle("PMT gain");
  ////////////////////////////////////////

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr1,"p");
  //mg->Add(gr1_s,"l");
  mg->Add(gr1_ss,"l");
  //mg->Add(gr1_deltaGain,"l");
  //mg->SetMaximum(1500000.0);
  mg->SetMinimum(1000.0);
  mg->Draw("AP");

  mg->GetXaxis()->SetTitle("Absolute value of high voltage on PMT, V");
  mg->GetYaxis()->SetTitle("PMT gain");

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(gr1, "Measurements","p");
  leg->AddEntry(gr1_s, "Fit ","p");
  //leg->Draw();
  //mg->GetXaxis()->SetTitle("Number of incident electrons");
  //mg->GetYaxis()->SetTitle("Number of p.e.");

  hv_tt[0] = 400;
  hv_tt[1] = 450;
  hv_tt[2] = 500;
  hv_tt[3] = 550;
  hv_tt[4] = 600;
  hv_tt[5] = 650;
  hv_tt[6] = 700;
  hv_tt[7] = 750;
  hv_tt[8] = 800;
  hv_tt[9] = 850;
  hv_tt[10] = 900;
  hv_tt[11] = 950;
  hv_tt[12] =1000;
  hv_tt[13] =1050;
  hv_tt[14] =1100;
  cout<<"a = "<<f2paramsFit[1]<<endl;
  cout<<"k = "<<f2paramsFit[0]<<endl;
  for(i = 0;i<15;i++){
    if(hv_tt[i]<500)
      cout<<hv_tt[i]<<" V  "<<fun2(&hv_tt[i], f2paramsFit)/1E3<<" * 10^3"<<endl;
    if(hv_tt[i]>450 && hv_tt[i]<650)
      cout<<hv_tt[i]<<" V  "<<fun2(&hv_tt[i], f2paramsFit)/1E4<<" * 10^4"<<endl;
    else if(hv_tt[i]>600 && hv_tt[i]<900)
      cout<<hv_tt[i]<<" V  "<<fun2(&hv_tt[i], f2paramsFit)/1E5<<" * 10^5"<<endl;
    else if(hv_tt[i]>800)
      cout<<hv_tt[i]<<" V  "<<fun2(&hv_tt[i], f2paramsFit)/1E6<<" * 10^6"<<endl;
  }
  
  return 0;
}

Double_t fun2(Double_t *x, Double_t *par) {
  //Double_t *p1 = &par[0];
  // Double_t *p2 = &par[1];
  //cout<<x[0]<<endl;
  Double_t p3 = 10;//n
  Double_t result = TMath::Power(par[1],p3)/TMath::Power((p3+1),par[0]*p3)*TMath::Power(x[0],p3*par[0]);
  return result;
}
