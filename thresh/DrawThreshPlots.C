void DrawThreshPlots(const char * fname="thresh.root")
{
  TFile *infile = new TFile(fname,"READ");
  TTree *th = (TTree*)infile->Get("th");
  gStyle->SetOptStat(0);
  Float_t imax = th->GetMaximum("index");

  TFile * outfile = new TFile("th_plots.root","RECREATE");

  // threshold list
  //////////////////////////////////////////////////
  const Int_t N_THRESH = 11;
  enum thresh_enum {kJP0,kJP1,kJP2,
                    kLBitSel,kLBS1,kLBS2,kLBS3,
                    kSBitSel,kSBS1,kSBS2,kSBS3};
  thresh_enum thresh;
  
  char thresh_str[N_THRESH][32];
  strcpy(thresh_str[kJP0],"FMS-JP0");
  strcpy(thresh_str[kJP1],"FMS-JP1");
  strcpy(thresh_str[kJP2],"FMS-JP2");
  strcpy(thresh_str[kLBitSel],"FMSlarge-BitSelect");
  strcpy(thresh_str[kLBS1],"FMSlarge-BS-th1");
  strcpy(thresh_str[kLBS2],"FMSlarge-BS-th2");
  strcpy(thresh_str[kLBS3],"FMSlarge-BS-th3");
  strcpy(thresh_str[kSBitSel],"FMSsmall-BitSelect");
  strcpy(thresh_str[kSBS1],"FMSsmall-BS-th1");
  strcpy(thresh_str[kSBS2],"FMSsmall-BS-th2");
  strcpy(thresh_str[kSBS3],"FMSsmall-BS-th3");

  char thresh_br[N_THRESH][32];
  strcpy(thresh_br[kJP0],"JP0");
  strcpy(thresh_br[kJP1],"JP1");
  strcpy(thresh_br[kJP2],"JP2");
  strcpy(thresh_br[kLBitSel],"LBitSel");
  strcpy(thresh_br[kLBS1],"LBS1");
  strcpy(thresh_br[kLBS2],"LBS2");
  strcpy(thresh_br[kLBS3],"LBS3");
  strcpy(thresh_br[kSBitSel],"SBitSel");
  strcpy(thresh_br[kSBS1],"SBS1");
  strcpy(thresh_br[kSBS2],"SBS2");
  strcpy(thresh_br[kSBS3],"SBS3");
  //////////////////////////////////////////////////


  // define graphs and set branch addresses
  TGraph * thresh_gr[11];
  Float_t fthresh[N_THRESH];
  for(Int_t i=0; i<N_THRESH; i++)
  {
    thresh_gr[i] = new TGraph();
    thresh_gr[i]->SetTitle(thresh_str[i]);
    thresh_gr[i]->SetName(thresh_br[i]);
    thresh_gr[i]->SetLineColor(kRed);
    thresh_gr[i]->SetLineWidth(2);
    th->SetBranchAddress(thresh_br[i],&(fthresh[i]));
  };
  Float_t findex;
  th->SetBranchAddress("index",&findex);


  // tree loop
  for(Int_t i=0; i<th->GetEntries(); i++)
  {
    th->GetEntry(i);
    for(Int_t j=0; j<N_THRESH; j++) 
      thresh_gr[j]->SetPoint(i,findex,fthresh[j]);
  };
    
  // draw pdf
  TCanvas *cc = new TCanvas("cc","cc",700,500);
  cc->SetGrid();
  char pdfname[32];
  for(Int_t i=0; i<N_THRESH; i++)
  {
    if(i==0) strcpy(pdfname,"th_plots.pdf(");
    else if (i+1==N_THRESH) strcpy(pdfname,"th_plots.pdf)");
    else strcpy(pdfname,"th_plots.pdf");
    thresh_gr[i]->Draw("AL");
    cc->Print(pdfname,"pdf");
    thresh_gr[i]->Write();
  };
};
