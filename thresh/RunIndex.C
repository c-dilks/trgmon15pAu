// takes runlog.root, generated in trgmon,
// and writes a text file with columns [run#] [index#] [sqrt(s)] [download?]
// (the index is a running integer corresponding to 
//  run number used in trgmon code)
//
//  download boolean is true if no file in thfiles/ corresponds
//  to the run number and thus needs to be downloaded


void RunIndex(const char * rootname="../200/runlog.root", Int_t sqrts=200)
{
  // list existing downloaded information
  system("ls thfiles | sed 's/th_//g' > tmplist");
  TTree *synctree = new TTree("synctree","synctree");
  synctree->ReadFile("tmplist","rn");
  char synctree_cut[32];
  Bool_t download_this;

  TFile *infile = new TFile(rootname,"READ");
  TTree *runtr = (TTree*)infile->Get("tr");

  Int_t runnum,index;
  runtr->SetBranchAddress("runnum",&runnum);
  runtr->SetBranchAddress("index",&index);

  system("touch index.dat");
  system("rm index.dat");
  system("touch index.dat");
  for(Int_t i=0; i<runtr->GetEntries(); i++)
  {
    runtr->GetEntry(i);
    sprintf(synctree_cut,"rn==%d",runnum);
    if(synctree->GetEntries(synctree_cut)) download_this=0;
    else download_this=1;
    gSystem->RedirectOutput("index.dat");
    printf("%d %d %d %d\n",runnum,index,sqrts,download_this);
    gSystem->RedirectOutput(0);
  };
};
