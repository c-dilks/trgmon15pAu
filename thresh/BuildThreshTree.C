// takes thresh.dat (table of trigger thresholds, see 
// thresh_columns.info for column information)
// and builds a tree

void BuildThreshTree(const char * datname="thresh.dat")
{
  TTree *th = new TTree("th","th");
  th->ReadFile(datname,"runnum:index:sqrts:JP0:JP1:JP2:LBitSel:LBS1:LBS2:LBS3:SBitSel:SBS1:SBS2:SBS3");
  TFile *tf = new TFile("thresh.root","RECREATE");
  th->Write("th");
};
