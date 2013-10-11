void runMyClass( bool isBoosted=0 ) {
	gROOT->LoadMacro("MyClass.C+");

	gROOT->LoadMacro("tdrstyle.C");
	setTDRStyle(); //plotting style

	// mkdir plots fold accoding time
	TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
	dir.ReplaceAll("runMyClass.C","");
	dir.ReplaceAll("/./","/");
	TDatime dt;

	TString boostedlable;
	if(isBoosted)boostedlable="boosted";
	else boostedlable="unboosted";
	char* plot_Dir_DateTime=Form("%sPlots/%i_%i_%i/%s_%i",dir.Data(), dt.GetYear(), dt.GetMonth(), dt.GetDay(), boostedlable.Data(), dt.GetTime());
	cout<<plot_Dir_DateTime<<endl;
	gSystem->mkdir(plot_Dir_DateTime, 1);



	TFile *file1;
	if (isBoosted){
		file1 = new TFile("full_boost_oct3_zmumujetsanalysisntuple.root");
	}else{
		file1 = new TFile("unboosted_full_zmumujetsanalysisntuple.root");
	}

	TFile *fout = new TFile(Form("%s/out_%s.root", plot_Dir_DateTime, boostedlable.Data()), "RECREATE");
	TTree *tree1 = (TTree*) file1->Get("ZJet");

	MyClass* analyzor;
	if (isBoosted){
		cout<<"boosted Z analysis"<<endl;
		analyzor =  new MyClass(tree1,"AK8","PF", isBoosted, plot_Dir_DateTime);
	}else{
		cout<<"unboosted Z analysis"<<endl;
		analyzor =  new MyClass(tree1,"AK5","PF", isBoosted, plot_Dir_DateTime);
	}
	analyzor->Loop();

	fout->Write();
	fout->Close();
}

