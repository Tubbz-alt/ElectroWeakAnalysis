void runMyClass( bool isBoosted=1 ) {
	TString jetlable="AK8";//AK5
	TString pflable="PF";//PFCHS
	TString boostedlable;
	if(isBoosted)boostedlable="boosted";
	else boostedlable="unboosted";	
	
	gROOT->LoadMacro("MyClass.C+");
	gROOT->LoadMacro("tdrstyle.C");
	setTDRStyle(); //plotting style

	// mkdir plots fold accoding time
	TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
	dir.ReplaceAll("runMyClass.C","");
	dir.ReplaceAll("/./","/");
	TDatime dt;

	char* plot_Dir_DateTime=Form("%sPlots/%i_%i_%i/%s_%s_%s_%i",dir.Data(), dt.GetYear(), dt.GetMonth(), dt.GetDay(), boostedlable.Data(),jetlable.Data(),pflable.Data(), dt.GetTime());
	cout<<plot_Dir_DateTime<<endl;
	gSystem->mkdir(plot_Dir_DateTime, 1);



	TFile *file1;
	if (isBoosted){
		//file1 = new TFile("full_boost_oct3_zmumujetsanalysisntuple.root");
		file1 = new TFile("boosted_full_v4_zmumujetsanalysisntuple.root");
	}else{
		//file1 = new TFile("unboosted_full_zmumujetsanalysisntuple.root");
		//file1 = new TFile("zmumujetsanalysisntuple.root");
		file1 = new TFile("unboosted_full_v1_zmumujetsanalysisntuple.root");
	}

	TFile *fout = new TFile(Form("%s/out_%s.root", plot_Dir_DateTime, boostedlable.Data()), "RECREATE");
	TTree *tree1 = (TTree*) file1->Get("ZJet");


	MyClass* analyzor;
	analyzor =  new MyClass(tree1,jetlable.Data(),pflable.Data(), isBoosted, plot_Dir_DateTime);
	analyzor->Loop();

	fout->Write();
	fout->Close();
}

