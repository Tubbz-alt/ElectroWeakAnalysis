void runMyClass() {
	gROOT->LoadMacro("MyClass.cc+");
	gROOT->LoadMacro("tdrstyle.C");
	setTDRStyle(); //plotting style

	runMyClass_signle(0, "AK5", "PF");
	/*runMyClass_signle(0, "AK5", "PFCHS");
	runMyClass_signle(0, "AK8", "PF");
	runMyClass_signle(0, "AK8", "PFCHS");

	runMyClass_signle(1, "AK5", "PF");
	runMyClass_signle(1, "AK5", "PFCHS");
	runMyClass_signle(1, "AK8", "PF");
	runMyClass_signle(1, "AK8", "PFCHS");*/
}

void runMyClass_signle( bool isBoosted, TString jetlable, TString pflable) {

	TString final_state="Dijets";//
	//TString final_state="ZJet";//

	TString boostedlable;
	if(isBoosted)boostedlable="boosted";
	else boostedlable="unboosted";	
	// mkdir plots fold accoding time
	TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
	dir.ReplaceAll("runMyClass.C","");
	dir.ReplaceAll("/./","/");
	TDatime dt;

	char* plot_Dir_DateTime=Form("%sPlots/%i_%i_%i/%s_%s_%s_%i",dir.Data(), dt.GetYear(), dt.GetMonth(), dt.GetDay(), boostedlable.Data(),jetlable.Data(),pflable.Data(), dt.GetTime());
	cout<<plot_Dir_DateTime<<endl;
	gSystem->mkdir(plot_Dir_DateTime, 1);

	TFile *file1;
	TTree *tree1;

	if (final_state.Contains("ZJet")){
		if (isBoosted){
			file1 = new TFile("boosted_full_v7_zmumujetsanalysisntuple.root");
		}else{
			//file1 = new TFile("unboosted_full_v7_zmumujetsanalysisntuple.root");
			//file1 = new TFile("full_zmumujetsanalysisntuple_highPU.root");
			file1 = new TFile("dijetsanalysisntuple_1000.root");
		}
		tree1 = (TTree*) file1->Get("ZJet");

	}else if (final_state.Contains("Dijets")) {

		file1 = new TFile("dijetsanalysisntuple_1000.root");
		tree1 = (TTree*) file1->Get("Dijets");

	}else{
		cout<<"Wrong final state: "<<final_state<<endl;
	}
	TFile *fout = new TFile(Form("%s/out_%s.root", plot_Dir_DateTime, boostedlable.Data()), "RECREATE");


	MyClass* analyzor;
	cout<<"final_state="<<final_state.Data()<<endl;
	analyzor =  new MyClass(tree1, final_state.Data(), jetlable.Data(), pflable.Data(), isBoosted, plot_Dir_DateTime);
	analyzor->Loop();

	fout->Write();
	fout->Close();

}

