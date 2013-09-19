void runMyClass() {
	gROOT->LoadMacro("MyClass.C");
	//TFile *file1 =new TFile("zmumujetsanalysisntuple.root");
	TFile *file1 =new TFile("full_zmumujetsanalysisntuple.root");


	TFile *fout = new TFile("out.root", "RECREATE");

	TTree *tree1 = (TTree*) file1->Get("ZJet");
	MyClass analyzor(tree1);
	analyzor.Loop();

	fout->Write();
	fout->Close();
}

