void mergePUFiles(){
	TFile* in1 = new TFile("Lumi_190456_208686MC_PU_S10_andData.root","READ");
	TFile* in2 = new TFile("HTauTau_MCandData_Pileup.root","READ");
	TFile* out = new TFile("Lumi_OfficialAndHtautau.root","RECREATE");


	// get official histos
	TH1D* off_mc = (TH1D*)in1->Get("MC_Summer12");
	TH1D* off_data = (TH1D*)in1->Get("h_190456_20868");
	TH1D* off_data_p5 = (TH1D*)in1->Get("h_190456_20868_p5");
	TH1D* off_data_m5 = (TH1D*)in1->Get("h_190456_20868_m5");

	// get htautau histos
	TH1D* htt_mc = (TH1D*)in2->Get("mc_pileup");
	TH1D* htt_data = (TH1D*)in2->Get("data_pileup");

	// rename
	off_mc->SetName("official_MC_Summer12");
	off_data->SetName("official_h_190456_20868");
	off_data_p5->SetName("official_h_190456_20868_p5");
	off_data_m5->SetName("official_h_190456_20868_m5");
	htt_mc->SetName("htautau_mc_pileup");
	htt_data->SetName("htautau_data_pileup");

	//out->Write();
	off_mc->Write();
	off_data->Write();
	off_data_p5->Write();
	off_data_m5->Write();
	htt_mc->Write();
	htt_data->Write();

	gROOT->ProcessLine(".q");
}
