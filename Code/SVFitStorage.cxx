#include "SVFitStorage.h"
#include "Parameters.h"
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "SVFitObject.h"
#include <iostream>
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "TChainIndex.h"
#include "TTreeIndex.h"

SVFitStorage::SVFitStorage():
	outfile_(0),
	outtree_(0),
	intree_(0),
	index_(0),
	treeName_("invalid"),
	suffix_(""),
	svfit_(0),
	b_RunNumber_(0),
	b_LumiNumber_(0),
	b_EventNumber_(0),
	b_svfit_(0),
	isConfigured_(false),
	intreeLoaded_(false){

	TString thelib= getenv ("DATAFORMATS_LIB");
	gSystem->Load(thelib.Data());

	inputFileName = "SVFitInput_temp";

	// Configure(...) MUST be called before instance of this class can be called
}

SVFitStorage::~SVFitStorage(){
	if (isConfigured_){
		SaveTree();

		// check state of output information
		bool storeOutFile = (outtree_->GetEntries() > 0);
		TString outfileName = outfile_->GetName();

		// cleaning up:
		// make sure to delete the TTree objects before closing/destroying the TFile they are associated to
		delete outtree_;
		delete intree_;
		delete outfile_;

		delete svfit_;

		//Store file on the grid
		if (storeOutFile){ // do not copy empty file to dCache
			StoreFile(outfileName , storageFileName_);
			Logger(Logger::Info) << outfileName.Data() << " saved to the grid " << storageFileName_.Data() << std::endl;
		}
	}
	Logger(Logger::Debug) << "Properly destroyed." << std::endl;
}

void SVFitStorage::Configure(TString datasetName, TString suffix /* ="" */){
	if (isConfigured_){
		Logger(Logger::Warning) << "SVFitStorage object has been configured before. Make sure to configure only once. Abort config..." << std::endl;
		return;
	}

	suffix_ = suffix;
	treeName_ = "Tree_" +  datasetName;

	if (suffix_ != ""){
		treeName_ = treeName_ + "_" + suffix_;
		inputFileName = inputFileName + suffix_;
	}

	// Specify file name of output file
	Parameters Par; // assumes configured in Analysis.cxx
	TString key = "OutputFileSVFit" + suffix_ + ":";
	Par.GetString(key, storageFileName_);
	TString outputFileLocal = "MySVFIT" + suffix_ + TString::Itoa(instance,10) + ".root";
	// Load output file
	outfile_ = TFile::Open(outputFileLocal, "RECREATE");
	if (!outfile_) {
		Logger(Logger::Error) << outputFileLocal << " could not be created" << std::endl;
		return;
	}

	// create input and output tree
	// make sure to set name and title to proper values below
	outtree_= new TTree("temp", "temp");
	intree_ = new TChain("temp2", "temp2");
	intree_->SetDirectory(0);

	// allocate memory for objects to be stored in tree
	svfit_ = new SVFitObject();

	// setup output tree
	outtree_->SetName(treeName_);
	outtree_->SetTitle(treeName_);
	outtree_->Branch("RunNumber", &RunNumber_);
	outtree_->Branch("LumiNumber", &LumiNumber_);
	outtree_->Branch("EventNumber", &EventNumber_);
	outtree_->Branch("svfit", &svfit_);

	isConfigured_ = true;

	// setup input tree
	LoadTree();
}

void SVFitStorage::LoadTree(){
	if ( !isConfigured_ ){
		Logger(Logger::Error) << "SVFitStorage must be configured before LoadTree can be called." << std::endl;
		return;
	}

	TString key = "InputFileSVFit" + suffix_ + ":";
	int nfiles = GetFile(key);
	if (nfiles == 0) {
		Logger(Logger::Warning) << "Key not found: " << key <<
				"\n\tNo input SVFit storage specified. Will calculate all SVFit values while running." << std::endl;
	} else {
		intree_->SetName(treeName_);
		intree_->SetTitle(treeName_);
		TDirectory *gdirectory_save = gDirectory;
		int nFilesLoaded = 0;
		for (int i = 0; i < nfiles; i++) {
			TString name = assemblyFileName(i);
			if ( isTreeInFile(name) ){ // check if tree exists in file (avoids TChain error)
				intree_->Add(name);
				nFilesLoaded++;
			}
		}
		if( nFilesLoaded == 0) {
			Logger(Logger::Info) << "None of the input files contains a tree " << treeName_
					<< "\n\t All SVFit values need to be calculated." << std::endl;
			intreeLoaded_ = false;
			gDirectory = gdirectory_save;
			gDirectory->cd();
			return;
		}
		if (intree_->LoadTree(0) < 0)
			Logger(Logger::Error) << "Input TChain was not loaded correctly." << std::endl;

		//Set branches
		intree_->SetBranchAddress("RunNumber", &RunNumber_, &b_RunNumber_);
		intree_->SetBranchAddress("LumiNumber", &LumiNumber_, &b_LumiNumber_);
		intree_->SetBranchAddress("EventNumber", &EventNumber_, &b_EventNumber_);
		intree_->SetBranchAddress("svfit", &svfit_, &b_svfit_);

		// Create and set index
		Logger(Logger::Debug) << "Building the index ..." << std::endl;
		index_ = new TTreeIndex(intree_,"(RunNumber<<13) + LumiNumber", "EventNumber");
		intree_->SetTreeIndex(index_);

		gDirectory = gdirectory_save;
		gDirectory->cd();

		Logger(Logger::Verbose) << "Input TTree " << treeName_ << " has been loaded." << std::endl;
		intreeLoaded_ = true;
	}
}

bool SVFitStorage::isTreeInFile(TString fileName){
	TFile* file = new TFile(fileName.Data(),"READ");
	if (!file || file->IsZombie()){
		Logger(Logger::Error) << "File " << fileName << " does not exist or is corrupted." << std::endl;
		delete file;
		return false;
	}
	bool hasTree = true;
	TObject* obj = file->Get(treeName_);
	if (!obj || !obj->InheritsFrom(TTree::Class())){
		Logger(Logger::Verbose) << "File " << fileName << " does not contain a tree named " << treeName_ << " -> Skip"<< std::endl;
		hasTree = false;
	}

	delete file;
	return hasTree;
}


void SVFitStorage::SaveTree(){
	if ( !isConfigured_ ){
		Logger(Logger::Error) << "SVFitStorage must be configured before SaveTree can be called." << std::endl;
		return;
	}
	if ( outtree_->GetEntries() < 1){
		Logger(Logger::Info) << "Output tree " << treeName_ << " contains no events, thus no output file is created." << std::endl;
		return;
	}

	//Save output
	TDirectory *gdirectory_save = gDirectory;
	outfile_->cd();
	outtree_->Write(treeName_);
	// do not close outfile_ here, because outtree_ is associated to it
	gDirectory = gdirectory_save;
	gDirectory->cd();
	Logger(Logger::Info) << "SVFit_Tree saved to " << outfile_->GetName() << std::endl;
}

void SVFitStorage::SaveEvent(Int_t RunNumber, Int_t LumiNumber, Int_t EventNumber, SVFitObject* svfit){
	if (!isConfigured_) {
		Logger(Logger::Error) << "SVFitStorage must be configured before SaveTree can be called." << std::endl;
		return;
	}

	//Fill event
	RunNumber_ = RunNumber;
	LumiNumber_ = LumiNumber;
	EventNumber_ = EventNumber;
	svfit_ = svfit;
	outtree_->Fill();
}

SVFitObject* SVFitStorage::GetEvent(UInt_t RunNumber, UInt_t LumiNumber, UInt_t EventNumber){
	if (!isConfigured_) {
		Logger(Logger::Error) << "SVFitStorage must be configured before GetEvent can be called." << std::endl;
		*svfit_ = SVFitObject(); // invalid object
		return svfit_;
	}
	if (!intreeLoaded_) {
		Logger(Logger::Verbose) << "No input tree loaded, thus GetEvent does not work." << std::endl;
		*svfit_ = SVFitObject(); // invalid object
		return svfit_;
	}
	if ( RunNumber != 1) { // true for data, false for MC
		if (RunNumber > 262143 || LumiNumber > 4095) {
			Logger(Logger::Error) << "In Data RunNumber must be smaller than 262144 and LumiNumber smaller than 4096" <<
					"\n\tTrying to access run " << RunNumber << " and lumi " << LumiNumber <<
					"\n\tExpect undefined behavior!" << std::endl;
			*svfit_ = SVFitObject(); // invalid object
			return svfit_;
		}
	}

	//Get tree entry using index and then get svfit
	Logger(Logger::Debug) << "Try to access run " << RunNumber << ", lumi " << LumiNumber << ", Event " << EventNumber <<
			"\n\ti.e. major = " << ((RunNumber << 13) + LumiNumber) << ", minor = " << EventNumber << std::endl;
	Long64_t local = index_->GetEntryNumberWithIndex((RunNumber << 13) + LumiNumber, EventNumber);

	// make sure that entry exists, otherwise return invalid result
	if (local < 0) {
		Logger(Logger::Debug) << "Index " << local << " not available in input tree." << std::endl;
		*svfit_ = SVFitObject(); // invalid object
		return svfit_;
	}
	else {
		Logger(Logger::Debug) << "Accessing index " << local << std::endl;
		intree_->GetEntry(local);
	}

	// check if correct event was loaded
	if (RunNumber_ != RunNumber) {
		Logger(Logger::Error) << "Event " << EventNumber << " was not loaded correctly."<< std::endl;
		*svfit_ = SVFitObject(); // rather return invalid object than wrong one
	}
	return svfit_;
}
