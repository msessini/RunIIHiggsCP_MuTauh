#include "Selection_Base.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

Selection_Base::Selection_Base(TString Name_, TString id_):  
  GeV(1000)
  ,Analysis(Name_)
  ,Name(Name_)
  ,systype(id_)
  ,isNtp(false)
  ,runtype(GRID)
  ,doDetails(false)
  ,Lumi(1)
{
  Name+="_";
  Name+=systype;
  Logger(Logger::Info) << "Configuring Selection " << Name << std::endl;
}

Selection_Base::~Selection_Base(){
	Logger(Logger::Info) << "done." << std::endl;
}



void Selection_Base::SetMode(int mode_){  
  mode=Selection_Base::ANALYSIS;
  if(mode_==Selection_Base::RECONSTRUCT) mode=Selection_Base::RECONSTRUCT;
}

void Selection_Base::SetRunType(int runtype_){
  runtype=Selection_Base::Local; 
  if(runtype_==Selection_Base::GRID) runtype=Selection_Base::GRID;
}
