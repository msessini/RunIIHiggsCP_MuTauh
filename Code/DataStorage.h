#ifndef DataStorage_h
#define DataStorage_h

#include <vector>
#include "TString.h"

class DataStorage {
 public:
  DataStorage();
  virtual ~DataStorage();

  int  GetFile(TString key);
  void StoreFile(TString File, TString savedFile);

  TString timeStamp(bool Date = true, bool Time = true);

 private:
  TString mydir;

  // list of temporary files to be deleted
  std::vector<TString> filesToDelete;

 protected:
  TString inputFileName;

  // Count number of instances of this class.
  // This is used to avoid duplication of local filenames for input and output files.
  static int instanceCounter;
  // Index of this instance of the class
  int instance;

  TString assemblyFileName(unsigned int idx_File);
};
#endif
