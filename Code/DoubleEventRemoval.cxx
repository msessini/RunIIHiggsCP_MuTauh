#define GoodRunList_cxx

#include "DoubleEventRemoval.h"

std::set<std::pair<int, int> > DoubleEventRemoval::fRunEventPair;

DoubleEventRemoval::DoubleEventRemoval()
{
}

DoubleEventRemoval::~DoubleEventRemoval()
{
}


bool DoubleEventRemoval::CheckDoubleEvents(int run, int event) {
  //static std::set<std::pair<int, int> > fRunEventPair;
   //check whether pair is already in set and processed, otherwise insert it in set
  return fRunEventPair.insert(std::make_pair(run,event)).second;
} 
