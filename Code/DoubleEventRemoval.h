#ifndef DoubleEventRemoval_h
#define DoubleEventRemoval_h

#include <set>

class DoubleEventRemoval{

 public:
  DoubleEventRemoval();

  bool    CheckDoubleEvents(int run, int event);

  ~DoubleEventRemoval();

 private:
  static std::set<std::pair<int, int> > fRunEventPair; 



};
#endif
