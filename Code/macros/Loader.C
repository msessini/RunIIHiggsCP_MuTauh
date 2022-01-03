#include<vector>
#ifdef __CINT__
#pragma link C++ class vector<vector<int> >;
#pragma link C++ class vector<vector<float> >;
#pragma link C++ class vector<vector<vector<float> > >;
#pragma link C++ class vector<vector<unsigned int> >;
#pragma link C++ class vector<vector<vector<int> > >;
#else
template class std::vector<std::vector<int> >;
template class std::vector<std::vector<float> >;
template class std::vector<std::vector<std::vector<float> > >;
template class std::vector<std::vector<unsigned int> >;
template class std::vector<std::vector<std::vector< int> > >;
#endif
