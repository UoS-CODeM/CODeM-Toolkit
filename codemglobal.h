#ifndef CODEMGLOBAL
#define CODEMGLOBAL

//#include <cmath>
#include <random>

namespace CODeM {

double randUni() { return (double)(std::rand()) / RAND_MAX; }

const double PI(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899);

/// Distribution constants
const int    DistMinNSamples(3);
const int    DistNSamples(500);
const double DistMinInterval(0.001);
const double DistPeakMinN(1.0);
const double DistPeakMaxN(50.0);
const double DistPeakMinNBasisFunc(30.0);
const double DistPeakMaxNBasisFunc(150.0);

} // namespace CODeM

#endif // CODEMGLOBAL

