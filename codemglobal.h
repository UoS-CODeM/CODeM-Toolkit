#ifndef CODEMGLOBAL
#define CODEMGLOBAL

#include <cmath>
#include <random>

namespace CODeM {

double randUni() { return (1.0 * std::rand() / RAND_MAX); }

const double PI(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899);

/// Distribution constants
const int   DistMinNSamples(3);
const int   DistNSamples(500);
const qreal DistMinInterval(0.001);
const qreal DistPeakMinN(1.0);
const qreal DistPeakMaxN(50.0);
const qreal DistPeakMinNBasisFunc(30.0);
const qreal DistPeakMaxNBasisFunc(150.0);

} // namespace CODeM

#endif // CODEMGLOBAL

