/****************************************************************************
**
** The MIT License (MIT)
**
** Copyright (c) 2016 The University of Sheffield (www.sheffield.ac.uk)
**
** Permission is hereby granted, free of charge, to any person obtaining a copy
** of this software and associated documentation files (the "Software"), to deal
** in the Software without restriction, including without limitation the rights
** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
** copies of the Software, and to permit persons to whom the Software is
** furnished to do so, subject to the following conditions:
**
** The above copyright notice and this permission notice shall be included in all
** copies or substantial portions of the Software.
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
** OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
** SOFTWARE
**
****************************************************************************/
#ifndef CODEMDISTRIBUTION_H
#define CODEMDISTRIBUTION_H

#include <core/CODeMGlobal.h>
#include <vector>


namespace CODeM{
class IDistribution;
//namespace Utils {
//class LinearInterpolator;
//}

class CODeMDistribution
{
public:
    CODeMDistribution(IDistribution*             d,
                      const std::vector<double>& oVec,
                      double                     lowerBound,
                      double                     upperBound,
                      const std::vector<double>& ideal,
                      const std::vector<double>& antiIdeal,
                      double                     dirPertRad,
                      double                     dirPertNorm);
    ~CODeMDistribution();

    std::vector<double> sampleDistribution();

private:
    // 2-norm direction
    void defineDirection(const std::vector<double> &oVec);

    IDistribution*       m_distribution;
    double               m_directionPertRadius;
    std::vector<double>  m_direction;
    std::vector<double>  m_ideal;
    std::vector<double>  m_antiIdeal;
    double               m_lb;
    double               m_ub;
    double               m_pNorm;
};

} //namespace CODeM

#endif // CODEMDISTRIBUTION_H
