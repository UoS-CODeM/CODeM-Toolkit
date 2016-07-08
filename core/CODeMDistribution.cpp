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
#include <core/CODeMDistribution.h>
#include <core/CODeMOperators.h>
#include <core/RandomDistributions.h>
#include <core/utils/ScalingUtils.h>

using namespace CODeM::Utils;
using std::vector;

namespace CODeM {

CODeMDistribution::CODeMDistribution(IDistribution*        d,
                                     const vector<double>& oVec,
                                     double                lowerBound,
                                     double                upperBound,
                                     const vector<double>& ideal,
                                     const vector<double>& antiIdeal,
                                     double                dirPertRad,
                                     double                dirPertNorm)
    : m_distribution(d),
      m_directionPertRadius(dirPertRad >= 0.0 ? dirPertRad : 0.0),
      m_ideal(ideal),
      m_antiIdeal(antiIdeal),
      m_lb(lowerBound),
      m_ub(upperBound),
      m_pNorm(dirPertNorm > 0.0 ? dirPertNorm : 1.0)
{    
    defineDirection(oVec);
}

CODeMDistribution::~CODeMDistribution()
{
    delete m_distribution;
}

vector<double> CODeMDistribution::sampleDistribution()
{
    if(m_distribution == 0) {
        return vector<double>(0);
    }
    double sFactor = m_distribution->sample();

    // scale to the interval [lb ub]
    sFactor = m_lb + sFactor*(m_ub-m_lb);

    // scale the 2-norm direction vector
    vector<double> samp = m_direction;
    scale(samp, sFactor);

    samp = directionPerturbation(samp, m_directionPertRadius, m_pNorm);

    scaleBackFromUnitBox(samp, m_ideal, m_antiIdeal);

    return samp;
}

void CODeMDistribution::defineDirection(const vector<double> &oVec)
{
    m_direction = oVec;
    normaliseToUnitBox(m_direction, m_ideal, m_antiIdeal);
    toUnitVec(m_direction);
}

} // namespace CODeM
