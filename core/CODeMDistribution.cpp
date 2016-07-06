/****************************************************************************
**
** Copyright (C) 2012-2015 The University of Sheffield (www.sheffield.ac.uk)
**
** This file is part of Liger.
**
** GNU Lesser General Public License Usage
** This file may be used under the terms of the GNU Lesser General
** Public License version 2.1 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU Lesser General Public License version 2.1 requirements
** will be met: http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
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
