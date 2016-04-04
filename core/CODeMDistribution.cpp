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
#include <core/Distributions/IDistribution.h>
#include <tigon/Utils/NormalisationUtils.h>

namespace CODeM {

CODeMDistribution::CODeMDistribution(IDistribution* d,
                                     const QVector<double> oVec,
                                     double lowerBound,
                                     double upperBound,
                                     const QVector<double> ideal,
                                     const QVector<double> antiIdeal,
                                     double dirPertRad,
                                     double dirPertNorm)
    : m_lb(lowerBound),
      m_ub(upperBound),
      m_pNorm(1)
{
    defineDistribution(d);
    defineIdealAndAntiIdeal(ideal, antiIdeal);
    defineDirection(oVec);
    defineDirectionPertRadius(dirPertRad);
    definePerturbationNorm(dirPertNorm);
}

CODeMDistribution::~CODeMDistribution()
{

}

QVector<double> CODeMDistribution::sampleDistribution()
{
    if(m_distribution.isNull()) {
        return QVector<double>(0);
    }
    double sFactor = m_distribution->sample();

    // scale to the interval [lb ub]
    sFactor = m_lb + sFactor*(m_ub-m_lb);

    // scale the 2-norm direction vector
    QVector<double> samp = m_direction;
    scale(samp,sFactor);

    samp = directionPerturbation(samp, m_directionPertRadius, m_pNorm);

    scaleBackFromUnitBox(samp, m_ideal, m_antiIdeal);

    return samp;
}

void CODeMDistribution::defineDirectionPertRadius(double r)
{
    if(r >= 0.0) {
        m_directionPertRadius = r;
    }
}

void CODeMDistribution::definePerturbationNorm(double p)
{
    if(p > 0.0) {
        m_pNorm = p;
    }
}

void CODeMDistribution::defineDirection(const QVector<double> oVec)
{
    m_direction = oVec;
    normaliseToUnitBox(m_direction, m_ideal, m_antiIdeal);
    toUnitVec(m_direction);
}

void CODeMDistribution::defineIdealAndAntiIdeal(const QVector<double> ideal,
                                                const QVector<double> antiIdeal)
{
    m_ideal = ideal;
    m_antiIdeal = antiIdeal;
}
void CODeMDistribution::defineDistribution(IDistribution* d)
{
    m_distribution = d;
}

} // namespace CODeM
