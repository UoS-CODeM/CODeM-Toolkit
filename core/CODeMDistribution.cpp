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
#include <tigon/Representation/Functions/CODeM/CODeMDistribution.h>
#include <tigon/Representation/Functions/CODeM/CODeMOperators.h>
#include <tigon/Representation/Distributions/IDistribution.h>
#include <tigon/Utils/NormalisationUtils.h>

using namespace Tigon;

namespace CODeM {

CODeMDistribution::CODeMDistribution(IDistributionSPtr d,
                                     const QVector<qreal> oVec,
                                     qreal lowerBound,
                                     qreal upperBound,
                                     const QVector<qreal> ideal,
                                     const QVector<qreal> antiIdeal,
                                     qreal dirPertRad,
                                     qreal dirPertNorm)
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

QVector<qreal> CODeMDistribution::sampleDistribution()
{
    if(m_distribution.isNull()) {
        return QVector<qreal>(0);
    }
    qreal sFactor = m_distribution->sample();

    // scale to the interval [lb ub]
    sFactor = m_lb + sFactor*(m_ub-m_lb);

    // scale the 2-norm direction vector
    QVector<qreal> samp = m_direction;
    scale(samp,sFactor);

    samp = directionPerturbation(samp, m_directionPertRadius, m_pNorm);

    scaleBackFromUnitBox(samp, m_ideal, m_antiIdeal);

    return samp;
}

void CODeMDistribution::defineDirectionPertRadius(qreal r)
{
    if(r >= 0.0) {
        m_directionPertRadius = r;
    }
}

void CODeMDistribution::definePerturbationNorm(qreal p)
{
    if(p > 0.0) {
        m_pNorm = p;
    }
}

void CODeMDistribution::defineDirection(const QVector<qreal> oVec)
{
    m_direction = oVec;
    normaliseToUnitBox(m_direction, m_ideal, m_antiIdeal);
    toUnitVec(m_direction);
}

void CODeMDistribution::defineIdealAndAntiIdeal(const QVector<qreal> ideal,
                                                const QVector<qreal> antiIdeal)
{
    m_ideal = ideal;
    m_antiIdeal = antiIdeal;
}
void CODeMDistribution::defineDistribution(IDistributionSPtr d)
{
    m_distribution = d;
}

} // namespace CODeM
