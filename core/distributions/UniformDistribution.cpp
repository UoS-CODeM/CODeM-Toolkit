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
#include <tigon/Representation/Distributions/UniformDistribution.h>
#include <tigon/Random/RandomGenerator.h>

namespace Tigon {
namespace Representation {

UniformDistribution::UniformDistribution()
{
    m_uniDist = 0;
    m_type = Tigon::UniformDistType;
    defineBoundaries(0.0, 1.0);
    defineResolution(m_ub-m_lb);
}

UniformDistribution::UniformDistribution(const UniformDistribution& dist)
    : IDistribution(dist)
{
    m_type = Tigon::UniformDistType;
    m_uniDist = new boost::math::uniform_distribution<qreal>(m_lb, m_ub);
}

UniformDistribution::UniformDistribution(qreal lb, qreal ub)
{
    m_uniDist = 0;
    m_type = Tigon::UniformDistType;
    defineBoundaries(lb, ub);
    defineResolution(m_ub-m_lb);
}

UniformDistribution::UniformDistribution(QVector<qreal> parameters)
{
    m_uniDist = 0;
    m_type = Tigon::UniformDistType;
    qreal lb = 0.0;
    qreal ub = 1.0;
    if(parameters.size() > 0) {
        lb = parameters[0];
        if((parameters.size() > 1) && (parameters[1] > lb)) {
            ub = parameters[1];
        } else {
            ub = lb + Tigon::DistMinInterval;
        }
    }
    defineBoundaries(lb, ub);
    defineResolution((m_ub-m_lb)/(Tigon::DistNSamples-1));
}

UniformDistribution::~UniformDistribution()
{
    delete m_uniDist;
    m_uniDist = 0;
}

UniformDistribution* UniformDistribution::clone() const
{
    return (new UniformDistribution(*this));
}

void UniformDistribution::defineBoundaries(qreal lb, qreal ub)
{
    if(lb >= ub) {
        if(lb == 0) {
            ub = Tigon::DistMinInterval;
        } else if(lb > 0) {
            ub = lb * (1 + Tigon::DistMinInterval);
        } else {
            lb = ub * (1 + Tigon::DistMinInterval);
        }
    }

    IDistribution::defineBoundaries(lb, ub);

    if(m_uniDist != 0) {
        delete m_uniDist;
        m_uniDist = 0;
    }

    m_uniDist = new boost::math::uniform_distribution<qreal>(m_lb, m_ub);
}


qreal UniformDistribution::sample()
{
    return TRAND.randUni(m_ub - m_lb, m_lb);
}

qreal UniformDistribution::mean()
{
    return (m_ub + m_lb) / 2.0;
}

qreal UniformDistribution::median()
{
    return boost::math::median(*m_uniDist);
}

qreal UniformDistribution::percentile(qreal p)
{
    if(p >= 1.0) {
        return m_ub;
    } else if(p <= 0.0) {
        return m_lb;
    }
    return boost::math::quantile(*m_uniDist, p);
}

qreal UniformDistribution::variance()
{
    return boost::math::variance(*m_uniDist);
}

qreal UniformDistribution::std()
{
    return boost::math::standard_deviation(*m_uniDist);
}

qreal UniformDistribution::pdf(qreal z)
{
    return boost::math::pdf(*m_uniDist, z);
}

qreal UniformDistribution::cdf(qreal z)
{
    return boost::math::cdf(*m_uniDist, z);
}

void UniformDistribution::generateZ()
{
    generateEquallySpacedZ();
}

void UniformDistribution::generatePDF()
{
    if(m_z.isEmpty()) {
        generateZ();
    }

    qreal probability = 1.0/(m_ub - m_lb);
    m_pdf = QVector<qreal>(m_nSamples, probability);
}

QVector<qreal> UniformDistribution::parameters()
{
    QVector<qreal> params;
    params << lowerBound() << upperBound();
    return params;
}

} // namespace Representation
} // namespace Tigon
