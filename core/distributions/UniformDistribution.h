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
#ifndef UNIFORMDISTRIBUTION_H
#define UNIFORMDISTRIBUTION_H
#include <tigon/tigon_global.h>
#include <tigon/Representation/Distributions/IDistribution.h>

// boost does not allow to use forward declerations
#include <boost/math/distributions/uniform.hpp>

namespace Tigon {
namespace Representation {

class LIGER_TIGON_EXPORT UniformDistribution : public IDistribution
{
public:
    UniformDistribution();
    UniformDistribution(const UniformDistribution& dist);
    UniformDistribution(qreal lb, qreal ub);
    UniformDistribution(QVector<qreal> parameters);
    virtual ~UniformDistribution();

    UniformDistribution* clone() const;

    void defineBoundaries(qreal lb, qreal ub);

    qreal sample();
    qreal mean();
    qreal median();
    qreal percentile(qreal p);
    qreal variance();
    qreal std();

    //to make overrides visible to the compiler
    using IDistribution::pdf;
    using IDistribution::cdf;
    qreal pdf(qreal       z);
    qreal cdf(qreal       z);

    void  generateZ();
    void  generatePDF();

    QVector<qreal> parameters();

private:
    boost::math::uniform_distribution<qreal>* m_uniDist;
};

} // namespace Representation
} // namespace Tigon

#endif // UNIFORMDISTRIBUTION_H
