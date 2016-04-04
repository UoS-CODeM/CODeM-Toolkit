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
#ifndef CODEMDISTRIBUTION_H
#define CODEMDISTRIBUTION_H
#include <tigon/tigon_global.h>
#include <tigon/tigonconstants.h>

// Qt Includes
#include <QVector>

namespace Tigon {
class LinearInterpolator;
}

using namespace Tigon::Representation;

namespace CODeM{


class LIGER_TIGON_EXPORT CODeMDistribution
{
public:
    CODeMDistribution(IDistributionSPtr d,
                      const QVector<qreal> oVec,
                      qreal lowerBound,
                      qreal upperBound,
                      const QVector<qreal> ideal,
                      const QVector<qreal> antiIdeal,
                      qreal dirPertRad,
                      qreal dirPertNorm);
    ~CODeMDistribution();

    QVector<qreal> sampleDistribution();

    void defineDirectionPertRadius(qreal r);
    void definePerturbationNorm(qreal p);
    // 2-norm direction
    void defineDirection(const QVector<qreal> oVec);
    void defineIdealAndAntiIdeal(const QVector<qreal> ideal,
                                 const QVector<qreal> antiIdeal);
    void defineDistribution(IDistributionSPtr d);


private:
    IDistributionSPtr    m_distribution;
    qreal                m_directionPertRadius;
    QVector<qreal>       m_direction;
    QVector<qreal>       m_ideal;
    QVector<qreal>       m_antiIdeal;
    qreal                m_lb;
    qreal                m_ub;
    qreal                m_pNorm;
};

} //namespace CODeM

#endif // CODEMDISTRIBUTION_H
