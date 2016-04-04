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

// Qt Includes
#include <QVector>

class LinearInterpolator;

namespace CODeM{

class CODeMDistribution
{
public:
    CODeMDistribution(IDistribution* d,
                      const QVector<double> oVec,
                      double lowerBound,
                      double upperBound,
                      const QVector<double> ideal,
                      const QVector<double> antiIdeal,
                      double dirPertRad,
                      double dirPertNorm);
    ~CODeMDistribution();

    QVector<double> sampleDistribution();

    void defineDirectionPertRadius(double r);
    void definePerturbationNorm(double p);
    // 2-norm direction
    void defineDirection(const QVector<double> oVec);
    void defineIdealAndAntiIdeal(const QVector<double> ideal,
                                 const QVector<double> antiIdeal);
    void defineDistribution(IDistribution* d);


private:
    IDistribution*    m_distribution;
    double                m_directionPertRadius;
    QVector<double>       m_direction;
    QVector<double>       m_ideal;
    QVector<double>       m_antiIdeal;
    double                m_lb;
    double                m_ub;
    double                m_pNorm;
};

} //namespace CODeM

#endif // CODEMDISTRIBUTION_H
