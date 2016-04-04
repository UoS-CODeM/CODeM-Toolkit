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

#include <core/Distributions/IDistribution.h>

// boost does not allow to use forward declerations
#include <boost/math/distributions/uniform.hpp>

class UniformDistribution : public IDistribution
{
public:
    UniformDistribution();
    UniformDistribution(const UniformDistribution& dist);
    UniformDistribution(double lb, double ub);
    UniformDistribution(QVector<double> parameters);
    virtual ~UniformDistribution();

    UniformDistribution* clone() const;

    void defineBoundaries(double lb, double ub);

    double sample();
    double mean();
    double median();
    double percentile(double p);
    double variance();
    double std();

    //to make overrides visible to the compiler
    using IDistribution::pdf;
    using IDistribution::cdf;
    double pdf(double       z);
    double cdf(double       z);

    void  generateZ();
    void  generatePDF();

    QVector<double> parameters();

private:
    boost::math::uniform_distribution<double>* m_uniDist;
};

#endif // UNIFORMDISTRIBUTION_H
