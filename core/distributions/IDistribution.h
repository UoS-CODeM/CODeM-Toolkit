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
#ifndef IDISTRIBUTION_H
#define IDISTRIBUTION_H



//#include <sstream>
//#include <string>

// Qt Includes
#include <QVector>
#include <QString>

class AbstractInterpolator;

class IDistribution
{
public:
    IDistribution();
    IDistribution(const IDistribution& dist);
    IDistribution(double sample);
    virtual ~IDistribution();

    virtual IDistribution* clone() const;

    Tigon::DistributionType type() const;

    virtual QVector<double> parameters();

    virtual double sample();
    virtual double mean();
    virtual double variance();
    virtual double median();
    virtual double percentile(double p);
    virtual double std();

    QVector<double>         pdf();
    QVector<double>         cdf();
    virtual QVector<double> pdf(const QVector<double> zVec);
    virtual QVector<double> cdf(const QVector<double> zVec);
    virtual double          pdf(double                   z);
    virtual double          cdf(double                   z);

    virtual void   defineBoundaries(double lb, double ub);
    void           defineResolution(double           dz);
    double          resolution()                   const;
    double          lowerBound()                   const;
    double          upperBound()                   const;
    void           defineZ(QVector<double>            z);
    QVector<double> zSamples();

    virtual void generateZ();
    virtual void generatePDF();
    virtual void calculateCDF();
    void generateEquallySpacedZ();
    void normalise();

    // basic math operations
    virtual void negate();
    virtual void add(double num);
    virtual void add(const IDistribution* other);
    virtual void subtract(double num);
    virtual void subtract(const IDistribution* other);
    virtual void multiply(double num);
    virtual void multiply(const IDistribution* other);
    virtual void divide(double num);
    virtual void divide(const IDistribution* other);
    void reciprocal();


protected:
    Tigon::DistributionType  m_type;
    double                    m_dz;
    double                    m_lb;
    double                    m_ub;
    QVector<double>           m_z;
    QVector<double>           m_pdf;
    QVector<double>           m_cdf;
    int                      m_nSamples;
    AbstractInterpolator*    m_quantileInterpolator;
    AbstractInterpolator*    m_pdfInterpolator;
    AbstractInterpolator*    m_cdfInterpolator;

};

#endif // IDISTRIBUTION_H
