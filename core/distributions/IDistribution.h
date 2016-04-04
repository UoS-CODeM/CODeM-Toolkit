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
#include <tigon/tigon_global.h>
#include <tigon/tigonconstants.h>

//#include <sstream>
//#include <string>

// Qt Includes
#include <QVector>
#include <QString>

namespace Tigon {
class AbstractInterpolator;
}

namespace Tigon {
namespace Representation {

class LIGER_TIGON_EXPORT IDistribution
{
public:
    IDistribution();
    IDistribution(const IDistribution& dist);
    IDistribution(qreal sample);
    virtual ~IDistribution();

    virtual IDistribution* clone() const;

    Tigon::DistributionType type() const;

    virtual QVector<qreal> parameters();

    virtual qreal sample();
    virtual qreal mean();
    virtual qreal variance();
    virtual qreal median();
    virtual qreal percentile(qreal p);
    virtual qreal std();

    QVector<qreal>         pdf();
    QVector<qreal>         cdf();
    virtual QVector<qreal> pdf(const QVector<qreal> zVec);
    virtual QVector<qreal> cdf(const QVector<qreal> zVec);
    virtual qreal          pdf(qreal                   z);
    virtual qreal          cdf(qreal                   z);

    virtual void   defineBoundaries(qreal lb, qreal ub);
    void           defineResolution(qreal           dz);
    qreal          resolution()                   const;
    qreal          lowerBound()                   const;
    qreal          upperBound()                   const;
    void           defineZ(QVector<qreal>            z);
    QVector<qreal> zSamples();

    virtual void generateZ();
    virtual void generatePDF();
    virtual void calculateCDF();
    void generateEquallySpacedZ();
    void normalise();

    // basic math operations
    virtual void negate();
    virtual void add(qreal num);
    virtual void add(const IDistributionSPtr other);
    virtual void subtract(qreal num);
    virtual void subtract(const IDistributionSPtr other);
    virtual void multiply(qreal num);
    virtual void multiply(const IDistributionSPtr other);
    virtual void divide(qreal num);
    virtual void divide(const IDistributionSPtr other);
    void reciprocal();


protected:
    Tigon::DistributionType  m_type;
    qreal                    m_dz;
    qreal                    m_lb;
    qreal                    m_ub;
    QVector<qreal>           m_z;
    QVector<qreal>           m_pdf;
    QVector<qreal>           m_cdf;
    int                      m_nSamples;
    AbstractInterpolator*    m_quantileInterpolator;
    AbstractInterpolator*    m_pdfInterpolator;
    AbstractInterpolator*    m_cdfInterpolator;

};

} // namespace Representation
} // namespace Tigon

Q_DECLARE_METATYPE(Tigon::Representation::IDistributionSPtr)

#endif // IDISTRIBUTION_H
