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
#ifndef SAMPLEDDISTRIBUTION_H
#define SAMPLEDDISTRIBUTION_H
#include <tigon/tigon_global.h>
#include <tigon/tigonconstants.h>
#include <tigon/Representation/Distributions/IDistribution.h>

// Qt Includes
#include <QVector>

namespace Tigon {
namespace Representation {

class LIGER_TIGON_EXPORT SampledDistribution : public IDistribution
{
public:
    SampledDistribution();
    SampledDistribution(const SampledDistribution& dist);
    SampledDistribution(QVector<qreal> samples,
                        QVector<qreal> weights = QVector<qreal>());
    virtual ~SampledDistribution();

    SampledDistribution* clone() const;

    qreal mean();
    qreal variance();
    qreal median();
    qreal percentile(qreal p);

    using IDistribution::pdf;
    using IDistribution::cdf;
    qreal pdf(qreal z);
    qreal cdf(qreal z);

    QVector<qreal> sampledVec() const;
    QVector<qreal> weights()    const;

    void addSamples(QVector<qreal> samples,
                    QVector<qreal> weights = QVector<qreal>());
    void addSample(qreal samp, qreal weight = 1.0);
    void removeSample(int idx);
    void clear();

    void generateZ();
    void generatePDF();
    void calculateCDF();

private:
    void update();

    QVector<qreal> m_samples;
    QVector<qreal> m_weights;
    qreal          m_mean;
    qreal          m_variance;
    bool           m_updated;
};

} // namespace Representation
} // namespace Tigon

#endif // SAMPLEDDISTRIBUTION_H
