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
#include <tigon/Representation/Distributions/SampledDistribution.h>
#include <tigon/Utils/ZeroOrderInterpolator.h>
#include <tigon/Utils/NormalisationUtils.h>
#include <tigon/Utils/TigonUtils.h>

#include <QtMath>

namespace Tigon {
namespace Representation {

SampledDistribution::SampledDistribution()
{
    m_type = SampledDistType;
    m_updated = false;
}

SampledDistribution::SampledDistribution(const SampledDistribution& dist)
    : IDistribution(dist)
{
    m_type = SampledDistType;
    m_samples = dist.m_samples;
    m_weights = dist.m_weights;
    m_mean = dist.m_mean;
    m_variance = dist.m_variance;
    m_updated = dist.m_updated;
}

SampledDistribution::SampledDistribution(QVector<qreal> samples,
                                         QVector<qreal> weights)
{
    m_type = SampledDistType;
    addSamples(samples, weights);
}

SampledDistribution::~SampledDistribution()
{

}

SampledDistribution* SampledDistribution::clone() const
{
    return (new SampledDistribution(*this));
}

qreal SampledDistribution::mean()
{
    if(!m_updated) {
        update();
    }
    return m_mean;
}

qreal SampledDistribution::variance()
{
    if(!m_updated) {
        update();
    }
    return m_variance;
}

QVector<qreal> SampledDistribution::sampledVec() const
{
    return m_samples;
}

qreal SampledDistribution::median()
{
    if(m_quantileInterpolator == 0) {
        m_quantileInterpolator = new ZeroOrderInterpolator(cdf(), zSamples());
    } else {
        m_quantileInterpolator->defineXY(cdf(), zSamples());
    }
    return m_quantileInterpolator->interpolate(0.5);
}

qreal SampledDistribution::percentile(qreal p)
{
    if(m_cdf.isEmpty() || m_cdf.size() != m_nSamples) {
        calculateCDF();
    }

    if(p>=1.0) {
        return m_ub;
    } else if(p<=0.0) {
        return m_lb;
    }

    if(m_quantileInterpolator == 0) {
        m_quantileInterpolator = new ZeroOrderInterpolator(cdf(), zSamples());
    } else {
        m_quantileInterpolator->defineXY(cdf(), zSamples());
    }
    return m_quantileInterpolator->interpolate(p);
}

qreal SampledDistribution::pdf(qreal z)
{
    int idx = zSamples().indexOf(z);
    if(idx >= 0) {
        return pdf().at(idx);
    } else {
        return 0.0;
    }
}

qreal SampledDistribution::cdf(qreal z)
{
    if(m_cdf.isEmpty() || m_cdf.size() != m_nSamples) {
        calculateCDF();
    }

    if(z <= m_lb) {
        return 0.0;
    } else if(z >= m_ub) {
        return 1.0;
    } else {
        if(m_cdfInterpolator == 0) {
            m_cdfInterpolator = new ZeroOrderInterpolator(zSamples(), cdf());
        } else {
            m_cdfInterpolator->defineXY(zSamples(), cdf());
        }
        return m_cdfInterpolator->interpolate(z);
    }
}


QVector<qreal> SampledDistribution::weights() const
{
    return m_weights;
}

void SampledDistribution::addSample(qreal samp, qreal weight)
{
    weight = qMax(weight, 0.0);
    int reSampIdx = m_samples.indexOf(samp);
    if(reSampIdx >= 0) {
        m_weights[reSampIdx] += weight;
    } else {
        m_samples.append(samp);
        m_weights.append(weight);
    }
    m_updated = false;
}

void SampledDistribution::addSamples(QVector<qreal> samples, QVector<qreal> weights)
{
    if(samples.size() == weights.size()) {
        for(int i=0; i<samples.size(); i++) {
            addSample(samples[i], weights[i]);
        }
    } else if(weights.isEmpty()) {
        for(int i=0; i<samples.size(); i++) {
            addSample(samples[i]);
        }
    } else {
        // ERROR -  different number of samples and weights
    }
}

void SampledDistribution::removeSample(int idx)
{
    if(m_samples.size() > idx && idx >= 0) {
        m_samples.remove(idx);
        m_weights.remove(idx);
    }
}

void SampledDistribution::clear()
{
    m_samples.clear();
    m_weights.clear();
    m_updated = false;
}

void SampledDistribution::generateZ()
{
    if(!(m_samples.isEmpty())) {
        m_nSamples = m_samples.size();

        // sort the samples and weights
        QVector<int> ind = indSort(m_samples);
        QVector<qreal> weights(m_nSamples);
        for(int i=0; i<m_nSamples; i++) {
            weights[i] = m_weights[ind[i]];
        }
        m_weights.swap(weights);

        m_z = m_samples;
        m_lb = m_z[0];
        m_ub = m_z[m_nSamples-1];
        defineResolution((m_ub-m_lb)/m_nSamples);
    } else {
        // ERROR - no samples
        return;
    }
}

void SampledDistribution::generatePDF()
{
    if(m_samples.isEmpty()) {
        // ERROR - no samples
        return;
    }
    if(m_z.isEmpty() || (m_nSamples != m_samples.size())) {
        generateZ();
    }

    m_pdf = m_weights;
    toUnitVec(m_pdf, 1.0);
}

void SampledDistribution::calculateCDF()
{
    if(m_samples.isEmpty()) {
        // ERROR - no samples
        return;
    }
    if(m_pdf.isEmpty() || (m_nSamples != m_samples.size())) {
        generatePDF();
    }

    m_cdf = m_pdf;
    for(int i=1; i<m_nSamples;i++) {
        m_cdf[i] += m_cdf[i-1];
    }
}

void SampledDistribution::update()
{
    if(m_samples.isEmpty()) {
        m_mean     = qQNaN();
        m_variance = qQNaN();
        return;
    }

    m_mean = 0.0;
    qreal totalWeights = 0.0;
    qreal tSquaredWeights = 0.0;
    // Calculate mean
    for(int i=0; i<m_samples.size(); i++) {
        m_mean += m_weights[i] * m_samples[i];
        totalWeights += m_weights[i];
        tSquaredWeights += m_weights[i] * m_weights[i];
    }
    m_mean /= totalWeights;

    // Calculate unbiased variance (from Wikipedia)
    m_variance  = 0.0;
    for(int i=0; i<m_samples.size(); i++) {
        m_variance += m_weights[i] * qPow(m_samples[i] - m_mean, 2.0);
    }
    m_variance /= (totalWeights - tSquaredWeights / totalWeights);

    m_updated = true;
}

} // namespace Representation
} // namespace Tigon
