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
#include <core/distributions/MergedDistribution.h>
#include <core/utils/LinearInterpolator.h>

namespace Tigon {
namespace Representation {

MergedDistribution::MergedDistribution()
{
    m_type = Tigon::MergedDistType;
}

MergedDistribution::MergedDistribution(const MergedDistribution& dist)
    : IDistribution(dist)
{
    m_type = Tigon::MergedDistType;
    int n = dist.m_distributions.size();
    m_distributions.resize(n);
    m_ratios.resize(n);
    for(int i=0; i<m_distributions.size(); i++) {
        m_distributions[i].reset(dist.m_distributions[i]->clone());
        m_ratios[i] = dist.m_ratios[i];
    }

}

MergedDistribution::~MergedDistribution()
{

}

MergedDistribution* MergedDistribution::clone() const
{
    return (new MergedDistribution(*this));
}
void MergedDistribution::generateZ()
{
    int nDistributions = m_distributions.size();

    if(nDistributions == 0) {
        return;
    }

    for(int i=0; i<nDistributions; i++) {
        addZSamplesOfOneDistribution(m_distributions[i]);
    }
}

void MergedDistribution::addZSamplesOfOneDistribution(IDistributionSPtr d)
{
    QVector<qreal> newZ = d->zSamples();
    // merge newZ and m_z into augZ
    QVector<qreal> augZ;
    int iNew = 0; // iterator for newZ
    int iAug = 0; // iterator for m_z

    //TODO: avoid adding samples from lower frequency distributions
    //      to regions with higher freqwency samples
    while(iNew < newZ.size() && iAug < m_nSamples) {
        if(newZ[iNew] == m_z[iAug]) {
            augZ << m_z[iAug++];
            iNew++;
        }
        else if(newZ[iNew] < m_z[iAug]) {
            augZ << newZ[iNew++];
        }
        else {
            augZ << m_z[iAug++];
        }
    }

    while(iNew < newZ.size()) {
        augZ << newZ[iNew++];
    }

    while(iAug < m_nSamples) {
        augZ << m_z[iAug++];
    }

    m_z.swap(augZ);
    m_lb = m_z.first();
    m_ub = m_z.last();
    m_nSamples = m_z.size();
    augZ.clear();
}

void MergedDistribution::generatePDF()
{
    int nDistributions = m_distributions.size();

    if(nDistributions == 0) {
        return;
    }

    generateZ();

    m_pdf = QVector<qreal>(m_nSamples, 0.0);

    for(int i=0; i<nDistributions; i++) {
        addOnePDF(m_distributions[i], m_ratios[i]);
    }
    normalise();
}

void MergedDistribution::addOnePDF(IDistributionSPtr d, qreal ratio)
{
    // call this function only after the samples are integrated into m_z
    QVector<qreal> newPDF = d->pdf();
    QVector<qreal> newZ   = d->zSamples();

    // find the range of the new pdf
    int firstIdx = 0;
    int lastIdx  = m_nSamples - 1;

    while(newZ.first() > m_z[firstIdx]) {
        ++firstIdx;
    }
    while(newZ.last() < m_z[lastIdx]) {
        --lastIdx;
    }

    int nNewSamples = lastIdx - firstIdx + 1;

    LinearInterpolation* pdfInterpolator = new LinearInterpolation(newZ, newPDF);
    newPDF = pdfInterpolator->interpolateV(m_z.mid(firstIdx, nNewSamples));

    for(int i=0; i<nNewSamples; i++) {
        m_pdf[firstIdx+i] += (newPDF[i] * ratio);
    }
}

void MergedDistribution::appendDistribution(IDistributionSPtr d)
{
    appendDistribution(d, 1.0);
}

void MergedDistribution::appendDistribution(IDistributionSPtr d, qreal ratio)
{
    m_distributions.append(d);
    m_ratios.append(std::max(ratio, 0.0));
    if(!m_pdf.isEmpty()) {
        generatePDF();
    } else if(m_nSamples > 0) {
        addZSamplesOfOneDistribution(d);
    }
}

void MergedDistribution::removeDistribution(IDistributionSPtr d)
{
    int idx = m_distributions.indexOf(d);
    if(idx >= 0) {
        removeDistribution(idx);
    }
}

void MergedDistribution::removeDistribution(int idx)
{
    m_distributions.remove(idx);
    m_ratios.remove(idx);
    if(!m_z.isEmpty()) {
        generateZ();
    }
    if(!m_pdf.isEmpty()) {
        generatePDF();
    }
}

void MergedDistribution::changeRatio(IDistributionSPtr d, qreal newRatio)
{
    int idx = m_distributions.indexOf(d);
    if(idx >= 0) {
        changeRatio(idx, newRatio);
    }
}

void MergedDistribution::changeRatio(int idx, qreal newRatio)
{
    m_ratios[idx] = std::max(newRatio, 0.0);
}

} // namespace Representation
} // namespace Tigon
