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
#include <core/RandomDistributions.h>
#include <core/utils/LinearInterpolator.h>
#include <random>
#include <complex>
#include <math.h>
#include <algorithm>


using CODeM::Utils::LinearInterpolator;
using namespace std;

namespace CODeM {

/// IDISTRIBUTION
IDistribution::IDistribution()
{
    m_nSamples = 0;
    m_lb = 0;
    m_ub = 1;
    m_dz = (m_ub-m_lb)/(DistMinNSamples - 1);
    m_quantileInterpolator = 0;
}

IDistribution::~IDistribution()
{
    if(m_quantileInterpolator != 0) {
        delete m_quantileInterpolator;
    }
}

double IDistribution::sample()
{
    if(m_quantileInterpolator == 0) {
        m_quantileInterpolator = new LinearInterpolator(cdf(), zSamples());
    } else {
        m_quantileInterpolator->defineXY(cdf(), zSamples());
    }
    // A value between 0-1: 0==>lb , 1==>ub
    double sample = m_quantileInterpolator->interpolate(randUni());
    return sample;
}

vector<double> IDistribution::pdf()
{
    if(m_pdf.empty() || m_pdf.size() != m_nSamples) {
        generatePDF();
    }
    return m_pdf;
}

vector<double> IDistribution::cdf()
{
    if(m_cdf.empty() || m_cdf.size() != m_nSamples) {
        calculateCDF();
    }
    return m_cdf;
}

void IDistribution::defineResolution(double dz)
{
    if(dz > 0) {
        // make sure the range is a multiple of m_dz, and m_dz <= dz
        m_dz = (m_ub-m_lb) / ceil((m_ub-m_lb)/dz);
    }
}

double IDistribution::resolution() const
{
    return m_dz;
}

void IDistribution::defineBoundaries(double lb, double ub)
{
    if(lb >= ub) {
        if(lb == 0) {
            ub = DistMinInterval;
        } else if(lb > 0) {
            ub = lb * (1 + DistMinInterval);
        } else {
            lb = ub * (1 + DistMinInterval);
        }
    }

    double oldRange = m_ub - m_lb;
    double newRange = ub - lb;
    double ratio    = newRange / oldRange;

    m_lb = lb;
    m_ub = ub;

    defineResolution(m_dz * ratio);

    if(!m_z.empty()) {
        for(int i=0; i<m_nSamples; i++) {
            m_z[i] = lb + ratio * (m_z[i] - m_lb);
        }
    }

    if(!m_pdf.empty()) {
        normalise();
    }
}

double IDistribution::lowerBound() const
{
    return m_lb;
}

double IDistribution::upperBound() const
{
    return m_ub;
}

vector<double> IDistribution::zSamples()
{
    if(m_z.empty()) {
        generateZ();
    }
    return m_z;
}

void IDistribution::calculateCDF()
{
    if(m_pdf.empty()) {
        generatePDF();
    }
    m_cdf.fill(0.0, m_nSamples);
    double cur  = 0.0;
    double next = 0.0;
    for(int i=0; i<m_nSamples-1; i++) {
        cur  = m_pdf[i];
        next = m_pdf[i+1];
        m_cdf[i+1] = m_cdf[i] + (cur+next)/2 * (m_z[i+1] - m_z[i]);
    }

    // normalise
    double factor = m_cdf.last();
    if(factor == 1.0) {
        return;
    } else if(factor == 0.0) {
        double probability = 1.0/(m_ub - m_lb);
        m_pdf = vector<double>(m_nSamples, probability);
        calculateCDF();
    } else {
        for(int i=0; i<m_cdf.size(); i++) {
            m_pdf[i] /= factor;
            m_cdf[i] /= factor;
        }
    }
}

void IDistribution::generateEquallySpacedZ()
{
    m_nSamples = (int)((m_ub-m_lb)/m_dz) + 1;
    m_z.resize(m_nSamples);
    double zz = m_lb;
    for(int i=0; i<m_z.size()-1; i++) {
        m_z[i] = zz;
        zz += m_dz;
    }
    m_z[m_z.size()-1] = m_ub;
}

void IDistribution::normalise()
{
    if(m_pdf.empty()) {
        return;
    }

    calculateCDF();
}


/// UNIFORM DISTRIBUTION
UniformDistribution::UniformDistribution()
{
    m_uniDist = 0;
    defineBoundaries(0.0, 1.0);
    defineResolution(m_ub-m_lb);
}

UniformDistribution::UniformDistribution(double lb, double ub)
{
    m_uniDist = 0;
    defineBoundaries(lb, ub);
    defineResolution(m_ub-m_lb);
}

UniformDistribution::~UniformDistribution()
{
    delete m_uniDist;
    m_uniDist = 0;
}

void UniformDistribution::defineBoundaries(double lb, double ub)
{
    if(lb >= ub) {
        if(lb == 0) {
            ub = DistMinInterval;
        } else if(lb > 0) {
            ub = lb * (1 + DistMinInterval);
        } else {
            lb = ub * (1 + DistMinInterval);
        }
    }

    IDistribution::defineBoundaries(lb, ub);

    if(m_uniDist != 0) {
        delete m_uniDist;
        m_uniDist = 0;
    }
}

double UniformDistribution::sample()
{
    return m_lb + randUni() * (m_ub - m_lb);
}

void UniformDistribution::generateZ()
{
    generateEquallySpacedZ();
}

void UniformDistribution::generatePDF()
{
    if(m_z.empty()) {
        generateZ();
    }

    double probability = 1.0/(m_ub - m_lb);
    m_pdf = vector<double>(m_nSamples, probability);
}

vector<double> UniformDistribution::parameters()
{
    vector<double> params;
    params << lowerBound() << upperBound();
    return params;
}


/// LINEAR DISTRIBUTION
LinearDistribution::LinearDistribution()
{
    defineResolution((m_ub-m_lb)/2.0);
    m_increase = true;
}

LinearDistribution::LinearDistribution(const LinearDistribution& dist)
    : IDistribution(dist)
{
    m_increase = dist.m_increase;
}

LinearDistribution::LinearDistribution(double lb, double ub)
{
    defineBoundaries(lb, ub);
    defineResolution((m_ub-m_lb)/(DistNSamples-1));
    m_increase = true;
}

LinearDistribution::LinearDistribution(vector<double> parameters)
{
    double lb = 0.0;
    double ub = 1.0;
    m_increase = true;
    if(parameters.size() > 0) {
        lb = parameters[0];
        if((parameters.size() > 1) && (parameters[1] > lb)) {
            ub = parameters[1];
            if((parameters.size() < 3) && (parameters[2] <= 0.0)) {
                m_increase =  false;
            }
        } else {
            ub = lb + DistMinInterval;
        }
    }
    defineBoundaries(lb, ub);
}

LinearDistribution::~LinearDistribution()
{

}

double LinearDistribution::sample()
{
    double r = randUni();
    double samp;
    if(m_increase) {
        samp = m_lb + sqrt(r)*(m_ub-m_lb);
    } else {
        samp = m_ub - sqrt(1-r)*(m_ub-m_lb);
    }
    return samp;
}

void LinearDistribution::generateZ()
{
    generateEquallySpacedZ();
}

void LinearDistribution::generatePDF()
{
    if(m_z.empty()) {
        generateZ();
    }

    double maxProbability = 2/(m_ub - m_lb);
    m_pdf = vector<double>(m_nSamples);
    if(isIncreasing()) {
        for(int i=0; i<m_nSamples; i++) {
            m_pdf[i] = maxProbability * (m_z[i] - m_lb) / (m_ub - m_lb);
        }
    } else {
        for(int i=0; i<m_nSamples; i++) {
            m_pdf[i] = maxProbability * (m_ub - m_z[i]) / (m_ub - m_lb);
        }
    }
}

bool LinearDistribution::isIncreasing() const
{
    return m_increase;
}

void LinearDistribution::defineIncreasing(bool a)
{
    bool oldDir = m_increase;
    m_increase = a;
    if(m_increase != oldDir && !m_pdf.empty()) {
        generatePDF();
    }
}

vector<double> LinearDistribution::parameters()
{
    vector<double> params;
    params << lowerBound() << upperBound() << (double)isIncreasing();
    return params;
}


/// PEAK DISTRIBUTION
PeakDistribution::PeakDistribution()
{
    defineTendencyAndLocality(0.5, 1.0);
}

PeakDistribution::PeakDistribution(double tendency, double locality)
{
    defineTendencyAndLocality(tendency, locality);
}

PeakDistribution::~PeakDistribution()
{

}

void PeakDistribution::defineTendencyAndLocality(double tendency, double locality)
{
    if(tendency < 0.0) {
        tendency = 0.0;
    } else if(tendency > 1.0) {
        tendency = 1.0;
    }
    m_tendency = tendency;

    if(locality < 0.0) {
        locality = 0.0;
    } else if(locality > 1.0) {
        locality = 1.0;
    }

    m_locality = locality;

    //TODO: define high resolution at the peak and low at the rest
    defineResolution(1.0/(m_locality+0.1)/(DistNSamples-1));
}

double PeakDistribution::tendency() const
{
    return m_tendency;
}

double PeakDistribution::locality() const
{
    return m_locality;
}

void PeakDistribution::generateZ()
{
    generateEquallySpacedZ();
}

void PeakDistribution::generatePDF()
{
    if(m_z.empty()) {
        generateZ();
    }

    m_pdf = vector<double>(m_nSamples);

    double shift = PI * m_tendency;
    double N = DistPeakMinN + m_locality
            * (DistPeakMaxN - DistPeakMinN);
    vector<complex> psiN(m_nSamples, complex(0,0));
    double nMax = max(3*N, DistPeakMinNBasisFunc);
    nMax = min(nMax, DistPeakMaxNBasisFunc);
    for(int n=1; n<=nMax; n++) {
        double cNn = sqrt( (pow(N, n) * exp(-N) / tgamma(n+1)));
        vector<double> psi = eigenFunction(n);
        for(int i=0; i<m_nSamples; i++) {
            complex j(0, 1);
            psiN[i] += cNn * exp(-j*shift*(n+0.5)) * psi[i];
        }
    }
    for(int i=0; i<m_nSamples; i++) {
        m_pdf[i] = real(conj(psiN[i]) * psiN[i]);
    }
    normalise();
}

vector<double> PeakDistribution::parameters()
{
    vector<double> params;
    params << m_tendency << m_locality;
    return params;
}

vector<double> PeakDistribution::eigenFunction(int n)
{
    double Lz = m_ub - m_lb;
    double An = pow(2.0/Lz, 0.5);

    vector<double> psi(m_nSamples, 0);
    for(int i=1; i<m_nSamples-1; i++) {
        psi[i] = An*sin(boost::math::constants::pi<double>()*n*(m_z[i]-m_lb)/Lz);
    }
    return psi;
}


/// MERGED DISTRIBUTION
MergedDistribution::MergedDistribution()
{

}

MergedDistribution::~MergedDistribution()
{

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

void MergedDistribution::addZSamplesOfOneDistribution(IDistribution* d)
{
    vector<double> newZ = d->zSamples();
    // merge newZ and m_z into augZ
    vector<double> augZ;
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

    m_pdf = vector<double>(m_nSamples, 0.0);

    for(int i=0; i<nDistributions; i++) {
        addOnePDF(m_distributions[i], m_ratios[i]);
    }
    normalise();
}

void MergedDistribution::addOnePDF(IDistribution* d, double ratio)
{
    // call this function only after the samples are integrated into m_z
    vector<double> newPDF = d->pdf();
    vector<double> newZ   = d->zSamples();

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

    LinearInterpolator pdfInterpolator(newZ, newPDF);
    newPDF = pdfInterpolator.interpolateV(m_z.mid(firstIdx, nNewSamples));

    for(int i=0; i<nNewSamples; i++) {
        m_pdf[firstIdx+i] += (newPDF[i] * ratio);
    }
}

void MergedDistribution::appendDistribution(IDistribution* d)
{
    appendDistribution(d, 1.0);
}

void MergedDistribution::appendDistribution(IDistribution* d, double ratio)
{
    m_distributions.push_back(d);
    m_ratios.push_back(std::max(ratio, 0.0));
    if(!m_pdf.empty()) {
        generatePDF();
    } else if(m_nSamples > 0) {
        addZSamplesOfOneDistribution(d);
    }
    generatePDF();
}

void MergedDistribution::removeDistribution(IDistribution* d)
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
    if(!m_z.empty()) {
        generateZ();
    }
    if(!m_pdf.empty()) {
        generatePDF();
    }
}

void MergedDistribution::changeRatio(IDistribution* d, double newRatio)
{
    int idx = m_distributions.indexOf(d);
    if(idx >= 0) {
        changeRatio(idx, newRatio);
    }
}

void MergedDistribution::changeRatio(int idx, double newRatio)
{
    m_ratios[idx] = std::max(newRatio, 0.0);
}


} // namespace CODeM
