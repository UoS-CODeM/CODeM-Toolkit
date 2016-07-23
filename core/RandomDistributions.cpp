/****************************************************************************
**
** The MIT License (MIT)
**
** Copyright (c) 2016 The University of Sheffield (www.sheffield.ac.uk)
**
** Permission is hereby granted, free of charge, to any person obtaining a copy
** of this software and associated documentation files (the "Software"), to deal
** in the Software without restriction, including without limitation the rights
** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
** copies of the Software, and to permit persons to whom the Software is
** furnished to do so, subject to the following conditions:
**
** The above copyright notice and this permission notice shall be included in all
** copies or substantial portions of the Software.
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
** OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
** SOFTWARE
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
    : m_lb(0.0),
      m_ub(1.0),
      m_dz(1.0),
      m_nSamples(0),
      m_quantileInterpolator(0),
      m_updated(false)
{

}

IDistribution::~IDistribution()
{
    if(m_quantileInterpolator != 0) {
        delete m_quantileInterpolator;
    }
}

double IDistribution::sample()
{
    if(!m_updated) {
        computeDistribution();
    }

    // A value between 0-1: 0==>lb , 1==>ub
    double sample = m_quantileInterpolator->interpolate(randUni());
    return sample;
}

vector<double> IDistribution::zSamples()
{
    if(!m_updated) {
        computeDistribution();
    }

    return m_z;
}

vector<double> IDistribution::pdf()
{
    if(!m_updated) {
        computeDistribution();
    }

    return m_pdf;
}

vector<double> IDistribution::cdf()
{
    if(!m_updated) {
        computeDistribution();
    }

    return m_cdf;
}

void IDistribution::defineResolution(double dz)
{
    if(dz > 0) {
        // make sure the range is a multiple of m_dz, and m_dz <= dz
        m_dz = (m_ub-m_lb) / ceil((m_ub-m_lb)/dz);
    }

    m_updated = false;
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

    m_updated = false;
}

double IDistribution::lowerBound() const
{
    return m_lb;
}

double IDistribution::upperBound() const
{
    return m_ub;
}

void IDistribution::computeDistribution()
{
    m_z.clear();
    m_pdf.clear();
    m_cdf.clear();

    generateZ();
    generatePDF();
    calculateCDF();

    if(m_quantileInterpolator == 0) {
        m_quantileInterpolator = new LinearInterpolator(m_cdf, m_z);
    } else {
        m_quantileInterpolator->defineXY(m_cdf, m_z);
    }

    m_updated = true;
}

void IDistribution::calculateCDF()
{
    m_cdf.assign(m_nSamples, 0.0);
    double cur  = 0.0;
    double next = 0.0;
    for(int i=0; i<m_nSamples-1; i++) {
        cur  = m_pdf[i];
        next = m_pdf[i+1];
        m_cdf[i+1] = m_cdf[i] + (cur+next)/2.0 * (m_z[i+1] - m_z[i]);
    }

    // normalise
    double factor = m_cdf.back();
    if(factor == 1.0) {
        return;
    } else if(factor == 0.0) {
        double probability = 1.0/(m_ub - m_lb);
        m_pdf.assign(m_nSamples, probability);
        calculateCDF();
    } else {
//        for(int i=0; i<m_cdf.size(); i++) {
//            m_pdf[i] /= factor;
//            m_cdf[i] /= factor;
//        }
        // Try this instead. TODO: remove the comment if it works
        transform(m_pdf.begin(), m_pdf.end(), m_pdf.begin(),
                  [factor](double p){return p/factor;});
        transform(m_cdf.begin(), m_cdf.end(), m_cdf.begin(),
                  [factor](double p){return p/factor;});
    }
}

void IDistribution::generateEquallySpacedZ()
{
    m_nSamples = (int)((m_ub - m_lb) / m_dz) + 1;
    m_z.resize(m_nSamples);
    double zz = m_lb;
    for(int i = 0; i < m_z.size() - 1; i++) {
        m_z[i] = zz;
        zz += m_dz;
    }
    m_z[m_z.size() - 1] = m_ub;
}

/// UNIFORM DISTRIBUTION
UniformDistribution::UniformDistribution()
{

}

UniformDistribution::UniformDistribution(double lb, double ub)
{
    defineBoundaries(lb, ub);
    defineResolution(m_ub - m_lb);
}

UniformDistribution::~UniformDistribution()
{

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
    double probability = 1.0/(m_ub - m_lb);
    m_pdf.assign(m_nSamples, probability);
}

/// LINEAR DISTRIBUTION
LinearDistribution::LinearDistribution()
    : m_ascend(true)
{
    defineResolution((m_ub-m_lb) / 2.0);
}

LinearDistribution::LinearDistribution(double lb, double ub, bool ascend)
    : m_ascend(ascend)
{
    defineBoundaries(lb, ub);
    defineResolution((m_ub-m_lb) / 2.0);
}

LinearDistribution::~LinearDistribution()
{

}

double LinearDistribution::sample()
{
    double r = randUni();
    double samp;
    if(m_ascend) {
        samp = m_lb + sqrt(r) * (m_ub - m_lb);
    } else {
        samp = m_ub - sqrt(1 - r) * (m_ub - m_lb);
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
    if(isAscend()) {
        for(int i=0; i<m_nSamples; i++) {
            m_pdf[i] = maxProbability * (m_z[i] - m_lb) / (m_ub - m_lb);
        }
    } else {
        for(int i=0; i<m_nSamples; i++) {
            m_pdf[i] = maxProbability * (m_ub - m_z[i]) / (m_ub - m_lb);
        }
    }
}

bool LinearDistribution::isAscend() const
{
    return m_ascend;
}

void LinearDistribution::defineAscend(bool a)
{
    if(m_ascend != a) {
        m_updated = false;
    }
    m_ascend = a;
}

/// PEAK DISTRIBUTION
PeakDistribution::PeakDistribution()
    : m_tendency(0.5),
      m_locality(1.0)
{
    defineResolution(1.0 / (DistNSamples - 1));
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
    m_pdf.resize(m_nSamples);

    double shift = PI * m_tendency;
    double N = DistPeakMinN + m_locality
            * (DistPeakMaxN - DistPeakMinN);
    vector<complex<double> > psiN(m_nSamples, complex<double>(0, 0));

    double nMax = max(3*N, DistPeakMinNBasisFunc);
    nMax = min(nMax, DistPeakMaxNBasisFunc);
    for(int n = 1; n <= nMax; n++) {
        double cNn = sqrt( (pow(N, n) * exp(-N) / factorial(n)));
        vector<double> psi = eigenFunction(n);
        for(int i=0; i<m_nSamples; i++) {
            complex<double> j(0.0, 1.0);
            psiN[i] += cNn * exp(-j * shift * (n + 0.5)) * psi[i];
        }
    }
    for(int i=0; i<m_nSamples; i++) {
        m_pdf[i] = real(conj(psiN[i]) * psiN[i]);
    }
}

vector<double> PeakDistribution::eigenFunction(int n)
{
    double Lz = m_ub - m_lb;
    double An = sqrt(2.0 / Lz);

    vector<double> psi(m_nSamples, 0.0);
    for(int i=1; i<m_nSamples-1; i++) {
        psi[i] = An * sin(PI * n * (m_z[i] - m_lb) / Lz);
    }
    return psi;
}


/// MERGED DISTRIBUTION
MergedDistribution::MergedDistribution()
{

}

MergedDistribution::~MergedDistribution()
{
    for(int i = 0; i < m_distributions.size(); i++) {
        delete m_distributions[i];
    }
}

void MergedDistribution::generateZ()
{
    int nDistributions = (int)m_distributions.size();

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
    if(newZ.empty())  {
        return;
    }
    vector<double>::iterator iNew = newZ.begin();
    vector<double>::iterator iExist = m_z.begin();
    // merge newZ and m_z into augZ
    vector<double> augZ;

//    int iNew = 0; // iterator for newZ
//    int iAug = 0; // iterator for m_z

    while((iNew != newZ.end()) && (iExist != m_z.end())) {
        if(*iNew == *iExist) {
            augZ.push_back(*iExist++);
            iNew++;
        }
        else if(*iNew < *iExist) {
            augZ.push_back(*iNew++);
        }
        else {
            augZ.push_back(*iExist++);
        }
    }

    while(iNew != newZ.end()) {
        augZ.push_back(*iNew++);
    }

    while(iExist != m_z.end()) {
        augZ.push_back(*iExist++);
    }

    m_z.swap(augZ);
    m_lb = m_z.front();
    m_ub = m_z.back();
    m_nSamples = (int)m_z.size();
    augZ.clear();
}

void MergedDistribution::generatePDF()
{
    size_t nDistributions = m_distributions.size();

    if(nDistributions == 0) {
        return;
    }

    m_pdf.resize(m_nSamples);

    for(int i=0; i<nDistributions; i++) {
        addOnePDF(m_distributions[i], m_ratios[i]);
    }
}

void MergedDistribution::addOnePDF(IDistribution* d, double ratio)
{
    // call this function only after the samples are integrated into m_z
    vector<double> newPDF = d->pdf();
    if(newPDF.empty()) {
        return;
    }
    vector<double> newZ   = d->zSamples();

    // find the range of the new pdf
    vector<double>::iterator first = m_z.begin();
    vector<double>::iterator last  = m_z.end();
    vector<double>::iterator pdfIter  = m_pdf.begin();

    while(newZ.front() > *first) {
        ++first;
        ++pdfIter;
    }
    while(newZ.back() < *(last - 1)) {
        --last;
    }

    //Interpolate the new pdf over the new samples in its range
    LinearInterpolator pdfInterpolator(newZ, newPDF);
    vector<double> tmp = pdfInterpolator.interpolateV(vector<double>(first, last));
    newPDF.swap(tmp);

    // add the new pdf times its weight to the existing pdf
    vector<double>::iterator newIter;
    for(newIter = newPDF.begin(); newIter == newPDF.end(); ++newIter, ++pdfIter) {
        *pdfIter += (*newIter * ratio);
    }
}

void MergedDistribution::appendDistribution(IDistribution* d)
{
    appendDistribution(d, 1.0);
}

void MergedDistribution::appendDistribution(IDistribution* d, double ratio)
{
    m_distributions.push_back(d);
    m_ratios.push_back(max(ratio, 0.0));
    computeDistribution();
}

} // namespace CODeM
