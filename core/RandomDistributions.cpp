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
#include <core/distributions/RandomDistributions.h>
#include <random>
#include <core/utils/LinearInterpolator.h>
//#include <tigon/Utils/TigonUtils.h>

namespace CODeM {
IDistribution::IDistribution()
{
    m_type = Tigon::GenericDistType;
    m_nSamples = 0;
    m_lb = 0;
    m_ub = 1;
    m_dz = (m_ub-m_lb)/(Tigon::DistMinNSamples - 1);
    m_pdfInterpolator      = 0;
    m_cdfInterpolator      = 0;
    m_quantileInterpolator = 0;
}

IDistribution::IDistribution(const IDistribution& dist)
{
    m_type = dist.m_type;
    m_dz = dist.m_dz;
    m_lb = dist.m_lb;
    m_ub = dist.m_ub;
    m_z = dist.m_z;
    m_pdf = dist.m_pdf;
    m_cdf = dist.m_cdf;
    m_nSamples = dist.m_nSamples;
    m_pdfInterpolator      = 0;
    m_cdfInterpolator      = 0;
    m_quantileInterpolator = 0;
}

IDistribution::IDistribution(double value)
{
    m_type = Tigon::GenericDistType;
    m_nSamples = 0;
    m_lb = 0;
    m_ub = 1;
    m_dz = (m_ub-m_lb)/(Tigon::DistMinNSamples - 1);
    m_pdfInterpolator      = 0;
    m_cdfInterpolator      = 0;
    m_quantileInterpolator = 0;

    defineBoundaries(value,value);
}

IDistribution::~IDistribution()
{
    if(m_pdfInterpolator != 0) {
        delete m_pdfInterpolator;
    }
    if(m_cdfInterpolator != 0) {
        delete m_cdfInterpolator;
    }
    if(m_quantileInterpolator != 0) {
        delete m_quantileInterpolator;
    }
}

IDistribution* IDistribution::clone() const
{
    return (new IDistribution(*this));
}

Tigon::DistributionType IDistribution::type() const
{
    return m_type;
}

vector<double> IDistribution::parameters()
{
    return vector<double>();
}

double IDistribution::sample()
{
    if(m_quantileInterpolator == 0) {
        m_quantileInterpolator = new LinearInterpolator(cdf(), zSamples());
    } else {
        m_quantileInterpolator->defineXY(cdf(), zSamples());
    }
    double r = TRAND.randUni();
    // A value between 0-1: 0==>lb , 1==>ub
    double sample = m_quantileInterpolator->interpolate(r);
    return sample;
}

double IDistribution::mean()
{
    if(m_pdf.isEmpty()) {
        calculateCDF();
    }
    double sum  = 0.0;
    for(int i=0; i<m_nSamples-1; i++) {
        double cur  = m_pdf[i] * m_z[i];
        double next = m_pdf[i+1] * m_z[i+1];
        sum += (cur+next)/2 * (m_z[i+1] - m_z[i]);
    }
    return sum;
}

double IDistribution::median()
{
    if(m_quantileInterpolator == 0) {
        m_quantileInterpolator = new LinearInterpolator(cdf(),zSamples());
    } else {
        m_quantileInterpolator->defineXY(cdf(),zSamples());
    }
    return m_quantileInterpolator->interpolate(0.5);
}

double IDistribution::percentile(double p)
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
        m_quantileInterpolator = new LinearInterpolator(cdf(),zSamples());
    } else {
        m_quantileInterpolator->defineXY(cdf(),zSamples());
    }
    return m_quantileInterpolator->interpolate(p);
}

double IDistribution::variance()
{
    double m = mean();
    double sum  = 0.0;
    for(int i=0; i<m_nSamples-1; i++) {
        double cur  = m_pdf[i] * m_z[i]*m_z[i];
        double next = m_pdf[i+1] * m_z[i+1]*m_z[i+1];
        sum += (cur+next)/2 * (m_z[i+1] - m_z[i]);
    }
    return sum - m*m;
}

double IDistribution::std()
{
    return sqrt(variance());
}

vector<double> IDistribution::pdf()
{
    if(m_pdf.isEmpty() || m_pdf.size() != m_nSamples) {
        generatePDF();
    }
    return m_pdf;
}

vector<double> IDistribution::cdf()
{
    if(m_cdf.isEmpty() || m_cdf.size() != m_nSamples) {
        calculateCDF();
    }
    return m_cdf;
}

vector<double> IDistribution::pdf(const vector<double> zVec)
{
    vector<double> ret(zVec.size());
    for(int i=0; i< zVec.size(); i++) {
        ret[i] = pdf(zVec[i]);
    }
    return ret;
}

vector<double> IDistribution::cdf(const vector<double> zVec)
{
    vector<double> ret(zVec.size());
    for(int i=0; i< zVec.size(); i++) {
        ret[i] = cdf(zVec[i]);
    }
    return ret;
}

double IDistribution::pdf(double z)
{
    if(m_pdf.isEmpty() || m_pdf.size() != m_nSamples) {
        generatePDF();
    }

    if(z < m_lb || z > m_ub) {
        return 0.0;
    }

    if(m_pdfInterpolator == 0) {
        m_pdfInterpolator = new LinearInterpolator(zSamples(), pdf());
    } else {
        m_pdfInterpolator->defineXY(zSamples(), pdf());
    }


    return m_pdfInterpolator->interpolate(z);
}

double IDistribution::cdf(double z)
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
            m_cdfInterpolator = new LinearInterpolator(zSamples(), cdf());
        } else {
            m_cdfInterpolator->defineXY(zSamples(), cdf());
        }
        return m_cdfInterpolator->interpolate(z);
    }
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
            ub = Tigon::DistMinInterval;
        } else if(lb > 0) {
            ub = lb * (1 + Tigon::DistMinInterval);
        } else {
            lb = ub * (1 + Tigon::DistMinInterval);
        }
    }

    // TODO: test scaling

    double oldRange = m_ub - m_lb;
    double newRange = ub - lb;
    double ratio    = newRange / oldRange;

    m_lb = lb;
    m_ub = ub;

    defineResolution(m_dz * ratio);

    if(!m_z.isEmpty()) {
        for(int i=0; i<m_nSamples; i++) {
            m_z[i] = lb + ratio * (m_z[i] - m_lb);
        }
    }

    if(!m_pdf.isEmpty()) {
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

void IDistribution::defineZ(vector<double> z)
{
    std::sort(z.begin(),z.end());
    QMutableVectorIterator<double> i(z);
    if(i.hasNext()) {
        i.next();
        while(i.hasNext()) {
            if(i.peekPrevious() == i.next()) {
                i.remove();
            }
        }
    }

    if(m_pdf.size() != z.size()) {
        m_pdf.clear();
    }

    if(z.size() >= 2) {
        m_z = z;
        m_lb = m_z.first();
        m_ub = m_z.last();
        m_nSamples = m_z.size();
    } else {
        m_z.clear();
        m_nSamples = 0;
        defineBoundaries(z[0],z[0]);
    }
}

vector<double> IDistribution::zSamples()
{
    if(m_z.isEmpty()) {
        generateZ();
    }
    return m_z;
}

void IDistribution::generateZ()

{
    generateEquallySpacedZ();
}

void IDistribution::generatePDF()
{
    // uniform distribution
    double probability = 1.0/(m_ub - m_lb);
    m_pdf.fill(probability, m_nSamples);
}

void IDistribution::generateEquallySpacedZ()
{
    m_nSamples = (int)((m_ub-m_lb)/m_dz) + 1;
    m_z.resize(m_nSamples);
    double zz=m_lb;
    for(int i=0; i<m_z.size()-1; i++) {
        m_z[i] = zz;
        zz += m_dz;
    }
    m_z[m_z.size()-1] = m_ub;
}

void IDistribution::calculateCDF()
{
    if(m_pdf.isEmpty()) {
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

void IDistribution::normalise()
{
    if(m_pdf.isEmpty()) {
        return;
    }

    calculateCDF();
}

void IDistribution::negate()
{
    double ub = m_ub;
    m_ub = -m_lb;
    m_lb = -ub;

    if(m_z.isEmpty()) {
        return;
    }
    vector<double> newZ(m_nSamples);
    for(int i=0; i<m_nSamples; i++) {
        newZ[i] = -m_z[m_nSamples-1-i];
    }
    m_z.swap(newZ);

    if(m_pdf.isEmpty()) {
        return;
    }
    vector<double> newPdf(m_nSamples);
    for(int i=0; i<m_nSamples; i++) {
        newPdf[i] = m_pdf[m_nSamples-1-i];
    }
    m_pdf.swap(newPdf);
    calculateCDF();
}

void IDistribution::add(double num)
{
    m_ub += num;
    m_lb += num;

    if(m_z.isEmpty()) {
        return;
    }
    for(int i=0; i<m_nSamples; i++) {
        m_z[i] += num;
    }

    if(!m_pdf.isEmpty()) {
        calculateCDF();
    }
}

void IDistribution::add(const IDistribution* other)
{
    double lbO = other->lowerBound();
    double ubO = other->upperBound();
    double lb = m_lb + lbO;
    double ub = m_ub + ubO;
    double dz = qMax(m_ub-m_lb, ubO-lbO)/Tigon::DistConvNSamples;

    int nSampT = (int)((m_ub-m_lb)/dz) + 1;
    int nSampO = (int)((ubO-lbO)/dz) + 1;
    int nSamples = nSampT + nSampO - 1;

    // corrections to dz
    dz = (ub-lb)/(nSamples-1);
    double dzT = (m_ub-m_lb) / (nSampT-1);
    double dzO = (ubO - lbO) / (nSampO-1);

    // evenly distributed samples for convoluted distribution
    vector<double> z(nSamples);
    double zz = lb;
    for(int i=0; i<nSamples-1; i++) {
        z[i] = zz;
        zz += dz;
    }
    z[nSamples-1] = ub;

    // and for original distributions
    vector<double> pdfT(nSampT);
    zz = m_lb;
    for(int i=0; i<nSampT-1; i++) {
        pdfT[i] = pdf(zz);
        zz += dzT;
    }
    pdfT[nSampT-1] = pdf(m_ub);

    vector<double> pdfO(nSampO);
    zz = lbO;
    for(int i=0; i<nSampO-1; i++) {
        pdfO[i] = other->pdf(zz);
        zz += dzO;
    }
    pdfO[nSampO-1] = other->pdf(ubO);

    // convolute the two pdfs
    m_pdf = conv(pdfT,pdfO);

    // update the distribution
    m_z.swap(z);
    m_ub = ub;
    m_lb = lb;
    m_dz = dz;
    m_nSamples = nSamples;

    normalise();
}

void IDistribution::subtract(double num)
{
    add(-num);
}

void IDistribution::subtract(const IDistribution* other)
{
    IDistribution* minusOther(other->clone());
    minusOther->negate();
    add(minusOther);
}

void IDistribution::multiply(double num)
{
    if(num == 0)
    {
        m_lb = 0.0;
        m_ub = Tigon::DistMinInterval;
        m_dz = Tigon::DistMinInterval/(Tigon::DistMinNSamples-1);
        m_z.clear();
        m_z << m_lb << (m_lb+m_ub)/2 << m_ub;
        m_pdf.fill(1.0,2);
        normalise();
        return;
    }

    if(num < 0) {
        negate();
        num = -num;
    }

    m_ub *= num;
    m_lb *= num;

    if(m_z.isEmpty()) {
        return;
    }
    for(int i=0; i<m_nSamples; i++) {
        m_z[i] *= num;
    }

    if(!m_pdf.isEmpty()) {
        calculateCDF();
    }
}

void IDistribution::multiply(const IDistribution* other)
{
    double lbO = other->lowerBound();
    double ubO = other->upperBound();
    double lb  = qMin(m_lb*lbO,qMin(m_lb*ubO,qMin(m_ub*lbO,m_ub*ubO)));
    double ub  = qMax(m_lb*lbO,qMax(m_lb*ubO,qMax(m_ub*lbO,m_ub*ubO)));

    int nSamples = Tigon::DistMultNSamples;
    double dz  = (ub-lb) / (nSamples-1);
    double dzT = (m_ub-m_lb) / (nSamples-1);

    // evenly distributed samples for joint multiplied distribution
    vector<double> z(nSamples);
    double zz = lb;
    for(int i=0; i<nSamples-1; i++) {
        z[i] = zz;
        zz += dz;
    }
    z[nSamples-1] = ub;

    // and for original distributions
    vector<double> zT(nSamples);
    zz = m_lb;
    for(int i=0; i<nSamples-1; i++) {
        zT[i] = zz;
        zz += dzT;
    }
    zT[nSamples-1] = m_ub;

    // multiply the two distributions
    vector<double> newPDF(nSamples);
    for(int i=0; i<nSamples; i++) {
        for(int j=0; j<nSamples; j++) {
            double zt = zT[j];
            if(qAbs(zt) >= dzT/2) {
                double zo = z[i]/zt;
                if(zo >= lbO && zo <= ubO) {
                    newPDF[i] += pdf(zt) * other->pdf(zo) / qAbs(zt);
                }
            } else {
                zt = -dzT/2;
                double zo = z[i]/zt;
                if(zo >= lbO && zo <= ubO) {
                    newPDF[i] += pdf(zt) * other->pdf(zo) / qAbs(zt) / 2.0;
                }
                zt = dzT/2;
                zo = z[i]/zt;
                if(zo >= lbO && zo <= ubO) {
                    newPDF[i] += pdf(zt) * other->pdf(zo) / qAbs(zt) / 2.0;
                }
            }
        }
        newPDF[i] *= dzT;
    }

    // update the distribution
    m_pdf.swap(newPDF);
    m_z.swap(z);
    m_ub = ub;
    m_lb = lb;
    m_dz = dz;
    m_nSamples = nSamples;

    normalise();
}

void IDistribution::divide(double num)
{
    if(num == 0)
    {
        return;
    }

    multiply(1.0/num);
}

void IDistribution::divide(const IDistribution* other)
{
    double lbO = other->lowerBound();
    double ubO = other->upperBound();
    double lb,ub;

    // Division by zero
    if((sgn(lbO) != sgn(ubO)) || (lbO == 0.0) || (ubO == 0.0)) {
        if((sgn(m_lb) != sgn(m_ub)) || (sgn(lbO) != sgn(ubO))){
            lb = Tigon::Lowest;
            ub =  Tigon::Highest;
        } else if(lbO == 0.0) {
            if(m_lb >= 0.0) {
                lb = 0.0;
                ub =  Tigon::Highest;
            } else {
                lb = Tigon::Lowest;
                ub =  0.0;
            }
        } else {
            if(m_lb >= 0.0) {
                lb = Tigon::Lowest;
                ub =  0.0;
            } else {
                lb = 0.0;
                ub =  Tigon::Highest;
            }
        }

        IDistribution::defineBoundaries(lb, ub);
        generateEquallySpacedZ();
        IDistribution::generatePDF();
        return;
    }

    lb  = qMin(m_lb/lbO,qMin(m_lb/ubO,qMin(m_ub/lbO,m_ub/ubO)));
    ub  = qMax(m_lb/lbO,qMax(m_lb/ubO,qMax(m_ub/lbO,m_ub/ubO)));


    int nSamples = Tigon::DistMultNSamples;
    double dz  = (ub-lb) / (nSamples-1);
    double dzO = (ubO-lbO) / (nSamples-1);

    // evenly distributed samples for joint divided distribution
    vector<double> z(nSamples);
    double zz = lb;
    for(int i=0; i<nSamples-1; i++) {
        z[i] = zz;
        zz += dz;
    }
    z[nSamples-1] = ub;

    // and for original distribution
    vector<double> zO(nSamples);
    zz = lbO;
    for(int i=0; i<nSamples-1; i++) {
        zO[i] = zz;
        zz += dzO;
    }
    zO[nSamples-1] = ubO;

    // divide the two distributions
    vector<double> newPDF(nSamples);
    for(int i=0; i<nSamples; i++) {
        for(int j=0; j<nSamples; j++) {
            double zo = zO[j];
                double zt = z[i]*zo;
                if(zt >= m_lb && zt <= m_ub) {
                    newPDF[i] += other->pdf(zo) * pdf(zt) * qAbs(zo);
                }
        }
        newPDF[i] *= dzO;
    }

    // update the distribution
    m_pdf.swap(newPDF);
    m_z.swap(z);
    m_ub = ub;
    m_lb = lb;
    m_dz = dz;
    m_nSamples = nSamples;

    normalise();
}

void IDistribution::reciprocal()
{
    double lb,ub;

    // Division by zero
    if( (sgn(m_lb) != sgn(m_ub)) || (m_lb == 0.0) || (m_ub == 0.0) ) {
        if(sgn(m_lb) != sgn(m_ub)) {
            lb = Tigon::Lowest;
            ub =  Tigon::Highest;
        } else if(m_lb == 0.0) {
            lb = 0.0;
            ub =  Tigon::Highest;
        } else {
            lb = Tigon::Lowest;
            ub =  0.0;
        }

        IDistribution::defineBoundaries(lb, ub);
        generateEquallySpacedZ();
        IDistribution::generatePDF();
        return;
    }

    lb = 1.0 / m_ub;
    ub = 1.0 / m_lb;

    int nSamples = Tigon::DistMultNSamples;
    double dz  = (ub-lb) / (nSamples-1);

    // evenly distributed samples for reciprocal distribution
    vector<double> z(nSamples);
    double zz = lb;
    for(int i=0; i<nSamples-1; i++) {
        z[i] = zz;
        zz += dz;
    }
    z[nSamples-1] = ub;

//    // and for original distribution
//    vector<double> zT(nSamples);
//    zz = m_lb;
//    for(int i=0; i<nSamples-1; i++) {
//        zT[i] = zz;
//        zz += dzT;
//    }
//    zT[nSamples-1] = m_ub;

    // invert the distribution
    vector<double> newPDF(nSamples);
    for(int i=0; i<nSamples; i++) {
        double zt = 1.0 / z[i];
        newPDF[i] = pdf(zt) * qAbs(zt*zt);
    }

    // update the distribution
    m_pdf.swap(newPDF);
    m_z.swap(z);
    m_ub = ub;
    m_lb = lb;
    m_dz = dz;
    m_nSamples = nSamples;

    normalise();
}


UniformDistribution::UniformDistribution()
{
    m_uniDist = 0;
    m_type = Tigon::UniformDistType;
    defineBoundaries(0.0, 1.0);
    defineResolution(m_ub-m_lb);
}

UniformDistribution::UniformDistribution(const UniformDistribution& dist)
    : IDistribution(dist)
{
    m_type = Tigon::UniformDistType;
    m_uniDist = new boost::math::uniform_distribution<double>(m_lb, m_ub);
}

UniformDistribution::UniformDistribution(double lb, double ub)
{
    m_uniDist = 0;
    m_type = Tigon::UniformDistType;
    defineBoundaries(lb, ub);
    defineResolution(m_ub-m_lb);
}

UniformDistribution::UniformDistribution(vector<double> parameters)
{
    m_uniDist = 0;
    m_type = Tigon::UniformDistType;
    double lb = 0.0;
    double ub = 1.0;
    if(parameters.size() > 0) {
        lb = parameters[0];
        if((parameters.size() > 1) && (parameters[1] > lb)) {
            ub = parameters[1];
        } else {
            ub = lb + Tigon::DistMinInterval;
        }
    }
    defineBoundaries(lb, ub);
    defineResolution((m_ub-m_lb)/(Tigon::DistNSamples-1));
}

UniformDistribution::~UniformDistribution()
{
    delete m_uniDist;
    m_uniDist = 0;
}

UniformDistribution* UniformDistribution::clone() const
{
    return (new UniformDistribution(*this));
}

void UniformDistribution::defineBoundaries(double lb, double ub)
{
    if(lb >= ub) {
        if(lb == 0) {
            ub = Tigon::DistMinInterval;
        } else if(lb > 0) {
            ub = lb * (1 + Tigon::DistMinInterval);
        } else {
            lb = ub * (1 + Tigon::DistMinInterval);
        }
    }

    IDistribution::defineBoundaries(lb, ub);

    if(m_uniDist != 0) {
        delete m_uniDist;
        m_uniDist = 0;
    }

    m_uniDist = new boost::math::uniform_distribution<double>(m_lb, m_ub);
}


double UniformDistribution::sample()
{
    return TRAND.randUni(m_ub - m_lb, m_lb);
}

double UniformDistribution::mean()
{
    return (m_ub + m_lb) / 2.0;
}

double UniformDistribution::median()
{
    return boost::math::median(*m_uniDist);
}

double UniformDistribution::percentile(double p)
{
    if(p >= 1.0) {
        return m_ub;
    } else if(p <= 0.0) {
        return m_lb;
    }
    return boost::math::quantile(*m_uniDist, p);
}

double UniformDistribution::variance()
{
    return boost::math::variance(*m_uniDist);
}

double UniformDistribution::std()
{
    return boost::math::standard_deviation(*m_uniDist);
}

double UniformDistribution::pdf(double z)
{
    return boost::math::pdf(*m_uniDist, z);
}

double UniformDistribution::cdf(double z)
{
    return boost::math::cdf(*m_uniDist, z);
}

void UniformDistribution::generateZ()
{
    generateEquallySpacedZ();
}

void UniformDistribution::generatePDF()
{
    if(m_z.isEmpty()) {
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



LinearDistribution::LinearDistribution()
{
    m_type = Tigon::LinearDistType;
    defineResolution((m_ub-m_lb)/(Tigon::DistNSamples-1));
    m_ascend = true;
}

LinearDistribution::LinearDistribution(const LinearDistribution& dist)
    : IDistribution(dist)
{
    m_type = Tigon::LinearDistType;
    m_ascend = dist.m_ascend;
}

LinearDistribution::LinearDistribution(double lb, double ub)
{
    m_type = Tigon::LinearDistType;
    defineBoundaries(lb, ub);
    defineResolution((m_ub-m_lb)/(Tigon::DistNSamples-1));
    m_ascend = true;
}

LinearDistribution::LinearDistribution(vector<double> parameters)
{
    m_type = Tigon::LinearDistType;
    double lb = 0.0;
    double ub = 1.0;
    m_ascend = true;
    if(parameters.size() > 0) {
        lb = parameters[0];
        if((parameters.size() > 1) && (parameters[1] > lb)) {
            ub = parameters[1];
            if((parameters.size() < 3) && (parameters[2] <= 0.0)) {
                m_ascend =  false;
            }
        } else {
            ub = lb + Tigon::DistMinInterval;
        }
    }
    defineBoundaries(lb, ub);
}

LinearDistribution::~LinearDistribution()
{

}

LinearDistribution* LinearDistribution::clone() const
{
    return (new LinearDistribution(*this));
}

double LinearDistribution::sample()
{
    double r = TRAND.randUni();
    double samp;
    if(m_ascend) {
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
    if(m_z.isEmpty()) {
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
    bool oldDir = m_ascend;
    m_ascend = a;
    if(m_ascend != oldDir && !m_pdf.isEmpty()) {
        generatePDF();
    }
}

vector<double> LinearDistribution::parameters()
{
    vector<double> params;
    params << lowerBound() << upperBound() << (double)isAscend();
    return params;
}



PeakDistribution::PeakDistribution()
{
    m_type = Tigon::PeakDistType;
    defineTendencyAndLocality(0.5, 1.0);
}

PeakDistribution::PeakDistribution(const PeakDistribution& dist)
    : IDistribution(dist)
{
    m_type = Tigon::PeakDistType;
    m_tendency = dist.m_tendency;
    m_locality = dist.m_locality;

}

PeakDistribution::PeakDistribution(double tendency, double locality)
{
    m_type = Tigon::PeakDistType;
    defineTendencyAndLocality(tendency, locality);
}

PeakDistribution::PeakDistribution(vector<double> parameters)
{
    m_type = Tigon::PeakDistType;
    double tendency = 0.5;
    double locality = 1.0;
    if(parameters.size() > 0) {
        tendency = parameters[0];
    }
    if(parameters.size() > 1) {
        locality = parameters[1];
    }
    defineTendencyAndLocality(tendency, locality);
}

PeakDistribution::~PeakDistribution()
{

}

PeakDistribution* PeakDistribution::clone() const
{
    return (new PeakDistribution(*this));
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
    defineResolution(1.0/(m_locality+0.1)/(Tigon::DistNSamples-1));
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
    using namespace std;

    if(m_z.isEmpty()) {
        generateZ();
    }

    m_pdf = vector<double>(m_nSamples);

    double shift = boost::math::constants::pi<double>() * m_tendency;
    double N = Tigon::DistPeakMinN + m_locality
            * (Tigon::DistPeakMaxN - Tigon::DistPeakMinN);
    vector<qcomplex> psiN(m_nSamples,qcomplex(0,0));
    double nMax = qMax(3*N,Tigon::DistPeakMinNBasisFunc);
    nMax = qMin(nMax,Tigon::DistPeakMaxNBasisFunc);
    for(double n=1.0; n<=nMax; n++) {
        double cNn = sqrt( (qPow(N, n) * exp(-N) / boost::math::factorial<double>(n)));
        vector<double> psi = eigenFunction(n);
        for(int i=0; i<m_nSamples; i++) {
            qcomplex j(0, 1);
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
    double An = qPow(2.0/Lz, 0.5);

    vector<double> psi(m_nSamples, 0);
    for(int i=1; i<m_nSamples-1; i++) {
        psi[i] = An*sin(boost::math::constants::pi<double>()*n*(m_z[i]-m_lb)/Lz);
    }
    return psi;
}

} // namespace CODeM
