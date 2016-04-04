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
#include <core/Distributions/PeakDistribution.h>
#include <complex>
#include <algorithm>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/constants/constants.hpp>
#include <qmath.h>

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

PeakDistribution::PeakDistribution(QVector<double> parameters)
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

    m_pdf = QVector<double>(m_nSamples);

    double shift = boost::math::constants::pi<double>() * m_tendency;
    double N = Tigon::DistPeakMinN + m_locality
            * (Tigon::DistPeakMaxN - Tigon::DistPeakMinN);
    QVector<qcomplex> psiN(m_nSamples,qcomplex(0,0));
    double nMax = qMax(3*N,Tigon::DistPeakMinNBasisFunc);
    nMax = qMin(nMax,Tigon::DistPeakMaxNBasisFunc);
    for(double n=1.0; n<=nMax; n++) {
        double cNn = sqrt( (qPow(N, n) * exp(-N) / boost::math::factorial<double>(n)));
        QVector<double> psi = eigenFunction(n);
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

QVector<double> PeakDistribution::parameters()
{
    QVector<double> params;
    params << m_tendency << m_locality;
    return params;
}

QVector<double> PeakDistribution::eigenFunction(int n)
{
    double Lz = m_ub - m_lb;
    double An = qPow(2.0/Lz, 0.5);

    QVector<double> psi(m_nSamples, 0);
    for(int i=1; i<m_nSamples-1; i++) {
        psi[i] = An*sin(boost::math::constants::pi<double>()*n*(m_z[i]-m_lb)/Lz);
    }
    return psi;
}
