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
#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <codemglobal.h>

#include <vector>

namespace CODeM {

namespace Utils {
class LinearInterpolator;
}

class IDistribution
{
public:
    IDistribution();
    virtual ~IDistribution();

    virtual double      sample();
    std::vector<double> pdf();
    std::vector<double> cdf();

    virtual void        defineBoundaries(double lb, double ub);
    void                defineResolution(double            dz);
    double              resolution()                     const;
    double              lowerBound()                     const;
    double              upperBound()                     const;
    std::vector<double> zSamples();

    virtual void generateZ() = 0;
    virtual void generatePDF() = 0;
    virtual void calculateCDF();
    void generateEquallySpacedZ();
    void normalise();

protected:
    double                      m_dz;
    double                      m_lb;
    double                      m_ub;
    std::vector<double>         m_z;
    std::vector<double>         m_pdf;
    std::vector<double>         m_cdf;
    int                         m_nSamples;
    Utils::LinearInterpolator*  m_quantileInterpolator;

};


class UniformDistribution : public IDistribution
{
public:
    UniformDistribution();
    UniformDistribution(double lb, double ub);
    virtual ~UniformDistribution();

    void defineBoundaries(double lb, double ub);

    double sample();

    void  generateZ();
    void  generatePDF();
};


class LinearDistribution : public IDistribution
{
public:
    LinearDistribution();
    LinearDistribution(double lb, double ub, bool increase);
    LinearDistribution(std::vector<double> parameters);
    virtual ~LinearDistribution();

    double sample();
    void   generateZ();
    void   generatePDF();

    bool isIncreasing()     const;
    void defineIncreasing(bool a);

    std::vector<double> parameters();

private:
    bool m_increase;
};


class PeakDistribution : public IDistribution
{
public:
    PeakDistribution();
    PeakDistribution(const PeakDistribution& dist);
    PeakDistribution(double tendency, double locality);
    PeakDistribution(vector<double> parameters);
    virtual ~PeakDistribution();

    void  defineTendencyAndLocality(double tendency, double locality);
    double tendency()  const;
    double locality()  const;

    void generateZ();
    void generatePDF();

    std::vector<double> parameters();

private:
    double m_tendency;
    double m_locality;

    std::vector<double> eigenFunction(int n);
};


class MergedDistribution : public IDistribution
{
public:
    MergedDistribution();
    MergedDistribution(const MergedDistribution& dist);
    virtual ~MergedDistribution();

    void generateZ();
    void generatePDF();

    void appendDistribution(IDistribution* d);
    void appendDistribution(IDistribution* d, double ratio);

    void removeDistribution(IDistribution* d);
    void removeDistribution(int          idx);


    void changeRatio(IDistribution* d, double newRatio);
    void changeRatio(int          idx, double newRatio);

private:
    std::vector<IDistribution*> m_distributions;
    std::vector<double>                m_ratios;

    void addZSamplesOfOneDistribution(IDistribution* d);
    void addOnePDF(IDistribution* d,       double ratio);
};

} // namespace CODeM
#endif // DISTRIBUTIONS_H
