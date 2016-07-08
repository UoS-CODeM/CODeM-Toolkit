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
#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <core/CODeMGlobal.h>

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
    std::vector<double> zSamples();
    std::vector<double> pdf();
    std::vector<double> cdf();

    virtual void        defineBoundaries(double lb, double ub);
    void                defineResolution(double            dz);
    double              resolution()                     const;
    double              lowerBound()                     const;
    double              upperBound()                     const;

protected:
    void                computeDistribution();
    virtual void        generateZ() = 0;
    virtual void        generatePDF() = 0;
    void                generateEquallySpacedZ();

    double              m_lb;
    double              m_ub;
    std::vector<double> m_z;
    std::vector<double> m_pdf;
    int                 m_nSamples;
    bool                m_updated;

private:
    void                calculateCDF();

    double                      m_dz;
    Utils::LinearInterpolator*  m_quantileInterpolator;
    std::vector<double>         m_cdf;

};


class UniformDistribution : public IDistribution
{
public:
    UniformDistribution();
    UniformDistribution(double lb, double ub);
    virtual ~UniformDistribution();

    double sample();

    void  generateZ();
    void  generatePDF();
};


class LinearDistribution : public IDistribution
{
public:
    LinearDistribution();
    LinearDistribution(double lb, double ub, bool ascend = true);
    virtual ~LinearDistribution();

    double sample();
    void   generateZ();
    void   generatePDF();

    bool isAscend()     const;
    void defineAscend(bool a);

private:
    bool m_ascend;
};


class PeakDistribution : public IDistribution
{
public:
    PeakDistribution();
    PeakDistribution(double tendency, double locality);
    virtual ~PeakDistribution();

    void  defineTendencyAndLocality(double tendency, double locality);
    double tendency()  const;
    double locality()  const;

    void generateZ();
    void generatePDF();

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

private:
    std::vector<IDistribution*> m_distributions;
    std::vector<double>                m_ratios;

    void addZSamplesOfOneDistribution(IDistribution* d);
    void addOnePDF(IDistribution* d,       double ratio);
};

} // namespace CODeM
#endif // DISTRIBUTIONS_H
