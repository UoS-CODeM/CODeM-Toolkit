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

protected:
    virtual void generateZ() = 0;
    virtual void generatePDF() = 0;
    void computeDistribution();
    virtual void calculateCDF();
    void generateEquallySpacedZ();
    void normalise();

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
