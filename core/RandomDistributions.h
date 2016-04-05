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

#include <vector>
using namespace std;

namespace CODeM {
class AbstractInterpolator;

class IDistribution
{
public:
    IDistribution();
    IDistribution(const IDistribution& dist);
    IDistribution(double sample);
    virtual ~IDistribution();

    virtual IDistribution* clone() const;

    Tigon::DistributionType type() const;

    virtual vector<double> parameters();

    virtual double sample();
    virtual double mean();
    virtual double variance();
    virtual double median();
    virtual double percentile(double p);
    virtual double std();

    vector<double>         pdf();
    vector<double>         cdf();
    virtual vector<double> pdf(const vector<double> zVec);
    virtual vector<double> cdf(const vector<double> zVec);
    virtual double          pdf(double                   z);
    virtual double          cdf(double                   z);

    virtual void   defineBoundaries(double lb, double ub);
    void           defineResolution(double           dz);
    double          resolution()                   const;
    double          lowerBound()                   const;
    double          upperBound()                   const;
    void           defineZ(vector<double>            z);
    vector<double> zSamples();

    virtual void generateZ();
    virtual void generatePDF();
    virtual void calculateCDF();
    void generateEquallySpacedZ();
    void normalise();

    // basic math operations
    virtual void negate();
    virtual void add(double num);
    virtual void add(const IDistribution* other);
    virtual void subtract(double num);
    virtual void subtract(const IDistribution* other);
    virtual void multiply(double num);
    virtual void multiply(const IDistribution* other);
    virtual void divide(double num);
    virtual void divide(const IDistribution* other);
    void reciprocal();


protected:
    Tigon::DistributionType  m_type;
    double                    m_dz;
    double                    m_lb;
    double                    m_ub;
    vector<double>           m_z;
    vector<double>           m_pdf;
    vector<double>           m_cdf;
    int                      m_nSamples;
    AbstractInterpolator*    m_quantileInterpolator;
    AbstractInterpolator*    m_pdfInterpolator;
    AbstractInterpolator*    m_cdfInterpolator;

};


class UniformDistribution : public IDistribution
{
public:
    UniformDistribution();
    UniformDistribution(const UniformDistribution& dist);
    UniformDistribution(double lb, double ub);
    UniformDistribution(vector<double> parameters);
    virtual ~UniformDistribution();

    UniformDistribution* clone() const;

    void defineBoundaries(double lb, double ub);

    double sample();
    double mean();
    double median();
    double percentile(double p);
    double variance();
    double std();

    //to make overrides visible to the compiler
    using IDistribution::pdf;
    using IDistribution::cdf;
    double pdf(double       z);
    double cdf(double       z);

    void  generateZ();
    void  generatePDF();

    vector<double> parameters();

private:
    boost::math::uniform_distribution<double>* m_uniDist;
};


class LinearDistribution : public IDistribution
{
public:
    LinearDistribution();
    LinearDistribution(const LinearDistribution& dist);
    LinearDistribution(double lb, double ub);
    LinearDistribution(vector<double> parameters);
    virtual ~LinearDistribution();

    LinearDistribution *clone() const;

    double sample();
    void  generateZ();
    void  generatePDF();

    bool isAscend()     const;
    void defineAscend(bool a);

    vector<double> parameters();

private:
    bool m_ascend;
};


class PeakDistribution : public IDistribution
{
public:
    PeakDistribution();
    PeakDistribution(const PeakDistribution& dist);
    PeakDistribution(double tendency, double locality);
    PeakDistribution(vector<double> parameters);
    virtual ~PeakDistribution();

    PeakDistribution* clone() const;

    void  defineTendencyAndLocality(double tendency, double locality);
    double tendency()  const;
    double locality()  const;

    void generateZ();
    void generatePDF();

    vector<double> parameters();

private:
    double m_tendency;
    double m_locality;

    vector<double> eigenFunction(int n);
};


class MergedDistribution : public IDistribution
{
public:
    MergedDistribution();
    MergedDistribution(const MergedDistribution& dist);
    virtual ~MergedDistribution();

    MergedDistribution* clone() const;

    void generateZ();
    void generatePDF();

    void appendDistribution(IDistribution* d);
    void appendDistribution(IDistribution* d, double ratio);

    void removeDistribution(IDistribution* d);
    void removeDistribution(int             idx);


    void changeRatio(IDistribution* d, double newRatio);
    void changeRatio(int             idx, double newRatio);

private:
    vector<IDistribution*> m_distributions;
    vector<double>             m_ratios;

    void addZSamplesOfOneDistribution(IDistribution* d);
    void addOnePDF(IDistribution* d,       double ratio);
};

} // namespace CODeM
#endif // DISTRIBUTIONS_H
