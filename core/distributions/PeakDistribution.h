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
#ifndef PEAKDISTRIBUTION_H
#define PEAKDISTRIBUTION_H


#include <core/Distributions/IDistribution.h>

//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/lagged_fibonacci.hpp>
//#include <sstream>
//#include <string>


// Qt Includes
#include <QVector>
#include <QString>

class PeakDistribution : public IDistribution
{
public:
    PeakDistribution();
    PeakDistribution(const PeakDistribution& dist);
    PeakDistribution(double tendency, double locality);
    PeakDistribution(QVector<double> parameters);
    virtual ~PeakDistribution();

    PeakDistribution* clone() const;

    void  defineTendencyAndLocality(double tendency, double locality);
    double tendency()  const;
    double locality()  const;

    void generateZ();
    void generatePDF();

    QVector<double> parameters();

private:
    double m_tendency;
    double m_locality;

    QVector<double> eigenFunction(int n);
};

#endif // PEAKDISTRIBUTION_H
