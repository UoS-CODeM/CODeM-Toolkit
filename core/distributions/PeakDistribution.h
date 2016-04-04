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
#include <tigon/tigon_global.h>
#include <tigon/tigonconstants.h>
#include <tigon/Representation/Distributions/IDistribution.h>

//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/lagged_fibonacci.hpp>
//#include <sstream>
//#include <string>


// Qt Includes
#include <QVector>
#include <QString>

namespace Tigon {
namespace Representation {

class LIGER_TIGON_EXPORT PeakDistribution : public IDistribution
{
public:
    PeakDistribution();
    PeakDistribution(const PeakDistribution& dist);
    PeakDistribution(qreal tendency, qreal locality);
    PeakDistribution(QVector<qreal> parameters);
    virtual ~PeakDistribution();

    PeakDistribution* clone() const;

    void  defineTendencyAndLocality(qreal tendency, qreal locality);
    qreal tendency()  const;
    qreal locality()  const;

    void generateZ();
    void generatePDF();

    QVector<qreal> parameters();

private:
    qreal m_tendency;
    qreal m_locality;

    QVector<qreal> eigenFunction(int n);
};

} // namespace Representation
} // namespace Tigon

#endif // PEAKDISTRIBUTION_H
