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
#ifndef MERGEDDISTRIBUTION_H
#define MERGEDDISTRIBUTION_H
#include <tigon/tigon_global.h>
#include <tigon/tigonconstants.h>
#include <tigon/Representation/Distributions/IDistribution.h>

// Qt Includes
#include <QVector>
#include <QString>

namespace Tigon {
namespace Representation {

class LIGER_TIGON_EXPORT MergedDistribution : public IDistribution
{
public:
    MergedDistribution();
    MergedDistribution(const MergedDistribution& dist);
    virtual ~MergedDistribution();

    MergedDistribution* clone() const;

    void generateZ();
    void generatePDF();

    void appendDistribution(IDistributionSPtr d);
    void appendDistribution(IDistributionSPtr d, qreal ratio);

    void removeDistribution(IDistributionSPtr d);
    void removeDistribution(int             idx);


    void changeRatio(IDistributionSPtr d, qreal newRatio);
    void changeRatio(int             idx, qreal newRatio);

private:
    QVector<IDistributionSPtr> m_distributions;
    QVector<qreal>             m_ratios;

    void addZSamplesOfOneDistribution(IDistributionSPtr d);
    void addOnePDF(IDistributionSPtr d,       qreal ratio);
};

} // namespace Representation
} // namespace Tigon

#endif // MERGEDDISTRIBUTION_H
