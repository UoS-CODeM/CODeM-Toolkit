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
#ifndef LINEARDISTRIBUTION_H
#define LINEARDISTRIBUTION_H
#include <tigon/tigon_global.h>
#include <tigon/tigonconstants.h>
#include <tigon/Representation/Distributions/IDistribution.h>

// Qt Includes
#include <QVector>
#include <QString>

namespace Tigon {
namespace Representation {

class LIGER_TIGON_EXPORT LinearDistribution : public IDistribution
{
public:
    LinearDistribution();
    LinearDistribution(const LinearDistribution& dist);
    LinearDistribution(qreal lb, qreal ub);
    LinearDistribution(QVector<qreal> parameters);
    virtual ~LinearDistribution();

    LinearDistribution *clone() const;

    qreal sample();
    void  generateZ();
    void  generatePDF();

    bool isAscend()     const;
    void defineAscend(bool a);

    QVector<qreal> parameters();

private:
    bool m_ascend;
};

} // namespace Representation
} // namespace Tigon

#endif // LINEARDISTRIBUTION_H
