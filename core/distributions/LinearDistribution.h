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


#include <core/Distributions/IDistribution.h>

// Qt Includes
#include <QVector>
#include <QString>

class LinearDistribution : public IDistribution
{
public:
    LinearDistribution();
    LinearDistribution(const LinearDistribution& dist);
    LinearDistribution(double lb, double ub);
    LinearDistribution(QVector<double> parameters);
    virtual ~LinearDistribution();

    LinearDistribution *clone() const;

    double sample();
    void  generateZ();
    void  generatePDF();

    bool isAscend()     const;
    void defineAscend(bool a);

    QVector<double> parameters();

private:
    bool m_ascend;
};

#endif // LINEARDISTRIBUTION_H
