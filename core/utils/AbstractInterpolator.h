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
#ifndef ABSTRACTINTERPOLATOR_H
#define ABSTRACTINTERPOLATOR_H
#include <tigon/tigon_global.h>
#include <tigon/tigonconstants.h>

// Qt Includes
#include <QtMath>
#include <QVector>

namespace Tigon {

class LIGER_TIGON_EXPORT AbstractInterpolator
{
public:
    explicit AbstractInterpolator(QVector<qreal> x, QVector<qreal> y, int m);
    virtual ~AbstractInterpolator();

    qreal interpolate(qreal xq);
    QVector<qreal> interpolateV(QVector<qreal> xq);
    virtual void defineXY(QVector<qreal> x, QVector<qreal> y);
    bool isConfigured();

protected:
    int locate(const qreal x);
    int hunt(const qreal x);
    virtual qreal baseInterpolate(int jlo, qreal x) = 0;
    virtual bool checkConfiguration();

    int n;
    int mm;
    int jsav;
    int cor;
    int dj;
    bool m_isConfigured;
    QVector<qreal> xx;
    QVector<qreal> yy;
};

} // namespace Tigon

#endif // ABSTRACTINTERPOLATOR_H
