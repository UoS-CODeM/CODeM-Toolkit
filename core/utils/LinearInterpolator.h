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
#ifndef LINEARINTERPOLATOR_H
#define LINEARINTERPOLATOR_H
#include <tigon/tigon_global.h>
#include <tigon/tigonconstants.h>
#include <tigon/Utils/AbstractInterpolator.h>

// Qt Includes
#include <QtMath>
#include <QVector>

namespace Tigon {

class LIGER_TIGON_EXPORT LinearInterpolator : public AbstractInterpolator
{
public:
    LinearInterpolator(QVector<qreal> xv, QVector<qreal> yv);
    ~LinearInterpolator();

protected:
    qreal baseInterpolate(int j, qreal x);
};

} // namespace Tigon

#endif // LINEARINTERPOLATOR_H
