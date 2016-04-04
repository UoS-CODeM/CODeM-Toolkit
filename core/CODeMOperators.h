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
#ifndef CODEMOPERATORS_H
#define CODEMOPERATORS_H
#include <tigon/tigon_global.h>
#include <tigon/tigonconstants.h>

// Qt Includes
#include <QVector>

namespace CODeM {

// Relations between UncertaintyKernel properties and uncertainty parameters
qreal LIGER_TIGON_EXPORT linearDecrease(qreal val);
qreal LIGER_TIGON_EXPORT skewedIncrease(qreal val, qreal alpha);
qreal LIGER_TIGON_EXPORT skewedDecrease(qreal val, qreal alpha);
qreal LIGER_TIGON_EXPORT lowOnValue(qreal val, qreal zeroVal, qreal width);
qreal LIGER_TIGON_EXPORT highOnValue(qreal val, qreal oneVal, qreal width);

QVector<qreal> LIGER_TIGON_EXPORT directionPerturbation(
        const QVector<qreal> oVec, qreal maxRadius, qreal pNorm=2);


} // namespace CODeM

#endif // CODEMOPERATORS_H
