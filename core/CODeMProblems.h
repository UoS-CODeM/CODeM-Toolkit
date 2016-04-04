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
#ifndef CODEMPROBLEMS_H
#define CODEMPROBLEMS_H
#include <tigon/tigon_global.h>
#include <tigon/tigonconstants.h>

// Qt Includes
#include <QVector>

using namespace Tigon::Representation;

namespace CODeM {

LIGER_TIGON_EXPORT QVector<qreal>           CODeM1(QVector<qreal> iVec,
                                                   int k, int nObj);
LIGER_TIGON_EXPORT QVector<QVector<qreal> > CODeM1(QVector<qreal> iVec,
                                                   int k, int nObj, int nSamp);
LIGER_TIGON_EXPORT QVector<QVector<qreal> > CODeM1Perturb(QVector<qreal> oVec,
                                                          int nSamp = 1);

LIGER_TIGON_EXPORT QVector<qreal>           CODeM2(QVector<qreal> iVec,
                                                   int k, int nObj);
LIGER_TIGON_EXPORT QVector<QVector<qreal> > CODeM2(QVector<qreal> iVec,
                                                   int k, int nObj, int nSamp);
LIGER_TIGON_EXPORT QVector<QVector<qreal> > CODeM2Perturb(QVector<qreal> oVec,
                                                          int nSamp = 1);

LIGER_TIGON_EXPORT QVector<qreal>           CODeM3(QVector<qreal> iVec,
                                                   int k, int nObj);
LIGER_TIGON_EXPORT QVector<QVector<qreal> > CODeM3(QVector<qreal> iVec,
                                                   int k, int nObj, int nSamp);
LIGER_TIGON_EXPORT QVector<QVector<qreal> > CODeM3Perturb(QVector<qreal> oVec,
                                                          int nSamp = 1);

LIGER_TIGON_EXPORT QVector<qreal>           CODeM4(QVector<qreal> iVec,
                                                   int k, int nObj);
LIGER_TIGON_EXPORT QVector<QVector<qreal> > CODeM4(QVector<qreal> iVec,
                                                   int k, int nObj, int nSamp);
LIGER_TIGON_EXPORT QVector<QVector<qreal> > CODeM4Perturb(QVector<qreal> oVec,
                                                          int nSamp = 1);

// CODeM5Perturb must have both decision and objective vectors defined
LIGER_TIGON_EXPORT QVector<qreal>           CODeM5(QVector<qreal> iVec,
                                                   int k, int nObj);
LIGER_TIGON_EXPORT QVector<QVector<qreal> > CODeM5(QVector<qreal> iVec,
                                                   int k, int nObj, int nSamp);
LIGER_TIGON_EXPORT QVector<QVector<qreal> > CODeM5Perturb(QVector<qreal> iVec,
                                                          QVector<qreal> oVec,
                                                          int nSamp = 1);

LIGER_TIGON_EXPORT QVector<qreal>           CODeM6(QVector<qreal> iVec,
                                                   int nObj);
LIGER_TIGON_EXPORT QVector<QVector<qreal> > CODeM6(QVector<qreal> iVec,
                                                   int nObj, int nSamp);
LIGER_TIGON_EXPORT QVector<QVector<qreal> > CODeM6Perturb(QVector<qreal> iVec,
                                                          QVector<qreal> oVec,
                                                          int nSamp = 1);

LIGER_TIGON_EXPORT QVector<qreal> deterministicOVec(int prob,
                                                    QVector<qreal> iVec,
                                                    int nObj, int k=0);

BoxConstraintsDataSPtr createBoxConstraints(int prob, int nVar);
} // namespace CODeM

#endif // CODEMPROBLEMS_H
