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



// Qt Includes
#include <QVector>

namespace CODeM {

QVector<double>           CODeM1(QVector<double> iVec,
                                                   int k, int nObj);
QVector<QVector<double> > CODeM1(QVector<double> iVec,
                                                   int k, int nObj, int nSamp);
QVector<QVector<double> > CODeM1Perturb(QVector<double> oVec,
                                                          int nSamp = 1);

QVector<double>           CODeM2(QVector<double> iVec,
                                                   int k, int nObj);
QVector<QVector<double> > CODeM2(QVector<double> iVec,
                                                   int k, int nObj, int nSamp);
QVector<QVector<double> > CODeM2Perturb(QVector<double> oVec,
                                                          int nSamp = 1);

QVector<double>           CODeM3(QVector<double> iVec,
                                                   int k, int nObj);
QVector<QVector<double> > CODeM3(QVector<double> iVec,
                                                   int k, int nObj, int nSamp);
QVector<QVector<double> > CODeM3Perturb(QVector<double> oVec,
                                                          int nSamp = 1);

QVector<double>           CODeM4(QVector<double> iVec,
                                                   int k, int nObj);
QVector<QVector<double> > CODeM4(QVector<double> iVec,
                                                   int k, int nObj, int nSamp);
QVector<QVector<double> > CODeM4Perturb(QVector<double> oVec,
                                                          int nSamp = 1);

// CODeM5Perturb must have both decision and objective vectors defined
QVector<double>           CODeM5(QVector<double> iVec,
                                                   int k, int nObj);
QVector<QVector<double> > CODeM5(QVector<double> iVec,
                                                   int k, int nObj, int nSamp);
QVector<QVector<double> > CODeM5Perturb(QVector<double> iVec,
                                                          QVector<double> oVec,
                                                          int nSamp = 1);

QVector<double>           CODeM6(QVector<double> iVec,
                                                   int nObj);
QVector<QVector<double> > CODeM6(QVector<double> iVec,
                                                   int nObj, int nSamp);
QVector<QVector<double> > CODeM6Perturb(QVector<double> iVec,
                                                          QVector<double> oVec,
                                                          int nSamp = 1);

QVector<double> deterministicOVec(int prob,
                                                    QVector<double> iVec,
                                                    int nObj, int k=0);

BoxConstraintsData* createBoxConstraints(int prob, int nVar);
} // namespace CODeM

#endif // CODEMPROBLEMS_H
