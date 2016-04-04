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
#ifndef DTLZPROBLEMS_H
#define DTLZPROBLEMS_H

//#include <tigon/tigon_global.h>
//#include <tigon/tigonconstants.h>
// Qt Includes
#include <QVector>

//namespace Tigon {
//namespace Operators {
namespace DTLZ {

//QVector<double > LIGER_TIGON_EXPORT DTLZ1(const QVector<double >& x, const int M);
QVector<double > DTLZ1(const QVector<double >& x, const int M);
QVector<double > DTLZ2(const QVector<double >& x, const int M);

// TODO: implement these
//QVector<double > DTLZ3(const QVector<double >& x, const int M);
//QVector<double > DTLZ4(const QVector<double >& x, const int M);
//QVector<double > DTLZ5(const QVector<double >& x, const int M);
//QVector<double > DTLZ6(const QVector<double >& x, const int M);
//QVector<double > DTLZ7(const QVector<double >& x, const int M);

}  // namespace DTLZ
//}  // namespace Tigon
//}  // namespace Operators

#endif // DTLZPROBLEMS_H
