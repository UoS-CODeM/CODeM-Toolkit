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

#include <vector>

namespace DTLZ {

std::vector<double > DTLZ1(const std::vector<double >& x, const int M);
std::vector<double > DTLZ2(const std::vector<double >& x, const int M);

}  // namespace DTLZ

#endif // DTLZPROBLEMS_H
