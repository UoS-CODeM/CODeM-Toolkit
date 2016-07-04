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
#ifndef CODEMDISTRIBUTION_H
#define CODEMDISTRIBUTION_H

#include <codemglobal.h>
#include <vector>


namespace CODeM{
class IDistribution;
//namespace Utils {
//class LinearInterpolator;
//}

class CODeMDistribution
{
public:
    CODeMDistribution(IDistribution*             d,
                      const std::vector<double>& oVec,
                      double                     lowerBound,
                      double                     upperBound,
                      const std::vector<double>& ideal,
                      const std::vector<double>& antiIdeal,
                      double                     dirPertRad,
                      double                     dirPertNorm);
    ~CODeMDistribution();

    std::vector<double> sampleDistribution();

private:
    // 2-norm direction
    void defineDirection(const std::vector<double> &oVec);

    IDistribution*       m_distribution;
    double               m_directionPertRadius;
    std::vector<double>  m_direction;
    std::vector<double>  m_ideal;
    std::vector<double>  m_antiIdeal;
    double               m_lb;
    double               m_ub;
    double               m_pNorm;
};

} //namespace CODeM

#endif // CODEMDISTRIBUTION_H
