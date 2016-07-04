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
#ifndef UNCERTAINTYKERNEL_H
#define UNCERTAINTYKERNEL_H

#include <codemglobal.h>

#include <vector>
using std::vector;

namespace CODeM {

namespace Utils {
class BoxConstraintsData;
}
using namespace Utils;
class UncertaintyKernel
{
public:
    UncertaintyKernel(const vector<double> &inputs,
                      const vector<double> &outputs,
                      BoxConstraintsData* box);
    UncertaintyKernel(const vector<double> &inputs,
                      const vector<double> &outputs,
                      BoxConstraintsData* box,
                      double lb,
                      double ub);
    UncertaintyKernel(const vector<double> &inputs,
                      const vector<double> &outputs,
                      BoxConstraintsData* box,
                      const vector<double> &ideal,
                      const vector<double> &antiIdeal);
    UncertaintyKernel(const vector<double> &inputs,
                      const vector<double> &outputs,
                      BoxConstraintsData* box,
                      double lb,
                      double ub,
                      const vector<double> &ideal,
                      const vector<double> &antiIdeal);
    UncertaintyKernel(const vector<double> &outputs,
                      double lb,
                      double ub,
                      const vector<double> &ideal,
                      const vector<double> &antiIdeal);
    ~UncertaintyKernel();

    double proximity();
    double symmetry();
    double oComponent(int idx) const;
    double dComponent(int idx) const;

    vector<double> direction() const;

private:
    void defineIdealAndAntiIdeal(const vector<double> &ideal,
                                 const vector<double> &antiIdeal);
    // use normalised 2-norm values in objective space
    void defineDirectedObjectiveBoundaries(double lb, double ub);
    // set the lb to 0 and the ub to the directed boxed interval length
    void defineDirectedObjectiveBoundaries();

    // The design and objective values of the IMapping are not normalised
    vector<double>      m_inputs;
    vector<double>      m_outputs;
    BoxConstraintsData* m_box;
    vector<double>      m_ideal;
    vector<double>      m_antiIdeal;
    vector<double>      m_direction;
    double              m_distance;
    // normalised 2-norm values
    double              m_ub;
    double              m_lb;

    void calcDirectionAndDistance();
};

} // namespace CODeM

#endif // UNCERTAINTYKERNEL_H
