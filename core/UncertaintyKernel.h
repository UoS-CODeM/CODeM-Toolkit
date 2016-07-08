/****************************************************************************
**
** The MIT License (MIT)
**
** Copyright (c) 2016 The University of Sheffield (www.sheffield.ac.uk)
**
** Permission is hereby granted, free of charge, to any person obtaining a copy
** of this software and associated documentation files (the "Software"), to deal
** in the Software without restriction, including without limitation the rights
** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
** copies of the Software, and to permit persons to whom the Software is
** furnished to do so, subject to the following conditions:
**
** The above copyright notice and this permission notice shall be included in all
** copies or substantial portions of the Software.
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
** OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
** SOFTWARE
**
****************************************************************************/
#ifndef UNCERTAINTYKERNEL_H
#define UNCERTAINTYKERNEL_H

#include <core/CODeMGlobal.h>

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
    UncertaintyKernel(const vector<double> &outputs,
                      double lb,
                      double ub,
                      const vector<double> &ideal,
                      const vector<double> &antiIdeal);
    UncertaintyKernel(const vector<double> &inputs,
                      const vector<double> &outputs,
                      double lb,
                      double ub,
                      const vector<double> &inputsLowerBounds,
                      const vector<double> &inputsUpperBounds,
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
    void defineInputsBounds(const vector<double> &lowerBounds,
                            const vector<double> &upperBounds);
    // use normalised 2-norm values in objective space
    void defineDirectedObjectiveBoundaries(double lb, double ub);
    // set the lb to 0 and the ub to the directed boxed interval length
    void defineDirectedObjectiveBoundaries();

    // The design and objective values of the IMapping are not normalised
    vector<double>      m_inputs;
    vector<double>      m_outputs;
    vector<double>      m_ideal;
    vector<double>      m_antiIdeal;
    vector<double>      m_inLowerBounds;
    vector<double>      m_inUpperBounds;
    vector<double>      m_direction;
    double              m_distance;
    // normalised 2-norm values
    double              m_dirUpperBound;
    double              m_dirLowerBound;

    void calcDirectionAndDistance();
};

} // namespace CODeM

#endif // UNCERTAINTYKERNEL_H
