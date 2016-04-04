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



// Qt Includes
#include <QVector>
class BoxConstraintsData;

namespace CODeM {

class UncertaintyKernel
{
public:
    UncertaintyKernel(QVector<double> inputs,
                      QVector<double> outputs,
                      BoxConstraintsData* box);
    UncertaintyKernel(QVector<double> inputs,
                      QVector<double> outputs,
                      BoxConstraintsData* box,
                      double lb,
                      double ub);
    UncertaintyKernel(QVector<double> inputs,
                      QVector<double> outputs,
                      BoxConstraintsData* box,
                      QVector<double> ideal,
                      QVector<double> antiIdeal);
    UncertaintyKernel(QVector<double> inputs,
                      QVector<double> outputs,
                      BoxConstraintsData* box,
                      double lb,
                      double ub,
                      QVector<double> ideal,
                      QVector<double> antiIdeal);
    UncertaintyKernel(QVector<double> outputs,
                      double lb,
                      double ub,
                      QVector<double> ideal,
                      QVector<double> antiIdeal);
    ~UncertaintyKernel();

    double proximity();
    double symmetry();
    double oComponent(int idx) const;
    double dComponent(int idx) const;

    QVector<double> direction() const;

private:
    void defineIdealAndAntiIdeal(QVector<double> ideal, QVector<double> antiIdeal);
    // use normalised 2-norm values in objective space
    void defineDirectedObjectiveBoundaries(double lb, double ub);
    // set the lb to 0 and the ub to the directed boxed interval length
    void defineDirectedObjectiveBoundaries();

    // The design and objective values of the IMapping are not normalised
    QVector<double>       m_inputs;
    QVector<double>      m_outputs;
    BoxConstraintsData*  m_box;
    QVector<double>        m_ideal;
    QVector<double>    m_antiIdeal;
    QVector<double>    m_direction;
    double              m_distance;
    // normalised 2-norm values
    double                    m_ub;
    double                    m_lb;

    void calcDirectionAndDistance();
};

} // namespace CODeM

#endif // UNCERTAINTYKERNEL_H
