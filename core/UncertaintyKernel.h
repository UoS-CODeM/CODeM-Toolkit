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
#include <tigon/tigon_global.h>
#include <tigon/tigonconstants.h>

// Qt Includes
#include <QVector>
namespace Tigon {
namespace Representation {
class BoxConstraintsData;
}
}

using namespace Tigon::Representation;

namespace CODeM {

class LIGER_TIGON_EXPORT UncertaintyKernel
{
public:
    UncertaintyKernel(QVector<qreal> inputs,
                      QVector<qreal> outputs,
                      BoxConstraintsDataSPtr box);
    UncertaintyKernel(QVector<qreal> inputs,
                      QVector<qreal> outputs,
                      BoxConstraintsDataSPtr box,
                      qreal lb,
                      qreal ub);
    UncertaintyKernel(QVector<qreal> inputs,
                      QVector<qreal> outputs,
                      BoxConstraintsDataSPtr box,
                      QVector<qreal> ideal,
                      QVector<qreal> antiIdeal);
    UncertaintyKernel(QVector<qreal> inputs,
                      QVector<qreal> outputs,
                      BoxConstraintsDataSPtr box,
                      qreal lb,
                      qreal ub,
                      QVector<qreal> ideal,
                      QVector<qreal> antiIdeal);
    UncertaintyKernel(QVector<qreal> outputs,
                      qreal lb,
                      qreal ub,
                      QVector<qreal> ideal,
                      QVector<qreal> antiIdeal);
    ~UncertaintyKernel();

    qreal proximity();
    qreal symmetry();
    qreal oComponent(int idx) const;
    qreal dComponent(int idx) const;

    QVector<qreal> direction() const;

private:
    void defineIdealAndAntiIdeal(QVector<qreal> ideal, QVector<qreal> antiIdeal);
    // use normalised 2-norm values in objective space
    void defineDirectedObjectiveBoundaries(qreal lb, qreal ub);
    // set the lb to 0 and the ub to the directed boxed interval length
    void defineDirectedObjectiveBoundaries();

    // The design and objective values of the IMapping are not normalised
    QVector<qreal>       m_inputs;
    QVector<qreal>      m_outputs;
    BoxConstraintsDataSPtr  m_box;
    QVector<qreal>        m_ideal;
    QVector<qreal>    m_antiIdeal;
    QVector<qreal>    m_direction;
    qreal              m_distance;
    // normalised 2-norm values
    qreal                    m_ub;
    qreal                    m_lb;

    void calcDirectionAndDistance();
};

} // namespace CODeM

#endif // UNCERTAINTYKERNEL_H
