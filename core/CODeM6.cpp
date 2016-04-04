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
#include <tigon/Representation/Functions/CODeM/CODeM6.h>
#include <tigon/Representation/Elements/IElement.h>
#include <tigon/Utils/IElementUtils.h>
#include <tigon/Representation/Functions/CODeM/CODeMProblems.h>

#include <QStringList>
#include <QtMath>

namespace Tigon {
namespace Representation {

CODeM6::CODeM6()
{
    isParallelizable();
    TP_defineNInputs(3);
    TP_defineNOutputs(2);
    QString name("CODeM6");
    QString description("CODeM6 benchmark problem. A scalable uncertain problem.\n"
               "Composed of a deterministic multimodal problem with a planar "
               "Pareto front (DTLZ1), and uncertainty in both radial "
               "and perpendicular directions.\n"
               "Radial uncertainty increase as the deterministic objective "
               "vector is moving away from the true pareto front, and as its "
               "components are more similar.");
    createFunctionProperties(name, description);
}

CODeM6::CODeM6(const CODeM6& func)
{
    isParallelizable();
    TP_defineNInputs(func.TP_nInputs());
    TP_defineNOutputs(func.TP_nOutputs());
    createFunctionProperties(func.name(), func.description());
}

CODeM6::~CODeM6()
{

}

void CODeM6::evaluate(QVector<IElementSPtr> inputs,
                               QVector<IElementSPtr> outputs)
{
    if((inputs.size() == m_nInputs) && (outputs.size() == m_nOutputs) &&
            (m_nInputs > m_nOutputs)) {
        QVector<qreal> iReal = IElementVecToRealVec(inputs);
        QVector<qreal> oReal = CODeM::CODeM6(iReal, m_nOutputs);

        for(int i=0; i<outputs.size(); i++) {
            outputs[i]->defineValue(oReal[i]);
        }
    }
}

void CODeM6::defineInputPrpts()
{
    QStringList          varNames;
    QStringList          varDescriptions;
    QVector<ElementType> typeVec;
    QStringList          varUnits;
    QVector<IElement>    lowerBounds(TP_nInputs(), IElement(RealType, 0.0));
    QVector<IElement>    upperBounds(TP_nInputs(), IElement(RealType, 1.0));

    for(int i = 0; i < TP_nInputs(); i++) {
        varNames.append("Input_Var_" + QString::number(i));
        typeVec.append(RealType);
        varUnits.append("");
        if(i < TP_nOutputs()-1) {
            varDescriptions.append("Direction related variable");
        } else {
            varDescriptions.append("Distance related variable");
        }
    }

    createInputProperties(varNames, varDescriptions, typeVec, varUnits,
                          lowerBounds, upperBounds);
}

void CODeM6::defineOutputPrpts()
{
    QStringList          varNames;
    QStringList          varDescriptions;
    QVector<ElementType> typeVec;
    QStringList          varUnits;

    for(int i = 0; i < TP_nOutputs(); i++) {
        varNames.append("Output_Var_" + QString::number(i));
        varDescriptions.append("Output_VarDesc_" + QString::number(i));
        typeVec.append(RealType);
        varUnits.append("");
    }
    createOutputProperties(varNames, varDescriptions, typeVec, varUnits);
}

} // namespace Representation
} // namespace Tigon
