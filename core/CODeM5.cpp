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
#include <tigon/Representation/Functions/CODeM/CODeM5.h>
#include <tigon/Representation/Elements/IElement.h>
#include <tigon/Utils/IElementUtils.h>
#include <tigon/Representation/Functions/CODeM/CODeMProblems.h>

#include <QStringList>
#include <QtMath>

namespace Tigon {
namespace Representation {

CODeM5::CODeM5()
{
    TP_defineNInputs(3);
    TP_defineNOutputs(2);
    defineNumDirVars(0);
    QString name("CODeM5");
    QString description("CODeM5 benchmark problem. A scalable uncertain problem.\n"
               "Composed of a deterministic non-separable problem with a concave "
               "Pareto front (WFG8), and uncertainty in both radial and "
               "perpendicular directions.\n"
               "The uncertainty is proportional to the first decision variable.");

    createFunctionProperties(name, description);
}

CODeM5::CODeM5(const CODeM5& func)
{
    TP_defineNInputs(func.TP_nInputs());
    TP_defineNOutputs(func.TP_nOutputs());
    defineNumDirVars(func.numDirVars());
    createFunctionProperties(func.name(), func.description());
}

CODeM5::~CODeM5()
{

}

void CODeM5::defineNumDirVars(int k)
{
    if(k<=0) {
        // Default value
        k = qRound(m_nInputs / 1.2);
    }
    if(k>=m_nInputs) {
        // Largest value of k
        k = m_nInputs - 1;
    }

    // Make sure k % (m_nOutputs-1) == 0
    int dff = k % (m_nOutputs-1);
    m_nDirVar = k - dff;
}

int CODeM5::numDirVars() const
{
    return m_nDirVar;
}

void CODeM5::evaluate(QVector<IElementSPtr> inputs,
                               QVector<IElementSPtr> outputs)
{
    if((inputs.size() == m_nInputs) && (outputs.size() == m_nOutputs) &&
            (m_nInputs > m_nOutputs)) {
        QVector<qreal> iReal = IElementVecToRealVec(inputs);
        QVector<qreal> oReal = CODeM::CODeM5(iReal, m_nDirVar, m_nOutputs);

        for(int i=0; i<outputs.size(); i++) {
            outputs[i]->defineValue(oReal[i]);
        }
    }
}

void CODeM5::defineInputPrpts()
{
    QStringList          varNames;
    QStringList          varDescriptions;
    QVector<ElementType> typeVec;
    QStringList          varUnits;
    QVector<IElement>    lowerBounds(TP_nInputs(), IElement(RealType, 0.0));
    QVector<IElement>    upperBounds;

    for(int i = 0; i < TP_nInputs(); i++) {
        varNames.append("Input_Var_" + QString::number(i));
        typeVec.append(RealType);
        varUnits.append("");
        if(i < m_nDirVar) {
            varDescriptions.append("Direction related variable");
        } else {
            varDescriptions.append("Distance related variable");
        }
        upperBounds.append(IElement(RealType, 2.0*(i+1.0)));
    }

    createInputProperties(varNames, varDescriptions, typeVec, varUnits,
                          lowerBounds, upperBounds);
}

void CODeM5::defineOutputPrpts()
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
