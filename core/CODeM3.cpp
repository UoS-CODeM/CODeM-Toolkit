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
#include <tigon/Representation/Functions/CODeM/CODeM3.h>
#include <tigon/Representation/Elements/IElement.h>
#include <tigon/Utils/IElementUtils.h>
#include <tigon/Representation/Functions/CODeM/CODeMProblems.h>

#include <QStringList>
#include <QtMath>

namespace Tigon {
namespace Representation {

CODeM3::CODeM3()
{
    TP_defineNInputs(3);
    TP_defineNOutputs(2);
    defineNumDirVars(0);
    QString name("CODeM3");
    QString description("CODeM3 benchmark problem. A scalable uncertain problem.\n"
               "Composed of a deterministic multimodal problem with a concave "
               "Pareto front (WFG4), and uncertainty in both radial and "
               "perpendicular directions.\n"
               "Radial uncertainty is lower close to the deterministic Pareto "
               "front and in regions with similar values of the objectives.\n"
               "Perependicular uncertainty is constant, exept for a region where the "
               "1st objective is close to 0.9");
    createFunctionProperties(name, description);
}

CODeM3::CODeM3(const CODeM3& func)
{
    TP_defineNInputs(func.TP_nInputs());
    TP_defineNOutputs(func.TP_nOutputs());
    defineNumDirVars(func.numDirVars());
    createFunctionProperties(func.name(), func.description());
}

CODeM3::~CODeM3()
{

}

void CODeM3::defineNumDirVars(int k)
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

int CODeM3::numDirVars() const
{
    return m_nDirVar;
}

void CODeM3::evaluate(QVector<IElementSPtr> inputs,
                               QVector<IElementSPtr> outputs)
{
    if((inputs.size() == m_nInputs) && (outputs.size() == m_nOutputs) &&
            (m_nInputs > m_nOutputs)) {
        QVector<qreal> iReal = IElementVecToRealVec(inputs);
        QVector<qreal> oReal = CODeM::CODeM3(iReal, m_nDirVar, m_nOutputs);

        for(int i=0; i<outputs.size(); i++) {
            outputs[i]->defineValue(oReal[i]);
        }
    }
}

void CODeM3::defineInputPrpts()
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

void CODeM3::defineOutputPrpts()
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
