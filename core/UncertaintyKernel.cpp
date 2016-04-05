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
#include <core/UncertaintyKernel.h>
#include <tigon/Representation/Mappings/IMapping.h>
#include <tigon/Representation/Elements/IElement.h>
#include <tigon/Representation/Constraints/BoxConstraintsData.h>
#include <tigon/Utils/NormalisationUtils.h>
#include <tigon/Utils/TigonUtils.h>
#include <math.h>

using std::vector;

namespace CODeM {

UncertaintyKernel::UncertaintyKernel(vector<double> inputs,
                                     vector<double> outputs,
                                     BoxConstraintsData* box)
{
    m_inputs  = inputs;
    m_outputs = outputs;
    m_box     = box->clone();
    defineIdealAndAntiIdeal(vector<double>(outputs.size(), 0.0),
                            vector<double>(outputs.size(), 1.0));
    calcDirectionAndDistance();
    defineDirectedObjectiveBoundaries();
}

UncertaintyKernel::UncertaintyKernel(vector<double> inputs,
                                     vector<double> outputs,
                                     BoxConstraintsData* box,
                                     double lb,
                                     double ub)
{
    m_inputs  = inputs;
    m_outputs = outputs;
    m_box     = box->clone();
    defineIdealAndAntiIdeal(vector<double>(outputs.size(), 0.0),
                            vector<double>(outputs.size(), 1.0));
    calcDirectionAndDistance();
    defineDirectedObjectiveBoundaries(lb, ub);
}

UncertaintyKernel::UncertaintyKernel(vector<double> inputs,
                                     vector<double> outputs,
                                     BoxConstraintsData* box,
                                     vector<double> ideal,
                                     vector<double> antiIdeal)
{
    m_inputs  = inputs;
    m_outputs = outputs;
    m_box     = box->clone();
    // default values in case of incorrect ideal and antiIdeal
    defineIdealAndAntiIdeal(vector<double>(outputs.size(), 0.0),
                            vector<double>(outputs.size(), 1.0));
    defineIdealAndAntiIdeal(ideal, antiIdeal);
    calcDirectionAndDistance();
    defineDirectedObjectiveBoundaries();
}

UncertaintyKernel::UncertaintyKernel(vector<double> inputs,
                                     vector<double> outputs,
                                     BoxConstraintsData* box,
                                     double lb,
                                     double ub,
                                     vector<double> ideal,
                                     vector<double> antiIdeal)
{
    m_inputs  = inputs;
    m_outputs = outputs;
    m_box     = box->clone();
    defineIdealAndAntiIdeal(vector<double>(outputs.size(), 0.0),
                            vector<double>(outputs.size(), 1.0));
    defineIdealAndAntiIdeal(ideal, antiIdeal);
    calcDirectionAndDistance();
    defineDirectedObjectiveBoundaries(lb, ub);
}

UncertaintyKernel::UncertaintyKernel(vector<double> outputs,
                                     double lb,
                                     double ub,
                                     vector<double> ideal,
                                     vector<double> antiIdeal)
{
    m_outputs = outputs;
    defineIdealAndAntiIdeal(vector<double>(outputs.size(), 0.0),
                            vector<double>(outputs.size(), 1.0));
    defineIdealAndAntiIdeal(ideal, antiIdeal);
    calcDirectionAndDistance();
    defineDirectedObjectiveBoundaries(lb, ub);
}


UncertaintyKernel::~UncertaintyKernel()
{

}

void UncertaintyKernel::calcDirectionAndDistance()
{
    if(m_outputs.size() < 2) {
        return;
    }

    m_direction = m_outputs;
    normaliseToUnitBox(m_direction, m_ideal, m_antiIdeal);

    // m_distance in 2-norm
    m_distance = magnitudeAndDirectionP(m_direction, 2.0);

    // normalise the direction vector to the k-1 simplex
    toUnitVec(m_direction, 1.0);
}

double UncertaintyKernel::proximity()
{
    if(m_distance <= m_lb) {
        return 0.0;
    } else if(m_distance >= m_ub) {
        return 1.0;
    } else {
        return (m_distance - m_lb) / (m_ub - m_lb);
    }
}

double UncertaintyKernel::symmetry()
{
    double euclideanDist = 0.0;
    for(int i=0; i < m_direction.size(); i++) {
        euclideanDist += m_direction[i] * m_direction[i];
    }
    euclideanDist = sqrt(euclideanDist);

    double symmetryVal = (1.0 - euclideanDist) /
            (1.0 - 1.0/sqrt(m_direction.size()));

    return pow(symmetryVal, 2.0);
}

double UncertaintyKernel::oComponent(int idx) const
{
    if(!isInRange(idx, m_outputs.size())) {
        return -1.0;
    }
    double oc = (m_outputs[idx]-m_ideal[idx]) / (m_antiIdeal[idx]-m_ideal[idx]);
    return oc;
}

double UncertaintyKernel::dComponent(int idx) const
{
    if(!isInRange(idx, m_inputs.size())) {
        return -1.0;
    }
    double lb = m_box->lowerBounds().at(idx).value<double>();
    double ub = m_box->upperBounds().at(idx).value<double>();
    double d  = m_inputs.at(idx);

    double dRatio = (d-lb)/(ub-lb);
    return dRatio;
}

vector<double> UncertaintyKernel::direction() const
{
    return m_direction;
}

void UncertaintyKernel::defineDirectedObjectiveBoundaries(double lb, double ub)
{
    if(ub > lb) {
        m_lb = lb;
        m_ub = ub;
    }
}

void UncertaintyKernel::defineDirectedObjectiveBoundaries()
{
    m_lb = 0.0;
    m_ub = directedBoxedIntervalLength(m_direction);
}

void UncertaintyKernel::defineIdealAndAntiIdeal(vector<double> ideal,
                                           vector<double> antiIdeal)
{
    if(ideal.size() == antiIdeal.size()) {
        for(int i=0; i<ideal.size(); i++) {
            if(antiIdeal[i] <= ideal[i]) {
                return;
            }
        }
        m_ideal     = ideal;
        m_antiIdeal = antiIdeal;
    }
}

} // namespace CODeM
