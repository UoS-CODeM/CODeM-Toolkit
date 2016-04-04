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
#include <core/CODeMProblems.h>
#include <core/UncertaintyKernel.h>
#include <core/CODeMOperators.h>
#include <core/CODeMDistribution.h>
#include <core/Distributions/MergedDistribution.h>
#include <core/Distributions/UniformDistribution.h>
#include <core/Distributions/PeakDistribution.h>
#include <tigon/Representation/Mappings/IMapping.h>
#include <tigon/Representation/Elements/IElement.h>
#include <tigon/Representation/Constraints/BoxConstraintsData.h>
#include <tigon/Utils/NormalisationUtils.h>
#include <libs/WFG/ExampleProblems.h>
#include <libs/DTLZ/DTLZProblems.h>

using namespace WFGT::Toolkit::Examples::Problems;

namespace CODeM {

QVector<double> CODeM1(QVector<double> iVec, int k, int nObj)

{
    // Evaluate the decision vector
    QVector<double> oVec =
            QVector<double>::fromStdVector(WFG4(iVec.toStdVector(), k, nObj));

    return CODeM1Perturb(oVec)[0];
}

QVector<QVector<double> > CODeM1(QVector<double> iVec,
                                int k, int nObj, int nSamp)
{
    // Evaluate the decision vector
    QVector<double> oVec =
            QVector<double>::fromStdVector(WFG4(iVec.toStdVector(), k, nObj));

    return CODeM1Perturb(oVec, nSamp);
}

QVector<QVector<double> > CODeM1Perturb(QVector<double> oVec, int nSamp)
{
    // Set the uncertainty kernel
    QVector<double> ideal(oVec.size(), 0);
    QVector<double> antiIdeal(oVec.size());
    for(int i=0; i<oVec.size(); i++) {
        antiIdeal[i] = 3*(i+1);
    }
    double lb = 2.0/3.0;
    double ub = 1.0;

    UncertaintyKernel uk(oVec, lb, ub, ideal, antiIdeal);

    // Evaluate the uncertainty parameters
    double peakTend, peakLoc, dirPertRad;

    peakTend = uk.proximity();
    peakLoc  = lowOnValue(uk.proximity(), 0.0, 0.05);

    double distanceNorm = 2.0;

    dirPertRad = 0.0;

    // Create the CODeM distribution
    PeakDistribution* d = PeakDistribution*(
                new PeakDistribution(peakTend, peakLoc));

    CODeMDistribution cd(d, oVec, lb, ub, ideal, antiIdeal, dirPertRad, distanceNorm);

    // Sample the distribution
    QVector<QVector<double> > samples;
    for(int i=0; i<nSamp; i++) {
        samples.append(cd.sampleDistribution());
    }
    return samples;
}

QVector<double> CODeM2(QVector<double> iVec, int k, int nObj)
{
    // Evaluate the decision vector
    QVector<double> oVec =
            QVector<double>::fromStdVector(WFG4(iVec.toStdVector(), k, nObj));

    return CODeM2Perturb(oVec)[0];
}

QVector<QVector<double> > CODeM2(QVector<double> iVec,
                                int k, int nObj, int nSamp)
{
    // Evaluate the decision vector
    QVector<double> oVec =
            QVector<double>::fromStdVector(WFG4(iVec.toStdVector(), k, nObj));

    return CODeM2Perturb(oVec, nSamp);
}

QVector<QVector<double> > CODeM2Perturb(QVector<double> oVec, int nSamp)
{
    // Set the uncertainty kernel
    QVector<double> ideal(oVec.size(), 0.0);
    QVector<double> antiIdeal(oVec.size());
    for(int i=0; i<oVec.size(); i++) {
        antiIdeal[i] = 3.0*(i+1);
    }
    double lb = 2.0/3.0;
    double ub = 1.0;

    UncertaintyKernel uk(oVec, lb, ub, ideal, antiIdeal);

    // Evaluate the uncertainty parameters
    double uniLB, uniUB, dirPertRad;

    uniLB  = uk.proximity();
    uniUB  = uniLB;

    double distanceNorm = 2.0;

    dirPertRad = 0.1 * uk.symmetry();

    // Create the CODeM distribution
    UniformDistribution* d = UniformDistribution*(
                new UniformDistribution(uniLB, uniUB));

    CODeMDistribution cd(d, oVec, lb, ub, ideal, antiIdeal, dirPertRad, distanceNorm);

    // Sample the distribution
    QVector<QVector<double> > samples;
    for(int i=0; i<nSamp; i++) {
        samples.append(cd.sampleDistribution());
    }
    return samples;
}

QVector<double> CODeM3(QVector<double> iVec, int k, int nObj)
{
    // Evaluate the decision vector
    QVector<double> oVec =
            QVector<double>::fromStdVector(WFG4(iVec.toStdVector(), k, nObj));

    return CODeM3Perturb(oVec)[0];
}

QVector<QVector<double> > CODeM3(QVector<double> iVec,
                                int k, int nObj, int nSamp)
{
    // Evaluate the decision vector
    QVector<double> oVec =
            QVector<double>::fromStdVector(WFG4(iVec.toStdVector(), k, nObj));

    return CODeM3Perturb(oVec, nSamp);
}

QVector<QVector<double> > CODeM3Perturb(QVector<double> oVec, int nSamp)
{
    // Set the uncertainty kernel
    QVector<double> ideal(oVec.size(), 0);
    QVector<double> antiIdeal(oVec.size());
    for(int i=0; i<oVec.size(); i++) {
        antiIdeal[i] = 3.0*(i+1);
    }
    double lb = 2.0/3.0;
    double ub = 1.0;


    UncertaintyKernel uk(oVec, lb, ub, ideal, antiIdeal);

    // Evaluate the uncertainty parameters
    double uniLB, uniLoc, uniUB, peakTend, peakLoc, dirPertRad;

    uniLB  = uk.proximity();
    uniLoc = skewedDecrease(uk.proximity(), 1.5);
    uniUB  = uniLB + (1-uniLoc) * (1.0-uniLB);

    peakTend = uk.proximity();
    peakLoc  = uk.symmetry();

    double distanceNorm = 2.0;

    dirPertRad = 0.04*lowOnValue(uk.oComponent(0), 0.45, 0.3);

    // Create the CODeM distribution
    MergedDistribution* d = MergedDistribution*::create();
    d->appendDistribution(UniformDistribution*(
                              new UniformDistribution(uniLB, uniUB)), 0.5);
    d->appendDistribution(PeakDistribution*(
                              new PeakDistribution(peakTend, peakLoc)), 0.5);

    CODeMDistribution cd(d, oVec, lb, ub, ideal, antiIdeal, dirPertRad, distanceNorm);

    // Sample the distribution
    QVector<QVector<double> > samples;
    for(int i=0; i<nSamp; i++) {
        samples.append(cd.sampleDistribution());
    }
    return samples;
}

QVector<double> CODeM4(QVector<double> iVec, int k, int nObj)
{
    // Evaluate the decision vector
    QVector<double> oVec =
            QVector<double>::fromStdVector(WFG6(iVec.toStdVector(), k, nObj));

    return CODeM4Perturb(oVec)[0];
}

QVector<QVector<double> > CODeM4(QVector<double> iVec,
                                int k, int nObj, int nSamp)
{
    // Evaluate the decision vector
    QVector<double> oVec =
            QVector<double>::fromStdVector(WFG6(iVec.toStdVector(), k, nObj));

    return CODeM4Perturb(oVec, nSamp);
}

QVector<QVector<double> > CODeM4Perturb(QVector<double> oVec, int nSamp)
{
    // Set the uncertainty kernel
    QVector<double> ideal(oVec.size(), 0.0);
    QVector<double> antiIdeal(oVec.size());
    for(int i=0; i<oVec.size(); i++) {
        antiIdeal[i] = 3.0*(i+1);
    }
    double lb = 2.0/3.0;
    double ub = 1.0;

    UncertaintyKernel uk(oVec, lb, ub, ideal, antiIdeal);

    // Evaluate the uncertainty parameters
    double peakTend, peakLoc, dirPertRad;

    peakTend = uk.proximity()+0.1;
    peakLoc  = 0.8;

    double distanceNorm = 2.0;

    dirPertRad = 0.2*linearDecrease(uk.symmetry())+0.01;

    // Create the CODeM distribution
    PeakDistribution* d = PeakDistribution*(
                new PeakDistribution(peakTend, peakLoc));

    CODeMDistribution cd(d, oVec, lb, ub, ideal, antiIdeal, dirPertRad, distanceNorm);

    // Sample the distribution
    QVector<QVector<double> > samples;
    for(int i=0; i<nSamp; i++) {
        samples.append(cd.sampleDistribution());
    }
    return samples;
}

QVector<double> CODeM5(QVector<double> iVec, int k, int nObj)
{
    // Evaluate the decision vector
    QVector<double> oVec =
            QVector<double>::fromStdVector(WFG8(iVec.toStdVector(), k, nObj));

    return CODeM5Perturb(iVec, oVec)[0];
}

QVector<QVector<double> > CODeM5(QVector<double> iVec,
                                int k, int nObj, int nSamp)
{
    // Evaluate the decision vector
    QVector<double> oVec =
            QVector<double>::fromStdVector(WFG8(iVec.toStdVector(), k, nObj));

    return CODeM5Perturb(iVec, oVec, nSamp);
}

QVector<QVector<double> > CODeM5Perturb(QVector<double> iVec,
                                       QVector<double> oVec,
                                       int nSamp)
{
    // Set the uncertainty kernel
    QVector<double> ideal(oVec.size(), 0.0);
    QVector<double> antiIdeal(oVec.size());
    for(int i=0; i<oVec.size(); i++) {
        antiIdeal[i] = 4.0*(i+1);
    }
    double lb = 2.0/4.0;
    double ub = 1.0;

    BoxConstraintsData* box = createBoxConstraints(5, iVec.size());

    UncertaintyKernel uk(iVec, oVec, box, lb, ub, ideal, antiIdeal);

    // Evaluate the uncertainty parameters
    double uniLB, uniUB, uniLoc, dirPertRad;

    uniLB  = uk.proximity();
    uniLoc = linearDecrease(uk.dComponent(0));
    uniUB  = uniLB + (1.0-uniLoc) * (1.0-uniLB);

    double distanceNorm = 2.0;

    dirPertRad = 0.1 * uk.dComponent(0);

    // Create the CODeM distribution
    UniformDistribution* d = UniformDistribution*(
                new UniformDistribution(uniLB, uniUB));

    CODeMDistribution cd(d, oVec, lb, ub, ideal, antiIdeal, dirPertRad, distanceNorm);

    // Sample the distribution
    QVector<QVector<double> > samples;
    for(int i=0; i<nSamp; i++) {
        samples.append(cd.sampleDistribution());
    }
    return samples;
}

QVector<double> CODeM6(QVector<double> iVec, int nObj)
{
    // Evaluate the decision vector
    QVector<double> oVec = DTLZ::DTLZ1(iVec, nObj);

    return CODeM6Perturb(iVec, oVec)[0];
}

QVector<QVector<double> > CODeM6(QVector<double> iVec,
                                int nObj, int nSamp)
{
    // Evaluate the decision vector
    QVector<double> oVec = DTLZ::DTLZ1(iVec, nObj);

    return CODeM6Perturb(iVec, oVec, nSamp);
}

QVector<QVector<double> > CODeM6Perturb(QVector<double> iVec, QVector<double> oVec,
                                       int nSamp)
{
    // Set the uncertainty kernel
    QVector<double> ideal(oVec.size(), 0.0);

    // DTLZ1 is modified so the 100 scale of the distance function
    // is not included
    double maxVal = 1.125 * iVec.size();
    QVector<double> antiIdeal(oVec.size(), maxVal);

    QVector<double> normVec(oVec);
    toUnitVec(normVec, 2.0);
    double sFactor = magnitudeAndDirectionP(normVec, 1);

    double ub = 1.0 / sFactor;
    double lb = 0.5 / maxVal / sFactor;

    UncertaintyKernel uk(oVec, lb, ub, ideal, antiIdeal);

    // Evaluate the uncertainty parameters
    double uniLB, uniUB, uniLoc, dirPertRad;

    uniLB  = uk.proximity();
    uniLoc = linearDecrease(uk.proximity()*uk.symmetry());
    uniUB  = uniLB + (1.0-uniLoc) * (1.0-uniLB);

    double distanceNorm = 1.0;

    dirPertRad = 0.2 * uk.oComponent(0);

    // Create the CODeM distribution
    UniformDistribution* d = UniformDistribution*(
                new UniformDistribution(uniLB, uniUB));

    CODeMDistribution cd(d, oVec, lb, ub, ideal, antiIdeal, dirPertRad, distanceNorm);

    // Sample the distribution
    QVector<QVector<double> > samples;
    for(int i=0; i<nSamp; i++) {
        samples.append(cd.sampleDistribution());
    }
    return samples;
}

QVector<double> deterministicOVec(int prob, QVector<double> iVec, int nObj, int k)
{
    QVector<double> oVec;

    switch(prob) {
    case 1: case 2: case 3:
        oVec = QVector<double>::fromStdVector(WFG4(iVec.toStdVector(), k, nObj));
        break;

    case 4:
        oVec = QVector<double>::fromStdVector(WFG6(iVec.toStdVector(), k, nObj));
        break;

    case 5:
        oVec = QVector<double>::fromStdVector(WFG8(iVec.toStdVector(), k, nObj));
        break;

    case 6:
        oVec = DTLZ::DTLZ1(iVec, nObj);
        break;
    default:
        break;
    }

    return oVec;
}

BoxConstraintsData* createBoxConstraints(int prob, int nVar)
{
    QVector<IElement> lowerBounds;
    QVector<IElement> upperBounds;

    switch(prob) {

    case 1: case 2: case 3: case 4: case 5:
        for(int i = 0; i < nVar; i++) {
            lowerBounds.append(IElement(RealType, 0.0));
            upperBounds.append(IElement(RealType, 2.0*(i+1.0)));
        }
        break;

    case 6:
        for(int i = 0; i < nVar; i++) {
            lowerBounds.append(IElement(RealType, 0.0));
            upperBounds.append(IElement(RealType, 1.0));
        }
        break;

    default:
        for(int i = 0; i < nVar; i++) {
            lowerBounds.append(IElement(RealType, 0.0));
            upperBounds.append(IElement(RealType, 1.0));
        }
        break;
    }

    BoxConstraintsData* box = BoxConstraintsData*::create();
    box->defineLowerBounds(lowerBounds);
    box->defineUpperBounds(upperBounds);
    return box;
}

} // namespace CODeM
