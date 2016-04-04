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
#ifndef CODEM5_H
#define CODEM5_H

#include <tigon/Representation/Functions/IFunction.h>
#include <tigon/tigon_global.h>


// Qt Includes
#include <QVector>

namespace Tigon {
namespace Representation {

class LIGER_TIGON_EXPORT CODeM5 : public IFunction
{
    Q_OBJECT

public:
    CODeM5();
    CODeM5(const CODeM5& func);
    virtual ~CODeM5();

    void evaluate(QVector<IElementSPtr> inputs,
                          QVector<IElementSPtr> outputs);

    void defineNumDirVars(int k);
    int  numDirVars()      const;

private:
    void defineInputPrpts();
    void defineOutputPrpts();

    int m_nDirVar;
};

} // namespace Representation
} // namespace Tigon

#endif // CODEM5_H
