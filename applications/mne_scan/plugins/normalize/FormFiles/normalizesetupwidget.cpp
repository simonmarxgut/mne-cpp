//=============================================================================================================
/**
 * @file     normalizesetupwidget.cpp
 * @author   Christoph Dinh <chdinh@nmr.mgh.harvard.edu>;
 *           Lorenz Esch <lesch@mgh.harvard.edu>;
 *           Viktor Klueber <Viktor.Klueber@tu-ilmenau.de>
 * @since    0.1.0
 * @date     February, 2013
 *
 * @section  LICENSE
 *
 * Copyright (C) 2013, Christoph Dinh, Lorenz Esch, Viktor Klueber. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that
 * the following conditions are met:
 *     * Redistributions of source code must retain the above copyright notice, this list of conditions and the
 *       following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
 *       the following disclaimer in the documentation and/or other materials provided with the distribution.
 *     * Neither the name of MNE-CPP authors nor the names of its contributors may be used
 *       to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * @brief    Definition of the NormalizeSetupWidget class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "normalizesetupwidget.h"
#include "../ui_normalizesetup.h"

#include "../normalize.h"
#include "normalizewidget.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QDebug>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace NORMALIZEPLUGIN;
using namespace Ui;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

NormalizeSetupWidget::NormalizeSetupWidget(Normalize* pNormalize, const QString& sSettingsPath, QWidget *parent)
: QWidget(parent)
, m_pNormalize(pNormalize)
{
    ui = new NormalizeSetupWidgetClass();
    ui->setupUi(this);

    NormalizeWidget* pNormalizeWidget = new NormalizeWidget(sSettingsPath,-1,-1,true);

    QVBoxLayout* settingsLayout = new QVBoxLayout;
    settingsLayout->addWidget(pNormalizeWidget);

    ui->m_qGroupBox_Options->setLayout(settingsLayout);

    connect(pNormalizeWidget,&NormalizeWidget::changeEnabled,m_pNormalize,&Normalize::changeEnabled);
    connect(pNormalizeWidget,&NormalizeWidget::changeModality,m_pNormalize,&Normalize::changeModality);
    connect(pNormalizeWidget,&NormalizeWidget::changeMethod,m_pNormalize,&Normalize::changeMethod);
    connect(pNormalizeWidget,&NormalizeWidget::changeIntervallLength,m_pNormalize,&Normalize::changeIntervallLength);
}

//=============================================================================================================

NormalizeSetupWidget::~NormalizeSetupWidget()
{
}
