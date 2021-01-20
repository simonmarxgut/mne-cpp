//=============================================================================================================
/**
 * @file     rereferencewidget.cpp
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
 * @brief    Definition of the RereferenceWidget class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "normalizewidget.h"
#include "../ui_normalizewidget.h"
#include "../normalize.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QSettings>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace NORMALIZEPLUGIN;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

NormalizeWidget::NormalizeWidget(const QString& sSettingsPath, double dDataSampFreq,
                                 int iDataBlockLength, bool bEnabled, int iMethod,
                                 int iModality, double dIntervallLength, QWidget *parent)
: QWidget(parent)
, ui(new Ui::NormalizeWidget)
, m_sSettingsPath(sSettingsPath)
, m_iDataBlockLength(iDataBlockLength)
, m_bEnabled(bEnabled)
, m_iModality(iModality)
, m_iMethod(iMethod)
, m_dIntervallLength(dIntervallLength)
{
    ui->setupUi(this);

    if(!(m_dIntervallLength > 0)) {
        m_dIntervallLength = 1.0;
        qWarning() << "Intervall length has to be larger than zero. Corrected to " << m_dIntervallLength;
    }

    if(m_iMethod == 0) {
        ui->radioButton_None->setChecked(true);
        ui->radioButton_Demean->setChecked(false);
        ui->radioButton_DemeanVar->setChecked(false);
    } else if(m_iMethod == 1) {
        ui->radioButton_None->setChecked(false);
        ui->radioButton_Demean->setChecked(true);
        ui->radioButton_DemeanVar->setChecked(false);
    } else {
        ui->radioButton_None->setChecked(false);
        ui->radioButton_Demean->setChecked(false);
        ui->radioButton_DemeanVar->setChecked(true);
    }
    ui->buttonGroupModality->setId(ui->radioButton_None,0);
    ui->buttonGroupModality->setId(ui->radioButton_Demean,1);
    ui->buttonGroupModality->setId(ui->radioButton_DemeanVar,2);
    connect(ui->buttonGroupModality,static_cast<void (QButtonGroup::*)(int)>(&QButtonGroup::buttonClicked),this,&NormalizeWidget::onClickedButtonModality);

    if((m_iModality > 0) || !(m_iDataBlockLength > 0)) {
        ui->radioButton_DataBlocks->setChecked(false);
        ui->radioButton_Seconds->setChecked(true);
    } else {
        ui->radioButton_DataBlocks->setChecked(true);
        ui->radioButton_Seconds->setChecked(false);
    }
    ui->buttonGroupMethod->setId(ui->radioButton_DataBlocks,0);
    ui->buttonGroupMethod->setId(ui->radioButton_Seconds,1);
    connect(ui->buttonGroupMethod,static_cast<void (QButtonGroup::*)(int)>(&QButtonGroup::buttonClicked),this,&NormalizeWidget::onClickedButtonMethod);

    if(m_iDataBlockLength > 0) {
        ui->doubleSpinBox_Seconds->setSingleStep(double(iDataBlockLength)/dDataSampFreq);
        if(m_iModality == 0) {
            m_dIntervallLength = std::round(m_dIntervallLength * dDataSampFreq / double(iDataBlockLength));
            ui->spinBox_DataBlocks->setValue(m_dIntervallLength);
            ui->doubleSpinBox_Seconds->setValue(ui->spinBox_DataBlocks->value()*iDataBlockLength/dDataSampFreq);
        } else {
            ui->spinBox_DataBlocks->setValue(std::round(m_dIntervallLength * dDataSampFreq / double(iDataBlockLength)));
            ui->doubleSpinBox_Seconds->setValue(m_dIntervallLength);
        }
        connect(ui->spinBox_DataBlocks,QOverload<int>::of(&QSpinBox::valueChanged),this,&NormalizeWidget::onChangeIntervallLengthDataBlocks);

    } else {
        ui->doubleSpinBox_Seconds->setSingleStep(0.1);
        ui->spinBox_DataBlocks->setEnabled(false);
        ui->doubleSpinBox_Seconds->setValue(1.0);
        ui->spinBox_DataBlocks->setValue(0);
        ui->radioButton_DataBlocks->setEnabled(false);
        ui->radioButton_Seconds->setEnabled(false);
    }

    connect(ui->doubleSpinBox_Seconds,QOverload<double>::of(&QDoubleSpinBox::valueChanged),this,&NormalizeWidget::onChangeIntervallLengthSeconds);

    ui->checkBox_Enabled->setChecked(m_bEnabled);
    connect(ui->checkBox_Enabled,&QCheckBox::stateChanged,this,&NormalizeWidget::onClickedCheckboxEnabled);

    loadSettings(m_sSettingsPath);
}

//=============================================================================================================

NormalizeWidget::~NormalizeWidget()
{
    saveSettings(m_sSettingsPath);

    delete ui;
}

//=============================================================================================================

void NormalizeWidget::saveSettings(const QString& settingsPath)
{
    if(settingsPath.isEmpty()) {
        return;
    }

    QSettings settings;

    settings.setValue(settingsPath + QString("/modality"), m_iModality);
    settings.setValue(settingsPath + QString("/method"), m_iMethod);
}

//=============================================================================================================

void NormalizeWidget::loadSettings(const QString& settingsPath)
{
    if(settingsPath.isEmpty()) {
        return;
    }

    QSettings settings;

    m_iModality = settings.value(settingsPath + QString("/modality"), m_iModality).toInt();
    m_iMethod = settings.value(settingsPath + QString("/method"), m_iMethod).toInt();
}

//=============================================================================================================

void NormalizeWidget::onClickedButtonModality(int value)
{
    m_iModality = value;

    emit changeModality(m_iModality);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void NormalizeWidget::onClickedButtonMethod(int value)
{
    m_iMethod = value;

    emit changeMethod(m_iMethod);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void NormalizeWidget::onClickedCheckboxEnabled(bool value)
{
    m_bEnabled = value;

    emit changeEnabled(m_bEnabled);
}

//=============================================================================================================

void NormalizeWidget::onChangeIntervallLengthDataBlocks(int value)
{
    m_dIntervallLength = value * m_iDataBlockLength / m_dDataSampFreq;

    ui->doubleSpinBox_Seconds->setValue(m_dIntervallLength);

    emit changeIntervallLength(m_dIntervallLength);
}

//=============================================================================================================

void NormalizeWidget::onChangeIntervallLengthSeconds(double value)
{
    m_dIntervallLength = value;

    if(m_iDataBlockLength > 0)
        ui->spinBox_DataBlocks->setValue(std::round(m_dIntervallLength * m_dDataSampFreq / m_iDataBlockLength));

    emit changeIntervallLength(m_dIntervallLength);
}
