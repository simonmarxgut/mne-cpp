//=============================================================================================================
/**
 * @file     rtspectrumviewsettings.cpp
 * @author   Lorenz Esch <lesch@mgh.harvard.edu>
 * @since    0.1.0
 * @date     July, 2018
 *
 * @section  LICENSE
 *
 * Copyright (C) 2018, Lorenz Esch. All rights reserved.
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
 * @brief    Definition of the RtSpectrumViewSettings Class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "rtspectrumviewsettings.h"

#include "ui_rtspectrumviewsettings.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QColorDialog>
#include <QSettings>
#include <QDebug>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace DISPLIB;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

RtSpectrumViewSettings::RtSpectrumViewSettings(const QString &sSettingsPath, double dWindowLength, double dWindowStepsize,
                                         QWidget *parent,
                                         Qt::WindowFlags f)
: QWidget(parent, f)
, ui(new Ui::RtSpectrumViewSettingsWidget)
, m_sSettingsPath(sSettingsPath)
{
    ui->setupUi(this);

    this->setWindowTitle("Spectrum View Settings");
    this->setMinimumWidth(330);
    this->setMaximumWidth(330);

    ui->doubleSpinBox_Maximum->setDecimals(10);
    ui->doubleSpinBox_Minimum->setDecimals(10);

    ui->doubleSpinBox_windowSize->setSingleStep(dWindowStepsize);
    ui->doubleSpinBox_windowSize->setValue(dWindowLength);

    connect(ui->doubleSpinBox_windowSize, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &RtSpectrumViewSettings::onTimeWindowChanged);
    connect(ui->checkBox_Colormap, &QCheckBox::stateChanged, this, &RtSpectrumViewSettings::onFixColormapChanged);
    connect(ui->doubleSpinBox_Maximum, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &RtSpectrumViewSettings::onColormapMaxChanged);
    connect(ui->doubleSpinBox_Minimum, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &RtSpectrumViewSettings::onColormapMinChanged);

    //loadSettings(m_sSettingsPath);
}

//=============================================================================================================

RtSpectrumViewSettings::~RtSpectrumViewSettings()
{
    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void RtSpectrumViewSettings::setWindowSize(int windowSize)
{
    ui->doubleSpinBox_windowSize->setValue(windowSize);

    timeWindowChanged(windowSize);
}

//=============================================================================================================

void RtSpectrumViewSettings::setColormapMax(double value)
{
    if(!ui->checkBox_Colormap->isChecked())
    {
        bool oldState = ui->doubleSpinBox_Maximum->blockSignals(true);
        ui->doubleSpinBox_Maximum->setValue(value);
        ui->doubleSpinBox_Maximum->blockSignals(oldState);
    }
}

//=============================================================================================================

void RtSpectrumViewSettings::setColormapMin(double value)
{
    if(!ui->checkBox_Colormap->isChecked())
    {
        bool oldState = ui->doubleSpinBox_Minimum->blockSignals(true);
        ui->doubleSpinBox_Minimum->setValue(value);
        ui->doubleSpinBox_Minimum->blockSignals(oldState);
    }
}

//=============================================================================================================

int RtSpectrumViewSettings::getWindowSize()
{
    return ui->doubleSpinBox_windowSize->value();
}

//=============================================================================================================

void RtSpectrumViewSettings::saveSettings(const QString& settingsPath)
{
    if(settingsPath.isEmpty()) {
        return;
    }

    QSettings settings;

    settings.setValue(settingsPath + QString("/viewWindowSize"), getWindowSize());
}

//=============================================================================================================

void RtSpectrumViewSettings::loadSettings(const QString& settingsPath)
{
    if(settingsPath.isEmpty()) {
        return;
    }

    QSettings settings;
    setWindowSize(settings.value(settingsPath + QString("/viewWindowSize"), 10).toDouble());
}

//=============================================================================================================

void RtSpectrumViewSettings::onTimeWindowChanged(double value)
{
    emit timeWindowChanged(value);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void RtSpectrumViewSettings::onMakeScreenshot()
{
    emit makeScreenshot(ui->m_comboBox_imageType->currentText());
}

//=============================================================================================================

void RtSpectrumViewSettings::onFixColormapChanged(bool value)
{
    emit fixColormapChanged(value);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void RtSpectrumViewSettings::onColormapMaxChanged(double value)
{
    emit colormapMaxChanged(value);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void RtSpectrumViewSettings::onColormapMinChanged(double value)
{
    emit colormapMinChanged(value);

    saveSettings(m_sSettingsPath);
}
