//=============================================================================================================
/**
 * @file     baandpowerettingsview.cpp
 * @author   Johannes Vorwerk <johannes.vorwerk@umit-tirol.at>
 *           Lorenz Esch <lesch@mgh.harvard.edu>;
 *           Christoph Dinh <chdinh@nmr.mgh.harvard.edu>
 * @version  dev
 * @date     May, 2020
 *
 * @section  LICENSE
 *
 * Copyright (C) 2020, Johannes Vorwerk, Lorenz Esch, Christoph Dinh. All rights reserved.
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
 * @brief    Definition of the BandPowerSettingsView Class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "bandpowersettingsview.h"
#include "ui_bandpowersettingsview.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QDoubleSpinBox>
#include <QLabel>
#include <QGridLayout>
#include <QSlider>
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

BandPowerSettingsView::BandPowerSettingsView(const QString& sSettingsPath, double dSampFreq, double dMin, double dMax,
                                             const QString& sSpectrumMethod, int iIntervallLength, QWidget *parent, Qt::WindowFlags f)
: QWidget(parent, f)
, ui(new Ui::BandPowerSettingsViewWidget)
, m_sSettingsPath(sSettingsPath)
, m_dSampFreq(dSampFreq)
, m_dMin(dMin)
, m_dMax(dMax)
, m_sSpectrumMethod(sSpectrumMethod)
, m_iIntervallLength(iIntervallLength)
, m_iDetrend(0)
{
    ui->setupUi(this);

    this->setWindowTitle("BandPower Settings");
    this->setMinimumWidth(330);
    this->setMaximumWidth(330);

    //loadSettings(m_sSettingsPath);

    if(m_dMax <= m_dMin)
    {
        if(m_dMin < 2)
            m_dMax = m_dMin + 1.0;
        else if(m_dMax > (m_dSampFreq/2.0 - 1.0))
            m_dMin = m_dMax - 1.0;
        else {
            m_dMin -= 0.5;
            m_dMax += 0.5;
        }
    }

    ui->doubleSpinBoxMin->setMinimum(1);
    ui->doubleSpinBoxMin->setMaximum(static_cast<int>(m_dSampFreq/2.0) - 1);
    ui->doubleSpinBoxMin->setMaximumWidth(100);
    ui->doubleSpinBoxMin->setSingleStep(1);
    ui->doubleSpinBoxMin->setDecimals(1);
    ui->doubleSpinBoxMin->setSuffix(" Hz");
    ui->doubleSpinBoxMin->setValue(m_dMin);

    ui->horizontalSliderMin->setMinimum(1);
    ui->horizontalSliderMin->setMaximum(static_cast<int>(m_dSampFreq/2.0) - 1);
    ui->horizontalSliderMin->setSingleStep(1);
    ui->horizontalSliderMin->setPageStep(10);
    ui->horizontalSliderMin->setValue(m_dMin);

    ui->doubleSpinBoxMax->setMinimum(2);
    ui->doubleSpinBoxMax->setMaximum(static_cast<int>(m_dSampFreq/2.0));
    ui->doubleSpinBoxMax->setMaximumWidth(100);
    ui->doubleSpinBoxMax->setSingleStep(1);
    ui->doubleSpinBoxMax->setDecimals(1);
    ui->doubleSpinBoxMax->setSuffix(" Hz");
    ui->doubleSpinBoxMax->setValue(m_dMax);

    ui->horizontalSliderMax->setMinimum(2);
    ui->horizontalSliderMax->setMaximum(static_cast<int>(m_dSampFreq/2.0));
    ui->horizontalSliderMax->setSingleStep(1);
    ui->horizontalSliderMax->setPageStep(10);
    ui->horizontalSliderMax->setValue(m_dMax);

    connect(ui->doubleSpinBoxMin,&QDoubleSpinBox::editingFinished,
            this,&BandPowerSettingsView::onUpdateSpinBoxMinScaling);
    connect(ui->horizontalSliderMin,&QSlider::sliderReleased,
            this,&BandPowerSettingsView::onUpdateSliderMinScaling);
    connect(ui->doubleSpinBoxMax,&QDoubleSpinBox::editingFinished,
            this,&BandPowerSettingsView::onUpdateSpinBoxMaxScaling);
    connect(ui->horizontalSliderMax,&QSlider::sliderReleased,
            this,&BandPowerSettingsView::onUpdateSliderMaxScaling);

    ui->radioButtonAR->setChecked(true);
    ui->radioButtonFFT->setChecked(false);
    ui->buttonGroupMethod->setId(ui->radioButtonAR,0);
    ui->buttonGroupMethod->setId(ui->radioButtonFFT,1);

    connect(ui->buttonGroupMethod,static_cast<void (QButtonGroup::*)(int)>(&QButtonGroup::buttonClicked),this,&BandPowerSettingsView::onClickedButtonMethod);

    ui->spinBoxIntervallLength->setValue(m_iIntervallLength);
    connect(ui->spinBoxIntervallLength,&QDoubleSpinBox::editingFinished,
            this,&BandPowerSettingsView::onUpdateSpinBoxIntervallLength);

    ui->radioButtonNoDetrend->setChecked(true);
    ui->radioButtonMean->setChecked(false);
    ui->radioButtonLinear->setChecked(false);
    ui->buttonGroupMethod->setId(ui->radioButtonNoDetrend,0);
    ui->buttonGroupMethod->setId(ui->radioButtonMean,1);
    ui->buttonGroupMethod->setId(ui->radioButtonLinear,2);

    connect(ui->buttonGroupDetrend,static_cast<void (QButtonGroup::*)(int)>(&QButtonGroup::buttonClicked),this,&BandPowerSettingsView::onClickedButtonDetrend);

    switch(m_iDetrend)
    {
    case 0:
        ui->radioButtonNoDetrend->setChecked(true);
        ui->radioButtonMean->setChecked(false);
        ui->radioButtonLinear->setChecked(false);
        break;
    case 1:
        ui->radioButtonNoDetrend->setChecked(false);
        ui->radioButtonMean->setChecked(true);
        ui->radioButtonLinear->setChecked(false);
        break;
    case 2:
        ui->radioButtonNoDetrend->setChecked(false);
        ui->radioButtonMean->setChecked(false);
        ui->radioButtonLinear->setChecked(true);
        break;
    }

}

//=============================================================================================================

BandPowerSettingsView::~BandPowerSettingsView()
{
    saveSettings(m_sSettingsPath);

    delete ui;
}

//=============================================================================================================

void BandPowerSettingsView::saveSettings(const QString& settingsPath)
{
    if(settingsPath.isEmpty()) {
        return;
    }

    // Store Settings
    QSettings settings;

    settings.setValue(settingsPath + QString("/min"), m_dMin);
    settings.setValue(settingsPath + QString("/max"), m_dMax);
    settings.setValue(settingsPath + QString("/spectrumMethod"), m_sSpectrumMethod);
    settings.setValue(settingsPath + QString("/intervallLength"), m_iIntervallLength);
    settings.setValue(settingsPath + QString("/detrend"), m_iDetrend);

}

//=============================================================================================================

void BandPowerSettingsView::loadSettings(const QString& settingsPath)
{
    if(settingsPath.isEmpty()) {
        return;
    }

    // Load Settings
    QSettings settings;

    m_dMin = settings.value(settingsPath + QString("/min"), m_dMin).toDouble();
    m_dMax = settings.value(settingsPath + QString("/max"), m_dMax).toDouble();
    m_sSpectrumMethod = settings.value(settingsPath + QString("/spectrumMethod"), "AR").toString();
    m_iIntervallLength = settings.value(settingsPath + QString("/intervallLength"), m_iIntervallLength).toInt();
    m_iDetrend = settings.value(settingsPath + QString("/detrend"), m_iDetrend).toInt();

}

//=============================================================================================================

void BandPowerSettingsView::onUpdateSpinBoxMinScaling()
{
    double value = ui->doubleSpinBoxMin->value();

    if(value == m_dMin)
        return;

    m_dMin = value;

    ui->horizontalSliderMin->blockSignals(true);
    ui->horizontalSliderMin->setValue(static_cast<int>(m_dMin));
    ui->horizontalSliderMin->blockSignals(false);

    if(m_dMin >= m_dMax) {
        m_dMax = m_dMin + 1;

        ui->doubleSpinBoxMax->blockSignals(true);
        ui->doubleSpinBoxMax->setValue(m_dMax);
        ui->doubleSpinBoxMax->blockSignals(false);

        ui->horizontalSliderMax->blockSignals(true);
        ui->horizontalSliderMax->setValue(static_cast<int>(m_dMax));
        ui->horizontalSliderMax->blockSignals(false);
    }

    emit changeMinMax(m_dMin,m_dMax);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void BandPowerSettingsView::onUpdateSpinBoxMaxScaling()
{
    double value = ui->doubleSpinBoxMax->value();

    if(value == m_dMax)
        return;

    m_dMax = value;

    ui->horizontalSliderMax->blockSignals(true);
    ui->horizontalSliderMax->setValue(static_cast<int>(m_dMax));
    ui->horizontalSliderMax->blockSignals(false);

    if(m_dMax <= m_dMin) {
        m_dMin = m_dMax - 1;

        ui->doubleSpinBoxMin->blockSignals(true);
        ui->doubleSpinBoxMin->setValue(m_dMin);
        ui->doubleSpinBoxMin->blockSignals(false);

        ui->horizontalSliderMin->blockSignals(true);
        ui->horizontalSliderMin->setValue(static_cast<int>(m_dMin));
        ui->horizontalSliderMin->blockSignals(false);
    }

    emit changeMinMax(m_dMin,m_dMax);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void BandPowerSettingsView::onUpdateSliderMinScaling()
{
    double value = ui->horizontalSliderMin->value();

    if(value == m_dMin)
        return;

    m_dMin = value;

    ui->doubleSpinBoxMin->blockSignals(true);
    ui->doubleSpinBoxMin->setValue(m_dMin);
    ui->doubleSpinBoxMin->blockSignals(false);

    if(m_dMin >= m_dMax) {
        m_dMax = m_dMin + 1;

        ui->doubleSpinBoxMax->blockSignals(true);
        ui->doubleSpinBoxMax->setValue(m_dMax);
        ui->doubleSpinBoxMax->blockSignals(false);

        ui->horizontalSliderMax->blockSignals(true);
        ui->horizontalSliderMax->setValue(static_cast<int>(m_dMax));
        ui->horizontalSliderMax->blockSignals(false);
    }

    emit changeMinMax(m_dMin,m_dMax);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void BandPowerSettingsView::onUpdateSliderMaxScaling()
{
    double value = ui->horizontalSliderMax->value();

    if(value == m_dMax)
        return;

    m_dMax = value;

    ui->doubleSpinBoxMax->blockSignals(true);
    ui->doubleSpinBoxMax->setValue(m_dMax);
    ui->doubleSpinBoxMax->blockSignals(false);

    if(m_dMax <= m_dMin) {
        m_dMin = m_dMax - 1;

        ui->doubleSpinBoxMin->blockSignals(true);
        ui->doubleSpinBoxMin->setValue(m_dMin);
        ui->doubleSpinBoxMin->blockSignals(false);

        ui->horizontalSliderMin->blockSignals(true);
        ui->horizontalSliderMin->setValue(static_cast<int>(m_dMin));
        ui->horizontalSliderMin->blockSignals(false);
    }

    emit changeMinMax(m_dMin,m_dMax);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void BandPowerSettingsView::onClickedButtonMethod(int value)
{
    if(value == 0)
        m_sSpectrumMethod = "AR";
    else
        m_sSpectrumMethod = "FFT";

    emit changeMethod(m_sSpectrumMethod);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void BandPowerSettingsView::onUpdateSpinBoxIntervallLength()
{
    double value = ui->spinBoxIntervallLength->value();

    if(value == m_iIntervallLength)
        return;

    m_iIntervallLength = value;

    emit changeIntervallLength(m_iIntervallLength);

    saveSettings(m_sSettingsPath);
}

void BandPowerSettingsView::onClickedButtonDetrend(int value)
{
    m_iDetrend = value;

    emit changeDetrend(m_iDetrend);

    saveSettings(m_sSettingsPath);
}
