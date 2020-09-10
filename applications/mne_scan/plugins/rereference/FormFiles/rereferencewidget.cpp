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

#include "rereferencewidget.h"
#include "../ui_rereferencewidget.h"

#include "../rereference.h"

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

using namespace REREFERENCEPLUGIN;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

RereferenceWidget::RereferenceWidget(const QString& sSettingsPath,
                                 QWidget *parent)
: QWidget(parent)
, ui(new Ui::RereferenceWidget)
, m_sSettingsPath(sSettingsPath)
, m_bEnabled(true)
, m_iModality(0)
, m_iMethod(0)
{
    ui->setupUi(this);

    ui->radioButton_EEG->setChecked(true);
    ui->radioButton_MEG->setChecked(false);
    ui->radioButton_EMEG->setChecked(false);
    ui->buttonGroupModality->setId(ui->radioButton_EEG,0);
    ui->buttonGroupModality->setId(ui->radioButton_MEG,1);
    ui->buttonGroupModality->setId(ui->radioButton_EMEG,2);
    connect(ui->buttonGroupModality,static_cast<void (QButtonGroup::*)(int)>(&QButtonGroup::buttonClicked),this,&RereferenceWidget::onClickedButtonModality);

    ui->radioButton_CAR->setChecked(true);
    ui->radioButton_Selection->setChecked(false);
    ui->radioButton_File->setChecked(false);
    ui->buttonGroupMethod->setId(ui->radioButton_CAR,0);
    ui->buttonGroupMethod->setId(ui->radioButton_Selection,1);
    ui->buttonGroupMethod->setId(ui->radioButton_File,2);
    connect(ui->buttonGroupMethod,static_cast<void (QButtonGroup::*)(int)>(&QButtonGroup::buttonClicked),this,&RereferenceWidget::onClickedButtonMethod);

    ui->checkBox_Enabled->setChecked(m_bEnabled);
    connect(ui->checkBox_Enabled,&QCheckBox::stateChanged,this,&RereferenceWidget::onClickedCheckboxEnabled);

    loadSettings(m_sSettingsPath);
}

//=============================================================================================================

RereferenceWidget::~RereferenceWidget()
{
    saveSettings(m_sSettingsPath);

    delete ui;
}

//=============================================================================================================

void RereferenceWidget::saveSettings(const QString& settingsPath)
{
    if(settingsPath.isEmpty()) {
        return;
    }

    QSettings settings;

    settings.setValue(settingsPath + QString("/modality"), m_iModality);
    settings.setValue(settingsPath + QString("/method"), m_iMethod);
}

//=============================================================================================================

void RereferenceWidget::loadSettings(const QString& settingsPath)
{
    if(settingsPath.isEmpty()) {
        return;
    }

    QSettings settings;

    m_iModality = settings.value(settingsPath + QString("/modality"), m_iModality).toInt();
    m_iMethod = settings.value(settingsPath + QString("/method"), m_iMethod).toInt();
}

//=============================================================================================================

void RereferenceWidget::onClickedButtonModality(int value)
{
    m_iModality = value;

    emit changeModality(m_iModality);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void RereferenceWidget::onClickedButtonMethod(int value)
{
    m_iMethod = value;

    emit changeMethod(m_iMethod);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void RereferenceWidget::onClickedCheckboxEnabled(bool value)
{
    m_bEnabled = value;

    emit changeEnabled(m_bEnabled);
}

//=============================================================================================================

void RereferenceWidget::changeFile()
{
    QString path = QFileDialog::getOpenFileName(this,
                                                "Select EEG rereferencing file",
                                                "resources/mne_scan/plugins/rereferencing/loc_files",
                                                 tr("Simple text files (*.*)"));

    if(path==NULL){
        path = ui->lineEdit_File->text();
    }

    ui->lineEdit_File->setText(path);
    //m_pEEGoSports->m_sElcFilePath = m_pUi->m_qLineEdit_EEGCap->text();
}

