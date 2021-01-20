//=============================================================================================================
/**
 * @file     arsettingsview.cpp
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
 * @brief    Definition of the ARSettingsView Class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "arsettingsview.h"
#include "ui_arsettingsview.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QSpinBox>
#include <QLabel>
#include <QGridLayout>
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

ARSettingsView::ARSettingsView(const QString& sSettingsPath, bool bChangeSamplingPoints, QWidget *parent, Qt::WindowFlags f)
    : QWidget(parent, f)
    , m_bChangeSamplingPoints(bChangeSamplingPoints)
    , ui(new Ui::ARSettingsViewWidget)
    , m_sSettingsPath(sSettingsPath)
    , m_iOrder(15)
    , m_iNumEvaluationPoints(100)
    , m_iSamplingPoints(1)
{
    ui->setupUi(this);

    this->setWindowTitle("AR Settings");
    this->setMinimumWidth(330);
    this->setMaximumWidth(330);

    loadSettings(m_sSettingsPath);

    ui->spinBoxAROrder->setMinimum(1);
    ui->spinBoxAROrder->setMaximum(999);
    ui->spinBoxAROrder->setMaximumWidth(100);
    ui->spinBoxAROrder->setSingleStep(1);
    ui->spinBoxAROrder->setValue(m_iOrder);

    connect(ui->spinBoxAROrder,&QSpinBox::editingFinished,this,&ARSettingsView::onUpdateSpinBoxAROrder);

    ui->spinBoxNumEvaluationPoints->setMinimum(1);
    ui->spinBoxNumEvaluationPoints->setMaximum(999);
    ui->spinBoxNumEvaluationPoints->setMaximumWidth(100);
    ui->spinBoxNumEvaluationPoints->setSingleStep(1);
    ui->spinBoxNumEvaluationPoints->setValue(m_iOrder);

    connect(ui->spinBoxNumEvaluationPoints,&QSpinBox::editingFinished,this,&ARSettingsView::onUpdateSpinBoxNumEvaluationPoints);

    if(m_bChangeSamplingPoints){
        //QLabel* labelSamplingPoints = new QLabel("Sampling Points");
        //labelSamplingPoints->setObjectName(QString::fromUtf8("label_3"));
        //ui->gridLayout->addWidget(labelSamplingPoints,2,0,1,1);

        ui->label_3->setText("Sampling Points");

        QSpinBox* spinBoxSamplingPoints = new QSpinBox;
        spinBoxSamplingPoints->setObjectName(QString::fromUtf8("spinBoxSamplingPoints"));
        ui->gridLayout->addWidget(spinBoxSamplingPoints,2,1,1,1);
        spinBoxSamplingPoints->setMinimum(1);
        spinBoxSamplingPoints->setMaximum(99);
        spinBoxSamplingPoints->setMaximumWidth(100);
        spinBoxSamplingPoints->setSingleStep(1);
        spinBoxSamplingPoints->setValue(m_iSamplingPoints);
        m_qSpinBoxSamplingPoints = spinBoxSamplingPoints;

        connect(spinBoxSamplingPoints,&QSpinBox::editingFinished,this,&ARSettingsView::onUpdateSpinBoxSamplingPoints);
    }

}

//=============================================================================================================

ARSettingsView::~ARSettingsView()
{
    saveSettings(m_sSettingsPath);

    //delete m_qSpinBoxSamplingPoints;
    delete ui;
}

//=============================================================================================================

void ARSettingsView::emitSignals()
{
    emit changeAROrder(m_iOrder);
    emit changeARNumEvaluationPoints(m_iNumEvaluationPoints);
}

//=============================================================================================================

void ARSettingsView::saveSettings(const QString& settingsPath)
{
    if(settingsPath.isEmpty()) {
        return;
    }

    // Store Settings
    QSettings settings;

    settings.setValue(settingsPath + QString("/order"), m_iOrder);
    settings.setValue(settingsPath + QString("/numEvaluationsPoints"), m_iNumEvaluationPoints);

}

//=============================================================================================================

void ARSettingsView::loadSettings(const QString& settingsPath)
{
    if(settingsPath.isEmpty()) {
        return;
    }

    // Load Settings
    QSettings settings;

    m_iOrder = settings.value(settingsPath + QString("/order"), m_iOrder).toInt();
    m_iNumEvaluationPoints = settings.value(settingsPath + QString("/numEvaluationPoints"), m_iNumEvaluationPoints).toInt();
}

//=============================================================================================================

void ARSettingsView::onUpdateSpinBoxAROrder()
{
    double value = ui->spinBoxAROrder->value();

    if(value == m_iOrder)
        return;

    m_iOrder = value;

    emit changeAROrder(m_iOrder);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void ARSettingsView::onUpdateSpinBoxNumEvaluationPoints()
{
    double value = ui->spinBoxNumEvaluationPoints->value();

    if(value == m_iNumEvaluationPoints)
        return;

    m_iNumEvaluationPoints = value;

    emit changeARNumEvaluationPoints(m_iNumEvaluationPoints);

    saveSettings(m_sSettingsPath);
}

//=============================================================================================================

void ARSettingsView::onUpdateSpinBoxSamplingPoints()
{
    double value = m_qSpinBoxSamplingPoints->value();

    if(value == m_iSamplingPoints)
        return;

    m_iSamplingPoints = value;

    emit changeARNumEvaluationPoints(m_iSamplingPoints);

    saveSettings(m_sSettingsPath);
}
