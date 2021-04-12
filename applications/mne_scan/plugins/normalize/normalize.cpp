//=============================================================================================================
/**
 * @file     normalize.cpp
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
 * @brief    Definition of the Normalize class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "normalize.h"

#include "FormFiles/normalizesetupwidget.h"
#include "FormFiles/normalizewidget.h"

#include <fiff/fiff.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <Eigen/Core>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace NORMALIZEPLUGIN;
using namespace SCSHAREDLIB;
using namespace SCMEASLIB;
using namespace FIFFLIB;
using namespace UTILSLIB;
using namespace Eigen;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

Normalize::Normalize()
: m_pNormalizeBuffer(CircularBuffer<Eigen::MatrixXd>::SPtr::create(64))
, m_iNChannels(-1)
, m_dDataSampFreq(-1)
, m_iDataBlockLength(-1)
, m_iNormalizeMethod(0)
, m_iNormalizeModality(0)
, m_bNormalizeEnabled(true)
, m_dIntervallLength(1)
{
}

//=============================================================================================================

Normalize::~Normalize()
{
    if(this->isRunning()) {
        stop();
    }
}


//=============================================================================================================

QSharedPointer<AbstractPlugin> Normalize::clone() const
{
    QSharedPointer<Normalize> pNormalizeClone(new Normalize);
    return pNormalizeClone;
}

//=============================================================================================================

void Normalize::init()
{
    QMutexLocker locker(&m_qMutex);

    // Input
    m_pNormalizeInput = PluginInputData<RealTimeMultiSampleArray>::create(this, "NormalizeIn", "Normalize input data");
    connect(m_pNormalizeInput.data(), &PluginInputConnector::notify,
            this, &Normalize::update, Qt::DirectConnection);
    m_inputConnectors.append(m_pNormalizeInput);

    // Output - Uncomment this if you don't want to send processed data (in form of a matrix) to other plugins.
    // Also, this output stream will generate an online display in your plugin
    m_pNormalizeOutput = PluginOutputData<RealTimeMultiSampleArray>::create(this, "ReferenceOut", "Reference output data");
    m_pNormalizeOutput->measurementData()->setName(this->getName());
    m_outputConnectors.append(m_pNormalizeOutput);

    //Delete Buffer - will be initailzed with first incoming data
    if(!m_pNormalizeBuffer.isNull())
        m_pNormalizeBuffer = CircularBuffer<Eigen::MatrixXd>::SPtr();
}

//=============================================================================================================

void Normalize::unload()
{
}

//=============================================================================================================

bool Normalize::start()
{
    //Start thread
    QThread::start();

    return true;
}

//=============================================================================================================

bool Normalize::stop()
{
    requestInterruption();
    wait(500);

    // Clear all data in the buffer connected to displays and other plugins
    m_pNormalizeOutput->measurementData()->clear();
    m_pNormalizeBuffer->clear();

    m_bPluginControlWidgetsInit = false;

    return true;
}

//=============================================================================================================

AbstractPlugin::PluginType Normalize::getType() const
{
    return _IAlgorithm;
}

//=============================================================================================================

QString Normalize::getName() const
{
    return "Normalize";
}

//=============================================================================================================

QWidget* Normalize::setupWidget()
{
    NormalizeSetupWidget* setupWidget = new NormalizeSetupWidget(this, QString("MNESCAN/%1/").arg(this->getName()));
    return setupWidget;
}

//=============================================================================================================

void Normalize::update(SCMEASLIB::Measurement::SPtr pMeasurement)
{
    QSharedPointer<RealTimeMultiSampleArray> pRTMSA = pMeasurement.dynamicCast<RealTimeMultiSampleArray>();

    if(pRTMSA) {
        //Check if buffer initialized
        if(!m_pNormalizeBuffer) {
            m_pNormalizeBuffer = CircularBuffer<Eigen::MatrixXd>::SPtr(new CircularBuffer<Eigen::MatrixXd>(64));
        }


        //Fiff information
        if(!m_pFiffInfo) {
            m_pFiffInfo = pRTMSA->info();

            m_pNormalizeOutput->measurementData()->initFromFiffInfo(m_pFiffInfo);
            m_pNormalizeOutput->measurementData()->setMultiArraySize(1);

            if(m_iNChannels == -1)
                m_iNChannels = m_pFiffInfo->nchan;
            if(m_dDataSampFreq == -1)
                m_dDataSampFreq = m_pFiffInfo->sfreq;

        }

        // Check if data is present
        if(pRTMSA->getMultiSampleArray().size() > 0) {

            if(!(m_iDataBlockLength > 0))
                m_iDataBlockLength = pRTMSA->getMultiSampleArray()[0].cols();

            if(!m_bPluginControlWidgetsInit) {
                initPluginControlWidgets();
            }

            for(unsigned char i = 0; i < pRTMSA->getMultiSampleArray().size(); ++i) {
                // Please note that we do not need a copy here since this function will block until
                // the buffer accepts new data again. Hence, the data is not deleted in the actual
                // Measurement function after it emitted the notify signal.
                while(!m_pNormalizeBuffer->push(pRTMSA->getMultiSampleArray()[i])) {
                    //Do nothing until the circular buffer is ready to accept new data again
                }
            }
        }
    }
}

//=============================================================================================================

void Normalize::initPluginControlWidgets()
{
    QList<QWidget*> plControlWidgets;

    // The plugin's control widget
    NormalizeWidget* pNormalizeWidget = new NormalizeWidget(QString("MNESCAN/%1/").arg(this->getName()),
                                                            m_dDataSampFreq, m_iDataBlockLength, m_bNormalizeEnabled,
                                                            m_iNormalizeMethod, m_iNormalizeModality, m_dIntervallLength);
    pNormalizeWidget->setObjectName("group_tab_Settings_Normalize");

    connect(pNormalizeWidget,&NormalizeWidget::changeEnabled,this,&Normalize::changeEnabled);
    connect(pNormalizeWidget,&NormalizeWidget::changeModality,this,&Normalize::changeModality);
    connect(pNormalizeWidget,&NormalizeWidget::changeMethod,this,&Normalize::changeMethod);
    connect(pNormalizeWidget,&NormalizeWidget::changeIntervallLength,this,&Normalize::changeIntervallLength);

    plControlWidgets.append(pNormalizeWidget);

    emit pluginControlWidgetsChanged(plControlWidgets, this->getName());

    m_bPluginControlWidgetsInit = true;
}

//=============================================================================================================

void Normalize::run()
{
    MatrixXd matData;
    MatrixXd matDataCleaned;

    // Wait for Fiff Info
    while(!m_pFiffInfo) {
        msleep(10);
    }

    while(!isInterruptionRequested()) {

        QVector<Eigen::VectorXd> pDataSums;
        QVector<Eigen::VectorXd> pDataSqSums;
        QVector<double> pDataLengthS;
        // Get the current data
        if(m_pNormalizeBuffer->pop(matData)) {

            if(m_iDataBlockLength==-1) {
                m_iDataBlockLength = matData.cols();
                matDataCleaned = Eigen::MatrixXd::Zero(matData.rows(),matData.cols());
                // calculate number of data blocks here?
            }

            if(true) { // add evaluation to leave out intervalls here
                pDataSums.push_back(matData.rowwise().sum());
                pDataSqSums.push_back(matData.cwiseAbs2().rowwise().sum());
                pDataLengthS.push_back(matData.cols()/m_dDataSampFreq);
            }

            Eigen::VectorXd dTempMean = std::accumulate(pDataSums.begin(),pDataSums.end(),Eigen::VectorXd::Zero(pDataSums.length()).eval());
            dTempMean *= 1.0/(m_iNDataBlocks*m_iDataBlockLength);
            Eigen::VectorXd dTempVar = std::accumulate(pDataSqSums.begin(),pDataSqSums.end(),Eigen::VectorXd::Zero(pDataSqSums.length()).eval()) - dTempMean.cwiseAbs2();
            dTempVar *= 1.0/(m_iNDataBlocks*m_iDataBlockLength);
            dTempVar.cwiseSqrt().cwiseInverse();

            if(m_bNormalizeEnabled) {
                if(m_iNormalizeMethod > 0)
                    matDataCleaned = matData.colwise() - dTempMean;
                if(m_iNormalizeMethod == 2)
                    matDataCleaned = dTempVar.asDiagonal()*matDataCleaned;
            } else
                matDataCleaned = matData;

            while(std::accumulate(pDataLengthS.begin(),pDataLengthS.end(),0.0) < m_dIntervallLength) {
                pDataLengthS.removeFirst();
                pDataSums.removeFirst();
                pDataSqSums.removeFirst();
            }

            //Send the data to the connected plugins and the online display
            //Unocmment this if you also uncommented the m_pOutput in the constructor above
            if(!isInterruptionRequested()) {
                m_pNormalizeOutput->measurementData()->setValue(matData);
            }
        }
    }
}

//=============================================================================================================

void Normalize::changeModality(int modality)
{
    QMutexLocker locker(&m_qMutex);

    m_iNormalizeModality = modality;
}

//=============================================================================================================

void Normalize::changeMethod(int method)
{
    QMutexLocker locker(&m_qMutex);

    m_iNormalizeMethod = method;
}

//=============================================================================================================

void Normalize::changeEnabled(bool enabled)
{
    QMutexLocker locker(&m_qMutex);

    m_bNormalizeEnabled = enabled;
}

//=============================================================================================================

void Normalize::changeIntervallLength(double length)
{
    QMutexLocker locker(&m_qMutex);

    m_dIntervallLength = length;
}
