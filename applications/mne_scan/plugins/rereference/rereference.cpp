//=============================================================================================================
/**
 * @file     rereference.cpp
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
 * @brief    Definition of the Rereference class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "rereference.h"

#include "FormFiles/rereferencesetupwidget.h"
#include "FormFiles/rereferencewidget.h"

#include <disp/viewers/channelselectionview.h>
#include <disp/viewers/helpers/channelinfomodel.h>

#include <utils/ioutils.h>

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

using namespace REREFERENCEPLUGIN;
using namespace SCSHAREDLIB;
using namespace SCMEASLIB;
using namespace DISPLIB;
using namespace FIFFLIB;
using namespace IOBUFFER;
using namespace UTILSLIB;
using namespace Eigen;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

Rereference::Rereference()
: m_pRereferenceBuffer(CircularBuffer<Eigen::MatrixXd>::SPtr::create(64))
, m_iNChannels(-1)
, m_iRereferenceModality(0)
, m_iRereferenceMethod(0)
, m_bRereferenceEnabled(true)
, m_sMatrixFilename("")
{
}

//=============================================================================================================

Rereference::~Rereference()
{
    if(this->isRunning()) {
        stop();
    }
}


//=============================================================================================================

QSharedPointer<IPlugin> Rereference::clone() const
{
    QSharedPointer<Rereference> pRereferenceClone(new Rereference);
    return pRereferenceClone;
}

//=============================================================================================================

void Rereference::init()
{
    QMutexLocker locker(&m_qMutex);

    // Input
    m_pRereferenceInput = PluginInputData<RealTimeMultiSampleArray>::create(this, "RereferenceIn", "Rereference input data");
    connect(m_pRereferenceInput.data(), &PluginInputConnector::notify,
            this, &Rereference::update, Qt::DirectConnection);
    m_inputConnectors.append(m_pRereferenceInput);

    // Output - Uncomment this if you don't want to send processed data (in form of a matrix) to other plugins.
    // Also, this output stream will generate an online display in your plugin
    m_pRereferenceOutput = PluginOutputData<RealTimeMultiSampleArray>::create(this, "ReferenceOut", "Reference output data");
    m_pRereferenceOutput->data()->setName(this->getName());
    m_outputConnectors.append(m_pRereferenceOutput);

    //Delete Buffer - will be initailzed with first incoming data
    if(!m_pRereferenceBuffer.isNull())
        m_pRereferenceBuffer = CircularBuffer<Eigen::MatrixXd>::SPtr();

    //add button
    m_pActionSelectSensors = new QAction(QIcon(":/images/selectSensors.png"), tr("Show the channel selection view for rereferencing"),this);
    m_pActionSelectSensors->setToolTip(tr("Show the channel selection view for rereferencing"));
    connect(m_pActionSelectSensors.data(), &QAction::triggered,
            this, &Rereference::showSensorSelectionWidget);
    addPluginAction(m_pActionSelectSensors);
    m_pActionSelectSensors->setVisible(true);
}

//=============================================================================================================

void Rereference::unload()
{
}

//=============================================================================================================

bool Rereference::start()
{
    //Start thread
    QThread::start();

    return true;
}

//=============================================================================================================

bool Rereference::stop()
{
    requestInterruption();
    wait(500);

    // Clear all data in the buffer connected to displays and other plugins
    m_pRereferenceOutput->data()->clear();
    m_pRereferenceBuffer->clear();

    m_bPluginControlWidgetsInit = false;

    return true;
}

//=============================================================================================================

IPlugin::PluginType Rereference::getType() const
{
    return _IAlgorithm;
}

//=============================================================================================================

QString Rereference::getName() const
{
    return "Rereference";
}

//=============================================================================================================

QWidget* Rereference::setupWidget()
{
    RereferenceSetupWidget* setupWidget = new RereferenceSetupWidget(this, QString("MNESCAN/%1/").arg(this->getName()));
    return setupWidget;
}

//=============================================================================================================

void Rereference::update(SCMEASLIB::Measurement::SPtr pMeasurement)
{
    QSharedPointer<RealTimeMultiSampleArray> pRTMSA = pMeasurement.dynamicCast<RealTimeMultiSampleArray>();

    if(pRTMSA) {
        //Check if buffer initialized
        if(!m_pRereferenceBuffer) {
            m_pRereferenceBuffer = CircularBuffer<Eigen::MatrixXd>::SPtr(new CircularBuffer<Eigen::MatrixXd>(64));
        }


        //Fiff information
        if(!m_pFiffInfo) {
            m_pFiffInfo = pRTMSA->info();

            m_pRereferenceOutput->data()->initFromFiffInfo(m_pFiffInfo);
            m_pRereferenceOutput->data()->setMultiArraySize(1);

            m_pMEGindex = m_pFiffInfo->pick_types(true,false,false);
            m_pEEGindex = m_pFiffInfo->pick_types(false,true,false);
            m_pMEEGindex = m_pFiffInfo->pick_types(true,true,false);

            qDebug() << "EEG Channels for rereferencing";
            for(int i=0; i<m_pEEGindex.size();++i)
                qDebug() << i << " " << m_pFiffInfo->ch_names.at(m_pEEGindex[i]);

            if(m_iNChannels == -1)
                m_iNChannels = m_pFiffInfo->nchan;

            updateRereferenceMatrix();
        }

        if(!m_bPluginControlWidgetsInit) {
            initPluginControlWidgets();
        }

        // Check if data is present
        if(pRTMSA->getMultiSampleArray().size() > 0) {
            //Init widgets
            if(m_iNChannels != pRTMSA->getMultiSampleArray().first().rows()) {
                qWarning() << "[Rereference::update] Data matrix has wrong format!";
            }

            for(unsigned char i = 0; i < pRTMSA->getMultiSampleArray().size(); ++i) {
                // Please note that we do not need a copy here since this function will block until
                // the buffer accepts new data again. Hence, the data is not deleted in the actual
                // Measurement function after it emitted the notify signal.
                while(!m_pRereferenceBuffer->push(pRTMSA->getMultiSampleArray()[i])) {
                    //Do nothing until the circular buffer is ready to accept new data again
                }
            }
        }
    }
}

//=============================================================================================================

void Rereference::initPluginControlWidgets()
{
    QList<QWidget*> plControlWidgets;

    // The plugin's control widget
    RereferenceWidget* pRereferenceWidget = new RereferenceWidget(QString("MNESCAN/%1/").arg(this->getName()));
    pRereferenceWidget->setObjectName("group_tab_Settings_Rereference");

    connect(pRereferenceWidget,&RereferenceWidget::changeEnabled,this,&Rereference::changeEnabled);
    connect(pRereferenceWidget,&RereferenceWidget::changeModality,this,&Rereference::changeModality);
    connect(pRereferenceWidget,&RereferenceWidget::changeMethod,this,&Rereference::changeMethod);

    plControlWidgets.append(pRereferenceWidget);

    //Channelselectionwidget

    //Init channel selection manager
    m_pChannelInfoModel = ChannelInfoModel::SPtr::create(m_pFiffInfo,
                                                         this);

    m_pChannelSelectionView = ChannelSelectionView::SPtr::create(QString("MNESCAN/%1/").arg(this->getName()),
                                                                 nullptr,
                                                                 m_pChannelInfoModel,
                                                                 Qt::Window);
    m_pChannelSelectionView->setWindowTitle(tr(QString("%1: Channel Selection Window").arg(this->getName()).toUtf8()));

    connect(m_pChannelSelectionView.data(), &ChannelSelectionView::loadedLayoutMap,
            m_pChannelInfoModel.data(), &ChannelInfoModel::layoutChanged);

    connect(m_pChannelInfoModel.data(), &ChannelInfoModel::channelsMappedToLayout,
            m_pChannelSelectionView.data(), &ChannelSelectionView::setCurrentlyMappedFiffChannels);

    connect(m_pChannelSelectionView.data(), &ChannelSelectionView::showSelectedChannelsOnly,
            this, &Rereference::evaluateSelectedChannelsOnly);

    m_pChannelInfoModel->layoutChanged(m_pChannelSelectionView->getLayoutMap());

    emit pluginControlWidgetsChanged(plControlWidgets, this->getName());

    m_bPluginControlWidgetsInit = true;
}

//=============================================================================================================

void Rereference::run()
{
    MatrixXd matData;

    // Wait for Fiff Info
    while(!m_pFiffInfo) {
        msleep(10);
    }

    while(!isInterruptionRequested()) {
        // Get the current data
        if(m_pRereferenceBuffer->pop(matData)) {

            if(m_bRereferenceEnabled)
                matData = m_pRereferenceMatrix * matData;

            //Send the data to the connected plugins and the online display
            //Unocmment this if you also uncommented the m_pOutput in the constructor above
            if(!isInterruptionRequested()) {
                m_pRereferenceOutput->data()->setValue(matData);
            }
        }
    }
}

//=============================================================================================================

void Rereference::evaluateSelectedChannelsOnly(const QStringList &selectedChannels)
{
    QMutexLocker locker(&m_qMutex);

    m_pSelectedChannels.clear();

    //Add selected channels to list
    for(int i = 0; i<m_pFiffInfo->ch_names.length(); i++) {
        QString channel = m_pFiffInfo->ch_names.value(i);

        if(selectedChannels.contains(channel)) {
            m_pSelectedChannels.append(i);
        }
    }

    if(m_iNChannels > 0)
        updateRereferenceMatrix();
}

void Rereference::updateRereferenceMatrix()
{
    m_pRereferenceMatrix = MatrixXd::Identity(m_iNChannels,m_iNChannels);

    switch(m_iRereferenceMethod) {
    case 0:
    {
        Eigen::RowVectorXd rereferenceRow = RowVectorXd::Zero(m_iNChannels);
        double rereferenceValue = 0.0;
        if(m_iRereferenceModality == 0 || m_iRereferenceModality == 2) {
            rereferenceValue = -1.0/m_pEEGindex.size();
            for(int i=0; i<m_pEEGindex.size(); ++i)
                rereferenceRow(m_pEEGindex[i]) = rereferenceValue;
        }
        if(m_iRereferenceModality == 1 || m_iRereferenceModality == 2) {
            rereferenceValue = -1.0/m_pMEGindex.size();
            for(int i=0; i<m_pMEGindex.size(); ++i)
                rereferenceRow(m_pMEGindex[i]) = rereferenceValue;
        }
        if(m_iRereferenceModality == 0 || m_iRereferenceModality == 2)
            for(int i=0; i<m_pEEGindex.size(); ++i)
                m_pRereferenceMatrix.row(i) += rereferenceRow;
        if(m_iRereferenceModality == 1 || m_iRereferenceModality == 2)
            for(int i=0; i<m_pMEGindex.size(); ++i)
                m_pRereferenceMatrix.row(i) += rereferenceRow;
        break;
    }
    case 1:
    {
        //check whether conversion to sparse matrix makes sense
        Eigen::RowVectorXd rereferenceRowEEG = RowVectorXd::Zero(m_iNChannels);
        Eigen::RowVectorXd rereferenceRowMEG = RowVectorXd::Zero(m_iNChannels);

        int iEEGreferences = 0;
        int iMEGreferences = 0;
        for(int i=0; i<m_pSelectedChannels.length(); ++i){
            if((m_pEEGindex.array() == m_pSelectedChannels.at(i)).any())
                iEEGreferences++;
            else if((m_pMEGindex.array() == m_pSelectedChannels.at(i)).any())
                iMEGreferences++;
        }

        double EEGfactor = 1.0/iEEGreferences;
        double MEGfactor = 1.0/iMEGreferences;

        for(int i=0; i<m_pSelectedChannels.length(); ++i){
            if((m_pEEGindex.array() == m_pSelectedChannels.at(i)).any())
                rereferenceRowEEG[m_pSelectedChannels.at(i)] = -1.0*EEGfactor;
            else if((m_pMEGindex.array() == m_pSelectedChannels.at(i)).any())
                rereferenceRowMEG[m_pSelectedChannels.at(i)] = -1.0*MEGfactor;
        }

        if(m_iRereferenceModality == 0 || m_iRereferenceModality == 2)
            for(int i=0; i<m_pEEGindex.size(); ++i)
                m_pRereferenceMatrix.row(i) += rereferenceRowEEG;
        if(m_iRereferenceModality == 1 || m_iRereferenceModality == 2)
            for(int i=0; i<m_pMEGindex.size(); ++i)
                m_pRereferenceMatrix.row(i) += rereferenceRowMEG;
        break;
    }
    case 2:
    {
        //check whether conversion to sparse matrix makes sense
        Eigen::MatrixXd tempRereference;
        if(!IOUtils::read_eigen_matrix(tempRereference,m_sMatrixFilename))
        {
            qWarning() << "[Rereference::updateRereferenceMatrix] Could not read rereference matrix file " << m_sMatrixFilename;
            break;
        }
        if(tempRereference.rows()!=tempRereference.cols())
        {
            qWarning() << "[Rereference::updateRereferenceMatrix] Rereference matrix has to be square";
                          break;
        }

        switch(m_iRereferenceModality) {
        case 0:
            if(m_pEEGindex.size() == tempRereference.cols())
            {
                for(int i=0; i<m_pEEGindex.size(); ++i)
                    for(int j=0; j<m_pEEGindex.size(); ++j)
                    m_pRereferenceMatrix(m_pEEGindex(i),m_pEEGindex(j)) = tempRereference(i,j);
            } else if(m_pMEEGindex.size() == tempRereference.rows())
            {
                for(int i=0; i<m_pEEGindex.size(); ++i)
                    m_pRereferenceMatrix.row(m_pEEGindex[i]) = tempRereference.row(i);

            } else {
                qWarning() << "[Rereference::updateRereferenceMatrix] EEG rereferencing selected, but rereferencing matrix does not match number of EEG channels, which is " << m_pEEGindex.size();
            }
            break;
        case 1:
            if(m_pMEGindex.size() == tempRereference.cols())
            {
                for(int i=0; i<m_pMEGindex.size(); ++i)
                    for(int j=0; j<m_pMEGindex.size(); ++j)
                    m_pRereferenceMatrix(m_pMEGindex(i),m_pMEGindex(j)) = tempRereference(i,j);
            } else if(m_pMEEGindex.size() == tempRereference.rows())
            {
                for(int i=0; i<m_pMEGindex.size(); ++i)
                    m_pRereferenceMatrix.row(m_pMEGindex[i]) = tempRereference.row(i);

            } else {
                qWarning() << "[Rereference::updateRereferenceMatrix] MEG rereferencing selected, but rereferencing matrix does not match number of MEG channels, which is " << m_pMEGindex.size();
            }
            break;
        case 2:
            if(m_pMEEGindex.size() != tempRereference.rows())
            {
                qWarning() << "[Rereference::updateRereferenceMatrix] EEG + MEG rereferencing selected, but rereferencing matrix does not match number of EEG and MEG channels, which is " << m_pMEEGindex.size();
            } else {
                for(int i=0; i<m_pMEEGindex.size(); ++i)
                    m_pRereferenceMatrix.row(m_pMEEGindex[i]) = tempRereference.row(i);
            }
            break;
        }
        break;
    }
    }

    qDebug() << "m_pRereferenceMatrix min:" << m_pRereferenceMatrix.minCoeff() << " max:" << m_pRereferenceMatrix.maxCoeff();
}

//=============================================================================================================

void Rereference::showSensorSelectionWidget()
{
    if(m_pChannelSelectionView->isActiveWindow()) {
        m_pChannelSelectionView->hide();
    } else {
        m_pChannelSelectionView->activateWindow();
        m_pChannelSelectionView->show();
    }
}

//=============================================================================================================

void Rereference::changeModality(int modality)
{
    QMutexLocker locker(&m_qMutex);

    m_iRereferenceModality = modality;

    if(m_iNChannels > 0)
        updateRereferenceMatrix();
}

//=============================================================================================================

void Rereference::changeMethod(int method)
{
    QMutexLocker locker(&m_qMutex);

    m_iRereferenceMethod = method;

    if(m_iNChannels > 0)
        updateRereferenceMatrix();
}

//=============================================================================================================

void Rereference::changeEnabled(bool enabled)
{
    QMutexLocker locker(&m_qMutex);

    m_bRereferenceEnabled = enabled;

    if(m_iNChannels > 0)
        updateRereferenceMatrix();
}

//=============================================================================================================

void Rereference::changeMatrixFilename(QString sMatrixFilename)
{
    QMutexLocker locker(&m_qMutex);

    m_sMatrixFilename = sMatrixFilename;

    if(m_iNChannels > 0)
        updateRereferenceMatrix();
}
