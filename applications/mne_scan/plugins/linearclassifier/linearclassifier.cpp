//=============================================================================================================
/**
 * @file     linearclassifier.cpp
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
 * @brief    Definition of the LinearClassifier class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "linearclassifier.h"

#include "FormFiles/linearclassifiersetupwidget.h"
#include "FormFiles/linearclassifierwidget.h"

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

using namespace LINEARCLASSIFIERPLUGIN;
using namespace SCSHAREDLIB;
using namespace SCMEASLIB;
//using namespace DISPLIB;
using namespace FIFFLIB;
using namespace UTILSLIB;
using namespace Eigen;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

LinearClassifier::LinearClassifier()
: m_pLinearClassifierBuffer(CircularBuffer<Eigen::MatrixXd>::SPtr::create(64))
{
}

//=============================================================================================================

LinearClassifier::~LinearClassifier()
{
    if(this->isRunning()) {
        stop();
    }
}


//=============================================================================================================

QSharedPointer<AbstractPlugin> LinearClassifier::clone() const
{
    QSharedPointer<LinearClassifier> pLinearClassifierClone(new LinearClassifier);
    return pLinearClassifierClone;
}

//=============================================================================================================

void LinearClassifier::init()
{
    QMutexLocker locker(&m_qMutex);

    // Input
    m_pLinearClassifierInput = PluginInputData<RealTimeMultiSampleArray>::create(this, "LinearClassifierIn", "LinearClassifier input data");
    connect(m_pLinearClassifierInput.data(), &PluginInputConnector::notify,
            this, &LinearClassifier::update, Qt::DirectConnection);
    m_inputConnectors.append(m_pLinearClassifierInput);

    // Output - Uncomment this if you don't want to send processed data (in form of a matrix) to other plugins.
    // Also, this output stream will generate an online display in your plugin
    m_pLinearClassifierOutput = PluginOutputData<RealTimeMultiSampleArray>::create(this, "LinearClassifierOut", "LinearClassifier output data");
    m_pLinearClassifierOutput->measurementData()->setName(this->getName());
    m_outputConnectors.append(m_pLinearClassifierOutput);

    //Delete Buffer - will be initailzed with first incoming data
    if(!m_pLinearClassifierBuffer.isNull())
        m_pLinearClassifierBuffer = CircularBuffer<Eigen::MatrixXd>::SPtr();
}

//=============================================================================================================

void LinearClassifier::unload()
{
}

//=============================================================================================================

bool LinearClassifier::start()
{
    //Start thread
    QThread::start();

    return true;
}

//=============================================================================================================

bool LinearClassifier::stop()
{
    requestInterruption();
    wait(500);

    // Clear all data in the buffer connected to displays and other plugins
    m_pLinearClassifierOutput->measurementData()->clear();
    m_pLinearClassifierBuffer->clear();

    m_bPluginControlWidgetsInit = false;

    return true;
}

//=============================================================================================================

AbstractPlugin::PluginType LinearClassifier::getType() const
{
    return _IAlgorithm;
}

//=============================================================================================================

QString LinearClassifier::getName() const
{
    return "LinearClassifier";
}

//=============================================================================================================

QWidget* LinearClassifier::setupWidget()
{
    LinearClassifierSetupWidget* setupWidget = new LinearClassifierSetupWidget(this, QString("MNESCAN/%1/").arg(this->getName()), false);
    return setupWidget;
}

//=============================================================================================================

void LinearClassifier::update(SCMEASLIB::Measurement::SPtr pMeasurement)
{
    QSharedPointer<RealTimeMultiSampleArray> pRTMSA = pMeasurement.dynamicCast<RealTimeMultiSampleArray>();

    if(pRTMSA) {
        //Check if buffer initialized
        if(!m_pLinearClassifierBuffer) {
            m_pLinearClassifierBuffer = CircularBuffer<Eigen::MatrixXd>::SPtr(new CircularBuffer<Eigen::MatrixXd>(64));
        }


        //Fiff information
        if(!m_pFiffInfo) {
            m_pFiffInfo = pRTMSA->info();

            m_pFiffInfo->filename = "";
            m_pFiffInfo->bads.clear();
            m_pFiffInfo->nchan = m_iNOutputChannels;

            QList<FiffChInfo> fakeChList;
            FiffChInfo fakeCh;
            fakeCh.ch_name = "LC ";
            fakeCh.kind = 502;
            fakeCh.range = -1;
            fakeCh.unit = -1;

            QStringList fakeChNames;

            for(int i=0; i<m_iNOutputChannels; ++i)
                {
                    fakeCh.ch_name = QString("LC %1").arg(i);
                    fakeChList.append(fakeCh);
                    fakeChNames.append(fakeCh.ch_name);
                }

            m_pFiffInfo->chs = fakeChList;
            m_pFiffInfo->ch_names = fakeChNames;
            m_pFiffInfo->file_id = FIFFLIB::FiffId::new_file_id(); //check if necessary

            m_pLinearClassifierOutput->measurementData()->initFromFiffInfo(m_pFiffInfo);
            m_pLinearClassifierOutput->measurementData()->setMultiArraySize(1);
        }

        if(!m_bPluginControlWidgetsInit) {
            initPluginControlWidgets();
        }

        // Check if data is present
        if(pRTMSA->getMultiSampleArray().size() > 0) {

            for(unsigned char i = 0; i < pRTMSA->getMultiSampleArray().size(); ++i) {
                // Please note that we do not need a copy here since this function will block until
                // the buffer accepts new data again. Hence, the data is not deleted in the actual
                // Measurement function after it emitted the notify signal.
                while(!m_pLinearClassifierBuffer->push(pRTMSA->getMultiSampleArray()[i])) {
                    //Do nothing until the circular buffer is ready to accept new data again
                }
            }
        }
    }
}

//=============================================================================================================

void LinearClassifier::initPluginControlWidgets()
{
    QList<QWidget*> plControlWidgets;

    // The plugin's control widget
    LinearClassifierWidget* pLinearClassifierWidget = new LinearClassifierWidget(this, QString("MNESCAN/%1/").arg(this->getName()), false);
    pLinearClassifierWidget->setObjectName("group_tab_Settings_LinearClassifier");

    plControlWidgets.append(pLinearClassifierWidget);

    emit pluginControlWidgetsChanged(plControlWidgets, this->getName());

    m_bPluginControlWidgetsInit = true;
}

//=============================================================================================================

void LinearClassifier::run()
{
    MatrixXd matData;

    // Wait for Fiff Info
    while(!m_pFiffInfo) {
        msleep(10);
    }

    while(!isInterruptionRequested()) {
        // Get the current data
        if(m_pLinearClassifierBuffer->pop(matData)) {

            MatrixXd outMatData;

            if(m_pLinearClassifierMatrix.cols()==matData.rows())
                outMatData = m_pLinearClassifierMatrix*matData;
            else if(m_pLinearClassifierMatrix.cols()<matData.rows()) {
                qWarning() << "Number of input channels larger than number of columns of classifier matrix; only using the first " << m_iNInputChannels << "inputs.";
                outMatData = m_pLinearClassifierMatrix*matData.topRows(m_iNInputChannels);
            } else {
                qWarning() << "Number of input channels smaller than number of columns of classifier matrix; only using the first " << m_iNInputChannels << "inputs.";
                outMatData = m_pLinearClassifierMatrix.leftCols(matData.rows())*matData;
            }

            //Send the data to the connected plugins and the online display
            //Unocmment this if you also uncommented the m_pOutput in the constructor above
            if(!isInterruptionRequested()) {
                m_pLinearClassifierOutput->measurementData()->setValue(matData);
            }
        }
    }
}

//=============================================================================================================

void LinearClassifier::setNInputChannels(int iNInputChannels)
{
    if(!this->isRunning())
        m_iNInputChannels = iNInputChannels;
    else
        qWarning() << "Cannot change number of input channels while plugin is running.";
}

//=============================================================================================================

void LinearClassifier::setNOutputChannels(int iNOutputChannels)
{
    if(!this->isRunning())
        m_iNOutputChannels = iNOutputChannels;
    else
        qWarning() << "Cannot change number of output channels while plugin is running.";
}

//=============================================================================================================

void LinearClassifier::setLinearClassifierMatrix(Eigen::MatrixXd matrix)
{
    if(matrix.cols() == m_iNInputChannels && matrix.rows()==m_iNOutputChannels)
        m_pLinearClassifierMatrix = matrix;
    else
        qWarning() << "Linear classifier matrix does not fit number of in- and output channels";
}
