//=============================================================================================================
/**
 * @file     bandpower.cpp
 * @author   Johannes Vorwerk <johannes.vorwerk@umit.at>;
 *           Gabriel B Motta <gabrielbenmotta@gmail.com>;
 *           Lorenz Esch <lesch@mgh.harvard.edu>
 * @version  dev
 * @date     April, 2020
 *
 * @section  LICENSE
 *
 * Copyright (C) 2020, Johannes Vorwerk, Gabriel B Motta, Lorenz Esch. All rights reserved.
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
 * @brief    Definition of the BandPower class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "bandpower.h"

#include "FormFiles/bandpowersetupwidget.h"

#include "disp/viewers/bandpowersettingsview.h"
#include "disp/viewers/arsettingsview.h"

#include <fiff/fiff_constants.h>
#include <fiff/fiff.h>

#include <disp/viewers/channelselectionview.h>
#include <disp/viewers/helpers/channelinfomodel.h>

#include <chrono>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace BANDPOWERPLUGIN;
using namespace SCMEASLIB;
using namespace UTILSLIB;
using namespace IOBUFFER;
using namespace DISPLIB;
//using namespace RTPROCESSINGLIB;
using namespace FIFFLIB;
using namespace SCSHAREDLIB;
using namespace Eigen;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

BandPower::BandPower()
    : m_bProcessData(false)
    , m_pBandPowerInput(NULL)
    , m_pBandPowerOutput(NULL)
    , m_pBandPowerBuffer(CircularBuffer<Eigen::MatrixXd>::SPtr())
    , m_iIntervallLengthFactor(5)
    , m_iResolutionFactorFFT(1)
    , m_iOrderAR(40)
    , m_iEvalsAR(100)
    , m_sSpectrumMethod("AR")
    , m_iDetrendMethod(0)
    , m_dDataSampFreq(-1)
    , m_dFreqMin(10)
    , m_dFreqMax(30)
    , m_iNTimeSteps(-1)
    , m_iNChannels(-1)
    , m_iBandPowerChannels(1)
    , m_iBandPowerBins(1)
{
}

//=============================================================================================================

BandPower::~BandPower()
{
    if(this->isRunning())
        stop();
}

//=============================================================================================================

QSharedPointer<IPlugin> BandPower::clone() const
{
    QSharedPointer<BandPower> pBandPowerClone(new BandPower);
    return pBandPowerClone;
}

//=============================================================================================================

void BandPower::init()
{
    QMutexLocker locker(&m_qMutex);

    // Input
    m_pBandPowerInput = PluginInputData<RealTimeMultiSampleArray>::create(this, "BandPowerIn", "BandPower input data");
    connect(m_pBandPowerInput.data(), &PluginInputConnector::notify, this, &BandPower::update, Qt::DirectConnection);
    m_inputConnectors.append(m_pBandPowerInput);

    // Output - Uncomment this if you don't want to send processed data (in form of a matrix) to other plugins.
    // Also, this output stream will generate an online display in your plugin
    m_pBandPowerOutput = PluginOutputData<RealTimeMultiSampleArray>::create(this, "BandPowerOut", "BandPower output data");
    m_pBandPowerOutput->data()->setName(this->getName());//Provide name to auto store widget settings
    m_pBandPowerOutput->data()->setMultiArraySize(1);

    m_outputConnectors.append(m_pBandPowerOutput);

    //Delete Buffer - will be initailzed with first incoming data
    if(!m_pBandPowerBuffer.isNull())
        m_pBandPowerBuffer = CircularBuffer<Eigen::MatrixXd>::SPtr();

    //add button
    m_pActionSelectSensors = new QAction(QIcon(":/images/selectSensors.png"), tr("Show the channel selection view for bandpower computation"),this);
    m_pActionSelectSensors->setToolTip(tr("Show the channel selection view for bandpower computation"));
    connect(m_pActionSelectSensors.data(), &QAction::triggered,
            this, &BandPower::showSensorSelectionWidget);
    addPluginAction(m_pActionSelectSensors);
    m_pActionSelectSensors->setVisible(true);
}

//=============================================================================================================

void BandPower::unload()
{
}

//=============================================================================================================

bool BandPower::start()
{
    //Start thread
    //QThread::start(); //Don't start here, so that we don't have to wait for pFiffInfo later

    return true;
}


//=============================================================================================================

bool BandPower::stop()
{
    requestInterruption();
    wait(500);

    m_pBandPowerOutput->data()->clear();
    m_pBandPowerBuffer->clear();

    return true;
}

//=============================================================================================================

IPlugin::PluginType BandPower::getType() const
{
    return _IAlgorithm;
}

//=============================================================================================================

QString BandPower::getName() const
{
    return "BandPower";
}

//=============================================================================================================

QWidget* BandPower::setupWidget()
{
    //QWidget* setupWidget = new SpectrumSetupWidget(this);//widget is later distroyed by CentralWidget - so it has to be created everytime new
    BandPowerSetupWidget* setupWidget = new BandPowerSetupWidget(this);//widget is later distroyed by CentralWidget - so it has to be created everytime new

    return setupWidget;
}

//=============================================================================================================

void BandPower::changeMinMaxFrequency(double minFrequency, double maxFrequency)
{
    QMutexLocker locker(&m_qMutex);

    m_dFreqMin = minFrequency;
    m_dFreqMax = maxFrequency;
}

//=============================================================================================================

void BandPower::changeSpectrumMethod(QString method)
{
    QMutexLocker locker(&m_qMutex);

    m_sSpectrumMethod = method;
}
//=============================================================================================================

void BandPower::changeIntervallLengthFactor(quint32 lengthfactor)
{
    QMutexLocker locker(&m_qMutex);

    m_iIntervallLengthFactor = lengthfactor;
}

//=============================================================================================================

void BandPower::changeResolutionFactor(quint32 resolutionfactor)
{
    QMutexLocker locker(&m_qMutex);

    m_iResolutionFactorFFT = resolutionfactor;
}

//=============================================================================================================

void BandPower::changeDetrendMethod(quint32 method)
{
    QMutexLocker locker(&m_qMutex);

    m_iDetrendMethod = method;
}

//=============================================================================================================

void BandPower::changeBandPowerChannels(quint32 channels)
{
    QMutexLocker locker(&m_qMutex);

    m_iBandPowerChannels = channels;

    qDebug() << "m_iBandPowerChannels " << m_iBandPowerChannels;
}

//=============================================================================================================

void BandPower::changeBandPowerBins(quint32 bins)
{
    QMutexLocker locker(&m_qMutex);

    m_iBandPowerBins = bins;

    qDebug() << "m_iBandPowerBins " << m_iBandPowerBins;
}

//=============================================================================================================

void BandPower::changeAROrder(quint32 order)
{
    QMutexLocker locker(&m_qMutex);

    if(order < 1){
        qDebug() << "AR order < 1 not allowed.";
        m_iOrderAR = 1;
    } else
        m_iOrderAR = order;
}

//=============================================================================================================

void BandPower::changeARNumEvaluationPoints(quint32 evaluationpoints)
{
    QMutexLocker locker(&m_qMutex);

    m_iEvalsAR = evaluationpoints;
}

//=============================================================================================================

void BandPower::evaluateSelectedChannelsOnly(const QStringList &selectedChannels)
{
    QMutexLocker locker(&m_qMutex);

    m_pSelectedChannels.clear();

    qDebug() << "selected channels";

    //Add selected channels to list
    for(int i = 0; i<m_pFiffInfo_orig->ch_names.length(); i++) {
        QString channel = m_pFiffInfo_orig->ch_names.value(i);

        if(selectedChannels.contains(channel)) {
            m_pSelectedChannels.append(i);
            qDebug() << i;
        }
    }
}

//=============================================================================================================

void BandPower::update(SCMEASLIB::Measurement::SPtr pMeasurement)
{
    QSharedPointer<RealTimeMultiSampleArray> pRTMSA = pMeasurement.dynamicCast<RealTimeMultiSampleArray>();

    if(pRTMSA) {
        //Check if buffer initialized
        if(!m_pBandPowerBuffer) {
            m_pBandPowerBuffer = CircularBuffer<Eigen::MatrixXd>::SPtr(new CircularBuffer<Eigen::MatrixXd>(64));
        }

        //Fiff information
        if(!m_pFiffInfo) {                

            m_pFiffInfo = FIFFLIB::FiffInfo::SPtr(new FIFFLIB::FiffInfo(*pRTMSA->info().data())); // pointer not working here...
            m_pFiffInfo_orig = pRTMSA->info();
            //m_pFiffInfo = FIFFLIB::FiffInfo::SPtr(new FIFFLIB::FiffInfo);

            //Init output - Uncomment this if you also uncommented the m_pDummyOutput in the constructor above
            //m_pBandPowerOutput->data()->initFromFiffInfo(m_pFiffInfo);

            m_dDataSampFreq = m_pFiffInfo->sfreq;

            m_pFiffInfo->filename = "";
            m_pFiffInfo->bads.clear();
            m_pFiffInfo->nchan = m_iBandPowerChannels * m_iBandPowerBins;

            QList<FiffChInfo> fakeChList;
            FiffChInfo fakeCh;
            /*fakeCh.ch_name = "AR";
            fakeCh.kind = 502;
            fakeCh.range = -1;
            //fakeCh.setMinValue(0.5e-20);
            //fakeCh.setMaxValue(2e-20);
            fakeCh.unit = -1;
            fakeChList.append(fakeCh);
            fakeCh.ch_name = "FFT";
            fakeChList.append(fakeCh);*/
            fakeCh.ch_name = "BP ";
            fakeCh.kind = 502;
            fakeCh.range = -1;
            //fakeCh.setMinValue(0.5e-20);
            //fakeCh.setMaxValue(2e-20);
            fakeCh.unit = -1;

            QStringList fakeChNames;

            for(int i=0; i<m_iBandPowerChannels; ++i)
                for(int j=0; j<m_iBandPowerBins; ++j)
                {
                    fakeCh.ch_name = QString("BP %1-%2").arg(i).arg(j);
                    fakeChList.append(fakeCh);
                    fakeChNames.append(fakeCh.ch_name);
                }

            m_pFiffInfo->chs = fakeChList;

            m_pFiffInfo->ch_names = fakeChNames;

            m_pFiffInfo->file_id = FIFFLIB::FiffId::new_file_id(); //check if necessary

            m_pFiffInfo->sfreq = pRTMSA->info()->sfreq/m_iNTimeSteps;

            m_pBandPowerOutput->data()->initFromFiffInfo(m_pFiffInfo);

            m_pBandPowerOutput->data()->setVisibility(true);
        }

        // Check if data is present
        if(pRTMSA->getMultiSampleArray().size() > 0) {
            //Init widgets
            if(m_iNChannels == -1) {
                m_iNChannels = pRTMSA->getMultiSampleArray().first().rows();
                initPluginControlWidgets();
                QThread::start();
            }
            if(m_iNTimeSteps == -1)
                m_iNTimeSteps = pRTMSA->getMultiSampleArray().first().cols();

            for(unsigned char i = 0; i < pRTMSA->getMultiSampleArray().size(); ++i) {
                // Please note that we do not need a copy here since this function will block until
                // the buffer accepts new data again. Hence, the data is not deleted in the actual
                // Measurement function after it emitted the notify signal.
                while(!m_pBandPowerBuffer->push(pRTMSA->getMultiSampleArray()[i])) {
                    //Do nothing until the circular buffer is ready to accept new data again
                }
            }
        }
    }
}

void BandPower::initPluginControlWidgets()
{
    QList<QWidget*> plControlWidgets;

    BandPowerSettingsView* pBandPowerSettingsView = new BandPowerSettingsView(QString("MNESCAN/%1/").arg(this->getName()),m_dDataSampFreq,
                                                                              m_dFreqMin, m_dFreqMax, m_sSpectrumMethod, m_iIntervallLengthFactor,
                                                                              m_iBandPowerChannels, m_iBandPowerBins, m_iDetrendMethod, true);
    pBandPowerSettingsView->setObjectName("group_tab_Settings_General");
    plControlWidgets.append(pBandPowerSettingsView);

    connect(pBandPowerSettingsView,&BandPowerSettingsView::changeMinMax,this,&BandPower::changeMinMaxFrequency);
    connect(pBandPowerSettingsView,&BandPowerSettingsView::changeMethod,this,&BandPower::changeSpectrumMethod);
    connect(pBandPowerSettingsView,&BandPowerSettingsView::changeIntervallLength,this,&BandPower::changeIntervallLengthFactor);
    connect(pBandPowerSettingsView,&BandPowerSettingsView::changeDetrend,this,&BandPower::changeDetrendMethod);
    connect(pBandPowerSettingsView,&BandPowerSettingsView::changeChannels,this,&BandPower::changeBandPowerChannels);
    connect(pBandPowerSettingsView,&BandPowerSettingsView::changeBins,this,&BandPower::changeBandPowerBins);

    ARSettingsView* pARSettingsView = new ARSettingsView(QString("MNESCAN/%1/").arg(this->getName()));

    pARSettingsView->setObjectName("group_tab_Settings_AR");
    plControlWidgets.append(pARSettingsView);

    connect(pARSettingsView, &ARSettingsView::changeAROrder,this,&BandPower::changeAROrder);
    connect(pARSettingsView, &ARSettingsView::changeARNumEvaluationPoints,this,&BandPower::changeARNumEvaluationPoints);

    //Channelselectionwidget

    //Init channel selection manager
    m_pChannelInfoModel = ChannelInfoModel::SPtr::create(m_pFiffInfo_orig,
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

    /*        connect(m_pChannelSelectionView.data(), &ChannelSelectionView::showSelectedChannelsOnly,
                m_pChannelDataView.data(), &RtFiffRawView::showSelectedChannelsOnly);

        connect(m_pChannelDataView.data(), &RtFiffRawView::channelMarkingChanged,
                m_pChannelSelectionView.data(), &ChannelSelectionView::updateBadChannels);

        connect(m_pChannelDataView.data(), &RtFiffRawView::selectedChannelsChanged,
                m_pChannelSelectionView.data(), &ChannelSelectionView::setUserSelection);

        connect(m_pChannelDataView.data(), &RtFiffRawView::selectedChannelsResetted,
                m_pChannelSelectionView.data(), &ChannelSelectionView::resetUserSelection);*/

    connect(m_pChannelSelectionView.data(), &ChannelSelectionView::showSelectedChannelsOnly,
            this, &BandPower::evaluateSelectedChannelsOnly);

    m_pChannelInfoModel->layoutChanged(m_pChannelSelectionView->getLayoutMap());

    emit pluginControlWidgetsChanged(plControlWidgets, this->getName());
}

//=============================================================================================================

void BandPower::run()
{
    // Length of block
    int iNSamples = m_iNTimeSteps * m_iIntervallLengthFactor;

    // Sampling frequency
    double dSampFreq = m_dDataSampFreq;

    qDebug() << "dSampFreq" << dSampFreq;

    // Determine spectrum resolution - calculation only for FFT, but we keep it for both methods here to keep the resolution
    Eigen::VectorXd FFTFreqs = Spectral::calculateFFTFreqs(iNSamples,dSampFreq);

    //
    // Create storage for data intervall
    //

    MatrixXd t_NSampleMat(m_iNChannels,iNSamples);

    if(m_iIntervallLengthFactor > 1)
        for(int i=0; i<m_iIntervallLengthFactor-1;++i){
            MatrixXd t_mat;
            while(!m_pBandPowerBuffer->pop(t_mat));
            if((t_mat.rows()==m_iNChannels)&&(t_mat.cols()==m_iNTimeSteps))
                t_NSampleMat.block(0,i*m_iNTimeSteps,m_iNChannels,m_iNTimeSteps) = t_mat;
            else
                qWarning() << "[BandPower::run] Data matrix has wrong format - not processing";
        }

    while(!isInterruptionRequested())
    {
        //Dispatch the inputs
        MatrixXd t_mat;
        while(!m_pBandPowerBuffer->pop(t_mat));

        m_qMutex.lock();

        if(t_NSampleMat.cols()/m_iNTimeSteps != m_iIntervallLengthFactor)
        {
            if (t_NSampleMat.cols()/m_iNTimeSteps > m_iIntervallLengthFactor) {
                t_NSampleMat = t_NSampleMat.rightCols(m_iIntervallLengthFactor*m_iNTimeSteps); //right-most m_iNTimeSteps columns don't matter here, since they are overwritten anyway
            } else {
                Eigen::MatrixXd t_NSampleMatTemp(t_NSampleMat.rows(),m_iIntervallLengthFactor*m_iNTimeSteps);
                t_NSampleMatTemp.leftCols(t_NSampleMat.cols()) = t_NSampleMat;
                for(int i=0; i < (m_iIntervallLengthFactor - t_NSampleMat.cols()/m_iNTimeSteps); i++)
                {
                    t_NSampleMatTemp.block(0,t_NSampleMat.cols() + i*m_iNTimeSteps, t_mat.rows(), t_mat.cols());
                    while(!m_pBandPowerBuffer->pop(t_mat));
                }
                t_NSampleMat = t_NSampleMatTemp;
            }
            iNSamples = m_iIntervallLengthFactor * m_iNTimeSteps;
            FFTFreqs = Spectral::calculateFFTFreqs(iNSamples,dSampFreq);
            qDebug() << "[BandPower::run] FFT resolution " << (FFTFreqs[1] - FFTFreqs[0]);
        }
        t_NSampleMat.rightCols(m_iNTimeSteps) = t_mat;
        MatrixXd t_SampleSubMat;

        if (m_pSelectedChannels.length() == 0)
            t_SampleSubMat = t_NSampleMat;
        else {
            if (m_pSelectedChannels.at(0) < t_NSampleMat.rows())
                t_SampleSubMat = t_NSampleMat.block(m_pSelectedChannels.at(0),0,1,t_NSampleMat.cols());
            else
            {
                t_SampleSubMat = t_NSampleMat;
                qDebug() << "m_pSelectedChannels" << m_pSelectedChannels.at(0) << "not in sample Matrix";
            }
            for (int i=1; i < m_pSelectedChannels.length(); ++i)
            {
                if (m_pSelectedChannels.at(i) < t_NSampleMat.rows()) {
                    t_SampleSubMat.conservativeResize(t_SampleSubMat.rows() + 1, NoChange);
                    t_SampleSubMat.row(i) = t_NSampleMat.row(m_pSelectedChannels.at(i));
                }
                else
                    qDebug() << "m_pSelectedChannels" << m_pSelectedChannels.at(i) << "not in sample Matrix";
            }
        }

        MatrixXd t_SampleMat_noav = Spectral::detrendData(t_SampleSubMat,m_iDetrendMethod);

        QVector<VectorXd> matSpectrum;
        //double stepwidth = 0;

        Eigen::MatrixXd bandpower(m_iBandPowerChannels*m_iBandPowerBins,1);

        if(m_sSpectrumMethod=="AR"){
            MatrixXcd ARSpectraWeights = Spectral::generateARSpectraWeights(m_dFreqMin/dSampFreq,m_dFreqMax/dSampFreq, m_iBandPowerBins, m_iEvalsAR, true);
            QVector<QPair<VectorXd, double>> ARCoeffs = Spectral::calculateARWeightsMEMMatrix(t_SampleMat_noav,m_iOrderAR,true);
            matSpectrum = Spectral::psdFromARSpectra(ARCoeffs,ARSpectraWeights,dSampFreq,true);
            //stepwidth = (m_dFreqMax - m_dFreqMin)/(matSpectrum.at(0).size() - 1); //check where the evaluation points lie
            //for (int i=0; i < matSpectrum.length(); ++i)
            //        meanbandpower(0,0) += 1e10*bandpowerFromSpectrumEntries(matSpectrum.at(i),stepwidth)/matSpectrum.length();
            if(m_iBandPowerChannels <= matSpectrum.length())
                qDebug() << "More channels selected than pre-defined! Only the first " << m_iBandPowerChannels << " are displayed!";
            for(int i=0; i<std::min(m_iBandPowerChannels,matSpectrum.length()); ++i){
                bandpower.block(i*m_iBandPowerBins,0,m_iBandPowerBins,1) = matSpectrum.at(i);
            }
            matSpectrum.clear();
        } else if(m_sSpectrumMethod=="FFT") {
            // Generate hanning window
            QPair<MatrixXd, VectorXd> tapers = Spectral::generateTapers(iNSamples, "hanning");
            MatrixXd matTaps = tapers.first;
            VectorXd vecTapWeights = tapers.second;

            // Compute Spectrum
            QVector<MatrixXcd> matTaperedSpectrum;
            matTaperedSpectrum = Spectral::computeTaperedSpectraMatrix(t_SampleMat_noav, matTaps, iNSamples, true);

            matSpectrum = Spectral::psdFromTaperedSpectra(matTaperedSpectrum, vecTapWeights, iNSamples, dSampFreq, false);

            // Select frequencies that fall within the band
            if(m_iBandPowerChannels <= matSpectrum.length())
                qDebug() << "[BandPower::run] More channels selected than pre-defined! Only the first " << m_iBandPowerChannels << " are displayed!";

            double binwidth = (m_dFreqMax-m_dFreqMin)/static_cast<double>(m_iBandPowerBins);

            if (binwidth < (FFTFreqs[1] - FFTFreqs[0]))
                qDebug() << "[BandPower::run] Selected bin width is smaller than FFT resolution";

            for(int i=0; i<std::min(m_iBandPowerChannels,matSpectrum.length()); ++i)
                for (int j=0; j<m_iBandPowerBins; ++j)
                {
                    bandpower(i*m_iBandPowerBins + j,1) = Spectral::bandpowerFromSpectrumEntriesOffset(FFTFreqs, matSpectrum.at(i), m_dFreqMin + j*binwidth, m_dFreqMin + (j+1)*binwidth);
                }
            matSpectrum.clear();
        }

        m_qMutex.unlock();

        //Send the data to the connected plugins and the online display

        m_pBandPowerOutput->data()->setValue(bandpower);

        qDebug() << "Power:" << bandpower(0,0);

        //move matrix entries one block to the left
        if(t_NSampleMat.cols()/m_iNTimeSteps > 1)
            for(int i=0; i<t_NSampleMat.cols()/m_iNTimeSteps-1; ++i)
                t_NSampleMat.block(0,i*m_iNTimeSteps,m_iNChannels,m_iNTimeSteps) = t_NSampleMat.block(0,(i+1)*m_iNTimeSteps,m_iNChannels,m_iNTimeSteps);
    }
}

//=============================================================================================================

void BandPower::showSensorSelectionWidget()
{
    if(m_pChannelSelectionView->isActiveWindow()) {
        m_pChannelSelectionView->hide();
    } else {
        m_pChannelSelectionView->activateWindow();
        m_pChannelSelectionView->show();
    }
}
