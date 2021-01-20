//=============================================================================================================
/**
 * @file     spectrum.cpp
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
 * @brief    Definition of the Spectrum class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "spectrum.h"

#include "FormFiles/spectrumsetupwidget.h"

//#include "../bandpower/bandpower.h"

#include "disp/viewers/bandpowersettingsview.h"
#include "disp/viewers/arsettingsview.h"

#include <disp/viewers/channelselectionview.h>
#include <disp/viewers/helpers/channelinfomodel.h>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace SPECTRUMPLUGIN;
//using namespace BANDPOWERPLUGIN;
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

Spectrum::Spectrum()
    : m_bProcessData(false)
    , m_pSpectrumInput(NULL)
    , m_pSpectrumOutput(NULL)
    , m_pSpectrumBuffer(CircularBuffer<Eigen::MatrixXd>::SPtr())
    , m_iIntervallLengthFactor(5)
    , m_iResolutionFactorFFT(1)
    , m_iOrderAR(15)
    , m_iEvalsAR(5)
    , m_sSpectrumMethod("AR")
    , m_iNTimeSteps(-1)
    , m_iNChannels(-1)
    , m_iSpectrumBins(10)
    , m_iDetrendMethod(0)
    , m_dFreqMin(8)
    , m_dFreqMax(30)
{
}

//=============================================================================================================

Spectrum::~Spectrum()
{
    if(this->isRunning())
        stop();
}

//=============================================================================================================

QSharedPointer<IPlugin> Spectrum::clone() const
{
    QSharedPointer<Spectrum> pSpectrumClone(new Spectrum);
    return pSpectrumClone;
}

//=============================================================================================================

void Spectrum::init()
{
    // Input
    m_pSpectrumInput = PluginInputData<RealTimeMultiSampleArray>::create(this, "SpectrumIn", "Spectrum input data");
    connect(m_pSpectrumInput.data(), &PluginInputConnector::notify, this, &Spectrum::update, Qt::DirectConnection);
    m_inputConnectors.append(m_pSpectrumInput);

    // Output - Uncomment this if you don't want to send processed data (in form of a matrix) to other plugins.
    // Also, this output stream will generate an online display in your plugin
    m_pSpectrumOutput = PluginOutputData<RealTimeSpectrum>::create(this, "SpectrumOut", "Spectrum output data");
    // m_pSpectrumOutput = PluginOutputData<RealTimeMultiSampleArray>::create(this, "SpectrumOut", "Spectrum output data");

    m_outputConnectors.append(m_pSpectrumOutput);

    //Delete Buffer - will be initailzed with first incoming data
    if(!m_pSpectrumBuffer.isNull())
        m_pSpectrumBuffer = CircularBuffer<Eigen::MatrixXd>::SPtr();

    //add button
    m_pActionSelectSensors = new QAction(QIcon(":/images/selectSensors.png"), tr("Show the channel selection view for spectrum computation"),this);
    m_pActionSelectSensors->setToolTip(tr("Show the channel selection view for spectrum computation"));
    connect(m_pActionSelectSensors.data(), &QAction::triggered,
            this, &Spectrum::showSensorSelectionWidget);
    addPluginAction(m_pActionSelectSensors);
    m_pActionSelectSensors->setVisible(true);
}

//=============================================================================================================

void Spectrum::unload()
{
}

//=============================================================================================================

bool Spectrum::start()
{
//    //Check if the thread is already or still running. This can happen if the start button is pressed immediately after the stop button was pressed. In this case the stopping process is not finished yet but the start process is initiated.
//    if(this->isRunning())
//        QThread::wait();

    return true;
}


//=============================================================================================================

bool Spectrum::stop()
{
    requestInterruption();
    wait(500);

    //m_pSpectrumOutput->data()->clear();
    m_pSpectrumBuffer->clear();

    return true;
}

//=============================================================================================================

IPlugin::PluginType Spectrum::getType() const
{
    return _IAlgorithm;
}

//=============================================================================================================

QString Spectrum::getName() const
{
    return "Spectrum";
}

//=============================================================================================================

QWidget* Spectrum::setupWidget()
{
    QWidget* setupWidget = new SpectrumSetupWidget(this);//widget is later distroyed by CentralWidget - so it has to be created everytime new
    //QWidget* setupWidget = new QWidget();//widget is later distroyed by CentralWidget - so it has to be created everytime new
    return setupWidget;
}

//=============================================================================================================

void Spectrum::changeMinMaxFrequency(double minFrequency, double maxFrequency)
{
    QMutexLocker locker(&m_qMutex);

    m_dFreqMin = minFrequency;
    m_dFreqMax = maxFrequency;

    calcFrequencies();
}

//=============================================================================================================

void Spectrum::changeSpectrumMethod(QString method)
{
    QMutexLocker locker(&m_qMutex);

    m_sSpectrumMethod = method;
}

//=============================================================================================================

void Spectrum::changeDetrendMethod(quint32 method)
{
    QMutexLocker locker(&m_qMutex);

    m_iDetrendMethod = method;
}

//=============================================================================================================

void Spectrum::changeIntervallLengthFactor(quint32 lengthfactor)
{
    QMutexLocker locker(&m_qMutex);

    m_iIntervallLengthFactor = lengthfactor;
}

//=============================================================================================================

void Spectrum::changeResolutionFactor(quint32 resolutionfactor)
{
    QMutexLocker locker(&m_qMutex);

    m_iResolutionFactorFFT = resolutionfactor;
}


//=============================================================================================================

void Spectrum::changeSpectrumBins(quint32 bins)
{
    QMutexLocker locker(&m_qMutex);

    m_iSpectrumBins = bins;

    calcFrequencies();

    qDebug() << "m_iBandPowerBins " << m_iSpectrumBins;
}

//=============================================================================================================

void Spectrum::changeAROrder(quint32 order)
{
    QMutexLocker locker(&m_qMutex);

    if(order < 1){
        qDebug() << "AR order < 1 not allowed.";
        m_iOrderAR = 1;
    } else
        m_iOrderAR = order;
}

//=============================================================================================================

void Spectrum::changeARNumEvaluationPoints(quint32 evaluationpoints)
{
    QMutexLocker locker(&m_qMutex);

    m_iEvalsAR = evaluationpoints;
}

//=============================================================================================================

void Spectrum::evaluateSelectedChannelsOnly(const QStringList &selectedChannels)
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

void Spectrum::calcFrequencies()
{
    m_pFrequencies = Eigen::VectorXd::Zero(m_iSpectrumBins);
    double delta = (m_dFreqMax - m_dFreqMin)/(m_iSpectrumBins - 1);
    m_pFrequencies[0] = m_dFreqMin;
    for (int i=1; i < m_iSpectrumBins; i++)
        m_pFrequencies[i] = m_pFrequencies[i-1] + delta;
}

//=============================================================================================================

void Spectrum::update(SCMEASLIB::Measurement::SPtr pMeasurement)
{
    QSharedPointer<RealTimeMultiSampleArray> pRTMSA = pMeasurement.dynamicCast<RealTimeMultiSampleArray>();

    if(pRTMSA) {
        //Check if buffer initialized
        if(!m_pSpectrumBuffer) {
            m_pSpectrumBuffer = CircularBuffer<Eigen::MatrixXd>::SPtr(new CircularBuffer<Eigen::MatrixXd>(64));
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
            m_pFiffInfo->nchan = 1;

            QList<FiffChInfo> fakeChList;
            FiffChInfo fakeCh;
            fakeCh.ch_name = "Spectrum";
            fakeCh.kind = 502;
            fakeCh.range = -1;
            //fakeCh.setMinValue(0.5e-20);
            //fakeCh.setMaxValue(2e-20);
            fakeCh.unit = -1;
            fakeChList.append(fakeCh);

            QStringList fakeChNames;
            fakeChNames.append(fakeCh.ch_name);

            m_pFiffInfo->chs = fakeChList;
            m_pFiffInfo->ch_names = fakeChNames;
            m_pFiffInfo->file_id = FIFFLIB::FiffId::new_file_id(); //check if necessary
            m_pFiffInfo->sfreq = pRTMSA->info()->sfreq/pRTMSA->getMultiSampleArray().first().cols();

            //Init output - Uncomment this if you also uncommented the m_pDummyOutput in the constructor above
            m_pSpectrumOutput->data()->initFromFiffInfo(m_pFiffInfo);
            // m_pSpectrumOutput->data()->setMultiArraySize(1);
            m_pSpectrumOutput->data()->setVisibility(true);
            m_pSpectrumOutput->data()->setName(this->getName());
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
                while(!m_pSpectrumBuffer->push(pRTMSA->getMultiSampleArray()[i])) {
                    //Do nothing until the circular buffer is ready to accept new data again
                }
            }
        }
    }
}

void Spectrum::initPluginControlWidgets()
{
    QList<QWidget*> plControlWidgets;

    BandPowerSettingsView* pSpectrumSettingsView = new BandPowerSettingsView(QString("MNESCAN/%1/").arg(this->getName()),m_dDataSampFreq,
                                                                              m_dFreqMin, m_dFreqMax, m_sSpectrumMethod, m_iIntervallLengthFactor,
                                                                              -1, m_iSpectrumBins, m_iDetrendMethod, true);
    pSpectrumSettingsView->setObjectName("group_tab_Settings_General");
    plControlWidgets.append(pSpectrumSettingsView);

    connect(pSpectrumSettingsView,&BandPowerSettingsView::changeMinMax,this,&Spectrum::changeMinMaxFrequency);
    connect(pSpectrumSettingsView,&BandPowerSettingsView::changeMethod,this,&Spectrum::changeSpectrumMethod);
    connect(pSpectrumSettingsView,&BandPowerSettingsView::changeIntervallLength,this,&Spectrum::changeIntervallLengthFactor);
    connect(pSpectrumSettingsView,&BandPowerSettingsView::changeDetrend,this,&Spectrum::changeDetrendMethod);
    connect(pSpectrumSettingsView,&BandPowerSettingsView::changeBins,this,&Spectrum::changeSpectrumBins);

    ARSettingsView* pARSettingsView = new ARSettingsView(QString("MNESCAN/%1/").arg(this->getName()));

    pARSettingsView->setObjectName("group_tab_Settings_AR");
    plControlWidgets.append(pARSettingsView);

    connect(pARSettingsView, &ARSettingsView::changeAROrder,this,&Spectrum::changeAROrder);
    connect(pARSettingsView, &ARSettingsView::changeARNumEvaluationPoints,this,&Spectrum::changeARNumEvaluationPoints);

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

    connect(m_pChannelSelectionView.data(), &ChannelSelectionView::showSelectedChannelsOnly,
            this, &Spectrum::evaluateSelectedChannelsOnly);

    m_pChannelInfoModel->layoutChanged(m_pChannelSelectionView->getLayoutMap());

    emit pluginControlWidgetsChanged(plControlWidgets, this->getName());
}


//=============================================================================================================

void Spectrum::run()
{
    m_qMutex.lock();
    // Length of block
    int iNSamples = m_iNTimeSteps * m_iIntervallLengthFactor;

    // Sampling frequency
    double dSampFreq = m_dDataSampFreq;
    qDebug() << "Sampling Frequency" << dSampFreq;

    // calc frequency vector
    calcFrequencies();

    qDebug() << "Frequencies" << m_dFreqMin << m_dFreqMax;

    // Determine spectrum resolution - calculation only for FFT, but we keep it for both methods here to keep the resolution
    Eigen::VectorXd FFTFreqs = Spectral::calculateFFTFreqs(iNSamples,dSampFreq);

    //
    // Create storage for data intervall
    //

    MatrixXd t_NSampleMat = Eigen::MatrixXd::Zero(m_iNChannels,iNSamples);

    if(m_iIntervallLengthFactor > 1)
        for(int i=0; i<m_iIntervallLengthFactor-1;++i){
            MatrixXd t_mat;
            while(!m_pSpectrumBuffer->pop(t_mat));
            if((t_mat.rows()==m_iNChannels)&&(t_mat.cols()==m_iNTimeSteps))
                t_NSampleMat.block(0,i*m_iNTimeSteps,m_iNChannels,m_iNTimeSteps) = t_mat;
            else
                qWarning() << "[Spectrum::run] Data matrix has wrong format - not processing";
        }
    m_qMutex.unlock();

    while(!isInterruptionRequested())
    {
        //Eigen::VectorXd FFTFreqs = Spectral::calculateFFTFreqs(iNfft,dSampFreq);

        //Dispatch the inputs
        MatrixXd t_mat;
        while(!m_pSpectrumBuffer->pop(t_mat));

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
                    while(!m_pSpectrumBuffer->pop(t_mat));
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

        Eigen::MatrixXd t_NSampleMat_noav = Spectral::detrendData(t_SampleSubMat,m_iDetrendMethod);

        QVector<VectorXd> matSpectrum;

        if(m_sSpectrumMethod=="AR"){
            MatrixXcd ARSpectraWeights = Spectral::generateARSpectraWeights(m_dFreqMin/dSampFreq,m_dFreqMax/dSampFreq, m_iSpectrumBins, m_iEvalsAR, true);
            QVector<QPair<VectorXd, double>> ARCoeffs = Spectral::calculateARWeightsMEMMatrix(t_NSampleMat_noav,m_iOrderAR,true);
            matSpectrum = Spectral::psdFromARSpectra(ARCoeffs,ARSpectraWeights,dSampFreq,true); // has to be extended for multiple rows!
        } else if(m_sSpectrumMethod=="FFT") {
            // Generate hanning window
            QPair<MatrixXd, VectorXd> tapers = Spectral::generateTapers(iNSamples, "hanning");
            MatrixXd matTaps = tapers.first;
            VectorXd vecTapWeights = tapers.second;

            // Compute Spectrum
            QVector<MatrixXcd> matTaperedSpectrum;
            matTaperedSpectrum = Spectral::computeTaperedSpectraMatrix(t_NSampleMat_noav, matTaps, iNSamples, true);

            QVector<VectorXd> matSpectrumFull = Spectral::psdFromTaperedSpectra(matTaperedSpectrum, vecTapWeights, iNSamples, dSampFreq, false);

            // Select frequencies that fall within the band
            double binwidth = (m_dFreqMax-m_dFreqMin)/static_cast<double>(m_iSpectrumBins);

            if (binwidth < (FFTFreqs[1] - FFTFreqs[0]))
                qDebug() << "[BandPower::run] Selected bin width is smaller than FFT resolution";

            for(int i=0; i<matSpectrumFull.length(); ++i)
                for (int j=0; j<m_iSpectrumBins; ++j)
                {
                    VectorXd powerTemp = VectorXd::Zero(m_iSpectrumBins);
                    powerTemp(j) = Spectral::bandpowerFromSpectrumEntriesOffset(FFTFreqs, matSpectrumFull.at(i), m_dFreqMin + j*binwidth, m_dFreqMin + (j+1)*binwidth);
                    matSpectrum.append(powerTemp);
                }
        }

        //Compute Mean Spectrum
        VectorXd matSpectrumMean = std::accumulate(matSpectrum.begin(),matSpectrum.end(),Eigen::VectorXd::Zero(matSpectrum.at(0).rows()).eval());
        //MatrixXd matSpectrumMean = matSpectrum.at(0);
        matSpectrumMean *= 1.0/matSpectrum.length();
        //matSpectrumMean.array().log();

        if(matSpectrumMean.rows() > 2)
            qDebug() << matSpectrumMean(0,0) << matSpectrumMean(1,0) << matSpectrumMean(2,0) << matSpectrumMean(3,0);

        MatrixXd matOut = MatrixXd::Zero(matSpectrumMean.rows(),2);

        matOut.col(0) = matSpectrumMean;
        matOut.col(1) = m_pFrequencies;

        m_qMutex.unlock();

        //Send the data to the connected plugins and the online display
        m_pSpectrumOutput->data()->setValue(matOut);

        //move matrix entries one block to the left
        if(t_NSampleMat.cols()/m_iNTimeSteps > 1)
            for(int i=0; i<t_NSampleMat.cols()/m_iNTimeSteps-1; ++i)
                t_NSampleMat.block(0,i*m_iNTimeSteps,t_NSampleMat.rows(),m_iNTimeSteps) = t_NSampleMat.block(0,(i+1)*m_iNTimeSteps,t_NSampleMat.rows(),m_iNTimeSteps);

    }
}


//=============================================================================================================

void Spectrum::showSensorSelectionWidget()
{
    if(m_pChannelSelectionView->isActiveWindow()) {
        m_pChannelSelectionView->hide();
    } else {
        m_pChannelSelectionView->activateWindow();
        m_pChannelSelectionView->show();
    }
}
