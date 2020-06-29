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
    QWidget* setupWidget = new QWidget();//widget is later distroyed by CentralWidget - so it has to be created everytime new
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

double BandPower::bandpowerFromSpectrumEntries(const Eigen::VectorXd &spectrumentries, double stepsize)
{
    double power;

    if (spectrumentries.size() % 2 == 0){
        Eigen::RowVectorXd weights = RowVectorXd::Ones(spectrumentries.size());
        for (int i=0; i<spectrumentries.size()/2; i++)
            weights(2*i+1) = 4.0;
        for (int i=0; i<spectrumentries.size()/2 - 1; i++)
            weights(2*i+2) = 2.0;
        power =  weights * spectrumentries;
        power *= (stepsize / 3.0);
    } else {
        // Use Simpson's rule for n-1 entries and trapezoid rule for the last intervall and average over 2 cases
        Eigen::RowVectorXd weights = RowVectorXd::Ones(spectrumentries.size() - 1);
        for (int i=0; i<spectrumentries.size()/2; i++)
            weights(2*i+1) = 4.0;
        for (int i=0; i<spectrumentries.size()/2 - 1; i++)
            weights(2*i+2) = 2.0;
        power = weights * spectrumentries.head(spectrumentries.size()-1);
        power = weights * spectrumentries.tail(spectrumentries.size()-1);
        power += stepsize/2.0 * (spectrumentries(0) + spectrumentries(1) + spectrumentries(spectrumentries.size()-2) + spectrumentries(spectrumentries.size()-1));
        power *= 0.5;
    }

    return power;
}

//=============================================================================================================

double BandPower::bandpowerFromSpectrumEntriesOffset(const Eigen::VectorXd &spectrumbins, const Eigen::VectorXd &spectrumentries, double minFreq, double maxFreq, double eps)
{
    if ((minFreq < spectrumbins(0)) || (maxFreq > spectrumbins(spectrumbins.size()-1)))
    {
            qDebug() << "Warning: BandPower::bandpowerFromSpectrumEntries - minFreq or maxFreq out of Spectrum";
    }

    double stepsize = spectrumbins(1) - spectrumbins(0);

    int minIndex, maxIndex;
    double deltaMin, deltaMax;

    Eigen::VectorXd spectrumbinsMin = spectrumbins;
    spectrumbinsMin.array() -= minFreq;
    spectrumbinsMin.cwiseAbs().minCoeff(&minIndex);
    deltaMin = minFreq - spectrumbins(minIndex);

    Eigen::VectorXd spectrumbinsMax = spectrumbins;
    spectrumbinsMax.array() -= maxFreq;
    spectrumbinsMax.cwiseAbs().minCoeff(&maxIndex);
    deltaMax = maxFreq - spectrumbins(maxIndex);

    if (std::abs(deltaMin) < eps && std::abs(deltaMax) < eps)
        // everything fits perfectly, nothing to do here
        return BandPower::bandpowerFromSpectrumEntries(spectrumentries,stepsize);
    else
    {
        if (deltaMin > 0)
            minIndex++;
        else
            deltaMin *= -1.0;

        if (deltaMax < 0)
        {
            maxIndex--;
            deltaMax *= -1.0;
        }
        double result = 0;

        int diffIndex = maxIndex - minIndex;

        if (diffIndex < 0) // minFreq and maxFreq in the same intervall (should not happen) apply trapezoid rule with linear interpolation
        {
            result = ((spectrumentries(minIndex) + (spectrumentries(minIndex-1)-spectrumentries(minIndex))*deltaMin/stepsize) +
                      (spectrumentries(maxIndex) + (spectrumentries(maxIndex+1)-spectrumentries(maxIndex))*deltaMax/stepsize))*(maxFreq-minFreq)/2.0;
        } else
        {
        // add the contributions from the incomplete bins
        if (minIndex > 0)
            result += (spectrumentries(minIndex) + (spectrumentries(minIndex-1)-spectrumentries(minIndex))*deltaMin/(2.0*stepsize))*deltaMin;

        if (maxIndex < spectrumentries.size() - 1)
            result += (spectrumentries(maxIndex) + (spectrumentries(maxIndex+1)-spectrumentries(maxIndex))*deltaMax/(2.0*stepsize))*deltaMax;

        if (diffIndex > 1) // use Simpsons's rule for more than two points
            result += BandPower::bandpowerFromSpectrumEntries(spectrumentries.segment(minIndex,maxIndex - minIndex + 1),stepsize);
        else if (diffIndex > 0) // use trapezoid rule for two points
            result += (spectrumentries(minIndex) + spectrumentries(maxIndex))*stepsize/2.0;
        }

        return result;
    }
}

//=============================================================================================================

Eigen::MatrixXd BandPower::detrendData(const Eigen::MatrixXd &data, int method)
{
    switch(method)
    {
    case 0:
        return data;
    case 1:
    {
        Eigen::MatrixXd outputData = data.colwise() - data.rowwise().mean();
        return outputData;
    }
    case 2:
        return data;
    default:
        qDebug() << "BandPower::detrendData - unknown detrend method" << method << "Returning data";
        return data;
    }
}

//=============================================================================================================

Eigen::MatrixXd BandPower::linearDetrend(const Eigen::MatrixXd &data)
{
    int n = data.cols();

    Eigen::VectorXd linspaced = VectorXd::LinSpaced(1,0,n-1);

    Eigen::VectorXd xy = data * linspaced;

    double  linspaced_sum = n * (n - 1) / 2,
            linspaced_sum2 = (n - 1) * n * (2 * n - 1) / 6;

    Eigen::VectorXd data_rowwise_sum = data.rowwise().sum();

    double denom_a = n*linspaced_sum2-linspaced_sum*linspaced_sum;

    Eigen::VectorXd a = (n*xy - linspaced_sum*data_rowwise_sum)/denom_a;
    Eigen::VectorXd b = (data_rowwise_sum - linspaced_sum*a)/n;

    Eigen::MatrixXd data_noav = data.colwise() - b;
    data_noav -= a*linspaced.transpose();

    return data_noav;
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
            //m_pFiffInfo = FIFFLIB::FiffInfo::SPtr(new FIFFLIB::FiffInfo);

            //Init output - Uncomment this if you also uncommented the m_pDummyOutput in the constructor above
            //m_pBandPowerOutput->data()->initFromFiffInfo(m_pFiffInfo);

            m_dDataSampFreq = m_pFiffInfo->sfreq;

            m_pFiffInfo->filename = "";
            m_pFiffInfo->bads.clear();
            m_pFiffInfo->nchan = 2;

            QList<FiffChInfo> fakeChList;
            FiffChInfo fakeCh;
            fakeCh.ch_name = "AR";
            fakeCh.kind = 502;
            fakeCh.range = -1;
            //fakeCh.setMinValue(0.5e-20);
            //fakeCh.setMaxValue(2e-20);
            fakeCh.unit = -1;
            fakeChList.append(fakeCh);
            fakeCh.ch_name = "FFT";
            fakeChList.append(fakeCh);
            m_pFiffInfo->chs = fakeChList;

            QStringList fakeChNames = {"AR", "FFT"};
            m_pFiffInfo->ch_names = fakeChNames;

            m_pFiffInfo->file_id = FIFFLIB::FiffId::new_file_id(); //check if necessary

            m_pFiffInfo->sfreq = pRTMSA->info()->sfreq/pRTMSA->getMultiSampleArray().first().cols();

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
    if(m_pFiffInfo) {
        QList<QWidget*> plControlWidgets;

        BandPowerSettingsView* pBandPowerSettingsView = new BandPowerSettingsView(QString("MNESCAN/%1/").arg(this->getName()),m_dDataSampFreq,
                                                                                  m_dFreqMin, m_dFreqMax, m_sSpectrumMethod, m_iIntervallLengthFactor);
        pBandPowerSettingsView->setObjectName("group_tab_Settings_General");
        plControlWidgets.append(pBandPowerSettingsView);

        connect(pBandPowerSettingsView,&BandPowerSettingsView::changeMinMax,this,&BandPower::changeMinMaxFrequency);
        connect(pBandPowerSettingsView,&BandPowerSettingsView::changeMethod,this,&BandPower::changeSpectrumMethod);
        //connect(pBandPowerSettingsView, &BandPowerSettingsView::changeIntervallLength,this,&BandPower::changeIntervallLengthFactor);
        connect(pBandPowerSettingsView,&BandPowerSettingsView::changeDetrend,this,&BandPower::changeDetrendMethod);

        ARSettingsView* pARSettingsView = new ARSettingsView(QString("MNESCAN/%1/").arg(this->getName()));

        pARSettingsView->setObjectName("group_tab_Settings_AR");
        plControlWidgets.append(pARSettingsView);

        connect(pARSettingsView, &ARSettingsView::changeAROrder,this,&BandPower::changeAROrder);
        connect(pARSettingsView, &ARSettingsView::changeARNumEvaluationPoints,this,&BandPower::changeARNumEvaluationPoints);

        //Channelselectionwidget

        //add button
        m_pActionSelectSensors = new QAction(QIcon(":/images/selectSensors.png"), tr("Show the channel selection view"),this);
        m_pActionSelectSensors->setToolTip(tr("Show the channel selection view"));
        connect(m_pActionSelectSensors.data(), &QAction::triggered,
                this, &RealTimeMultiSampleArrayWidget::showSensorSelectionWidget);
        addDisplayAction(m_pActionSelectSensors);
        m_pActionSelectSensors->setVisible(true);

        //Init channel selection manager
        m_pChannelInfoModel = ChannelInfoModel::SPtr::create(m_pFiffInfo,
                                                             this);

        m_pChannelSelectionView = ChannelSelectionView::SPtr::create(QString("RTMSAW/%1").arg(sRTMSAWName),
                                                                     this,
                                                                     m_pChannelInfoModel,
                                                                     Qt::Window);
        m_pChannelSelectionView->setWindowTitle(tr(QString("%1: Channel Selection Window").arg(sRTMSAWName).toUtf8()));

        connect(m_pChannelSelectionView.data(), &ChannelSelectionView::loadedLayoutMap,
                m_pChannelInfoModel.data(), &ChannelInfoModel::layoutChanged);

        connect(m_pChannelInfoModel.data(), &ChannelInfoModel::channelsMappedToLayout,
                m_pChannelSelectionView.data(), &ChannelSelectionView::setCurrentlyMappedFiffChannels);

        connect(m_pChannelSelectionView.data(), &ChannelSelectionView::showSelectedChannelsOnly,
                m_pChannelDataView.data(), &RtFiffRawView::showSelectedChannelsOnly);

        connect(m_pChannelDataView.data(), &RtFiffRawView::channelMarkingChanged,
                m_pChannelSelectionView.data(), &ChannelSelectionView::updateBadChannels);

        connect(m_pChannelDataView.data(), &RtFiffRawView::selectedChannelsChanged,
                m_pChannelSelectionView.data(), &ChannelSelectionView::setUserSelection);

        connect(m_pChannelDataView.data(), &RtFiffRawView::selectedChannelsResetted,
                m_pChannelSelectionView.data(), &ChannelSelectionView::resetUserSelection);

        m_pChannelInfoModel->layoutChanged(m_pChannelSelectionView->getLayoutMap());

        emit pluginControlWidgetsChanged(plControlWidgets, this->getName());
    }
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
    int iNfft = static_cast<int> (iNSamples / m_iResolutionFactorFFT);

    //
    // Create storage for data intervall
    //

    MatrixXd t_NSampleMat(m_iNChannels,iNSamples);

    if(m_iIntervallLengthFactor > 1)
        for(int i=0; i<m_iIntervallLengthFactor-1;++i){
            MatrixXd t_mat;
            while(!m_pBandPowerBuffer->pop(t_mat));
            //check for size!
            t_NSampleMat.block(0,i*m_iNTimeSteps,m_iNChannels,m_iNTimeSteps) = t_mat;
        }

    while(!isInterruptionRequested())
    {
        //Dispatch the inputs
        MatrixXd t_mat;
        while(!m_pBandPowerBuffer->pop(t_mat));

        m_qMutex.lock();

        t_NSampleMat.rightCols(m_iNTimeSteps) = t_mat;

        MatrixXd t_NSampleMat_noav = detrendData(t_NSampleMat,m_iDetrendMethod);

        QVector<VectorXd> matSpectrum;
        double stepwidth = 0;

//        if(m_sSpectrumMethod=="AR"){
            MatrixXcd ARSpectraWeights = Spectral::generateARSpectraWeights(m_dFreqMin/dSampFreq,m_dFreqMax/dSampFreq, m_iEvalsAR, 1, true);

            QVector<QPair<VectorXd, double>> ARCoeffs = Spectral::calculateARWeightsMEMMatrix(t_NSampleMat_noav,m_iOrderAR,true);

            matSpectrum = Spectral::psdFromARSpectra(ARCoeffs,ARSpectraWeights,dSampFreq,true); // has to be extended for multiple rows!

            stepwidth = (m_dFreqMax - m_dFreqMin)/(matSpectrum.at(0).size() - 1); //check where the evaluation points lie

            MatrixXd meanbandpower(2,1);
            meanbandpower(0,0) = 1e12*bandpowerFromSpectrumEntries(matSpectrum.at(31),stepwidth);

            matSpectrum.clear();

//        } else if(m_sSpectrumMethod=="FFT") {
            // Generate hanning window
            QPair<MatrixXd, VectorXd> tapers = Spectral::generateTapers(iNSamples, "hanning");
            MatrixXd matTaps = tapers.first;
            VectorXd vecTapWeights = tapers.second;

            // Compute Spectrum
            QVector<MatrixXcd> matTaperedSpectrum;
            matTaperedSpectrum = Spectral::computeTaperedSpectraMatrix(t_NSampleMat_noav, matTaps, iNfft, true);

            matSpectrum = Spectral::psdFromTaperedSpectra(matTaperedSpectrum, vecTapWeights, iNfft, dSampFreq, false);

            // Select frequencies that fall within the band
            Eigen::VectorXd FFTFreqs = Spectral::calculateFFTFreqs(iNfft,dSampFreq);

//        }

            qDebug() << FFTFreqs(0) << FFTFreqs(FFTFreqs.size()-1) << m_dFreqMin << m_dFreqMax;

        meanbandpower(1,0) = 1e12*bandpowerFromSpectrumEntriesOffset(FFTFreqs,matSpectrum.at(31),m_dFreqMin,m_dFreqMax);

        m_qMutex.unlock();

        //Send the data to the connected plugins and the online display
        m_pBandPowerOutput->data()->setValue(meanbandpower);

        qDebug() << "Power:" << meanbandpower(0,0) << meanbandpower(1,0);

        //move matrix entries one block to the left
        if(m_iIntervallLengthFactor > 1)
            for(int i=0; i<m_iIntervallLengthFactor-1; ++i)
                t_NSampleMat.block(0,i*m_iNTimeSteps,m_iNChannels,m_iNTimeSteps) = t_NSampleMat.block(0,(i+1)*m_iNTimeSteps,m_iNChannels,m_iNTimeSteps);
    }
}
