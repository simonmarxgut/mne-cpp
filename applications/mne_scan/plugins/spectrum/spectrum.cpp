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

//#include "FormFiles/spectrumsetupwidget.h"

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace SPECTRUMPLUGIN;
using namespace SCMEASLIB;
using namespace UTILSLIB;
using namespace IOBUFFER;
//using namespace DISPLIB;
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
    , m_iBinsAR(256)
    , m_iEvalPerBinAR(5)
    , m_sSpectrumMethod("AR")
    , m_iNTimeSteps(-1)
    , m_iNChannels(-1)
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
    //QWidget* setupWidget = new SpectrumSetupWidget(this);//widget is later distroyed by CentralWidget - so it has to be created everytime new
    QWidget* setupWidget = new QWidget();//widget is later distroyed by CentralWidget - so it has to be created everytime new
    return setupWidget;
}

//=============================================================================================================

void Spectrum::changeIntervallLengthFactor(quint32 lengthfactor)
{
    m_iIntervallLengthFactor = lengthfactor;
    //if(m_pRtCov) {
    //    m_pRtCov->setSamples(m_iEstimationSamples);
    //}
}

//=============================================================================================================

void Spectrum::changeResolutionFactor(quint32 resolutionfactor)
{
    m_iResolutionFactorFFT = resolutionfactor;
    //if(m_pRtCov) {
    //    m_pRtCov->setSamples(m_iEstimationSamples);
    //}
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
            m_pFiffInfo = pRTMSA->info();

            //Init output - Uncomment this if you also uncommented the m_pDummyOutput in the constructor above
            m_pSpectrumOutput->data()->initFromFiffInfo(m_pFiffInfo);
            // m_pSpectrumOutput->data()->setMultiArraySize(1);
            m_pSpectrumOutput->data()->setVisibility(true);
        }

        // Check if data is present
        if(pRTMSA->getMultiSampleArray().size() > 0) {
            //Init widgets
            if(m_iNChannels == -1) {
                m_iNChannels = pRTMSA->getMultiSampleArray().first().rows();
                //initPluginControlWidgets();
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


//=============================================================================================================

void Spectrum::run()
{
    // Length of block
    int iNSamples = m_iNTimeSteps * m_iIntervallLengthFactor;

    // Sampling frequency
    double dSampFreq = m_pFiffInfo->sfreq;

    // Determine spectrum resolution - calculation only for FFT, but we keep it for both methods here to keep the resolution
    int iNfft = static_cast<int> (iNSamples / m_iResolutionFactorFFT);

    //
    // Create storage for data intervall
    //

    MatrixXd t_NSampleMat(m_iNChannels,iNSamples);

    if(m_iIntervallLengthFactor > 1)
        for(int i=0; i<m_iIntervallLengthFactor-1;++i){
            MatrixXd t_mat;
            while(!m_pSpectrumBuffer->pop(t_mat));
            //check for size!
            t_NSampleMat.block(0,i*m_iNTimeSteps,m_iNChannels,m_iNTimeSteps) = t_mat;
        }

    qDebug() << m_sSpectrumMethod;

    while(!isInterruptionRequested())
    {
        Eigen::VectorXd FFTFreqs = Spectral::calculateFFTFreqs(iNfft,dSampFreq);

        //Dispatch the inputs
        MatrixXd t_mat;
        while(!m_pSpectrumBuffer->pop(t_mat));

        t_NSampleMat.rightCols(m_iNTimeSteps) = t_mat;

        Eigen::MatrixXd t_NSampleMat_noav = t_NSampleMat;

        for (int i = 0; i < t_NSampleMat_noav.rows(); i++)
            t_NSampleMat_noav.row(i).array() = t_NSampleMat_noav.row(i).array() - t_NSampleMat_noav.row(i).mean();


        QVector<VectorXd> matSpectrum;

        if(m_sSpectrumMethod=="AR"){
            MatrixXcd ARSpectraWeights = Spectral::generateARSpectraWeights(FFTFreqs.minCoeff()/dSampFreq,FFTFreqs.maxCoeff()/dSampFreq, FFTFreqs.size(), m_iEvalPerBinAR, true);
            QVector<QPair<VectorXd, double>> ARCoeffs = Spectral::calculateARWeightsMEMMatrix(t_NSampleMat_noav,m_iOrderAR,true);
            matSpectrum = Spectral::psdFromARSpectra(ARCoeffs,ARSpectraWeights,dSampFreq,true); // has to be extended for multiple rows!

        } else if(m_sSpectrumMethod=="FFT") {
            // Generate hanning window
            QPair<MatrixXd, VectorXd> tapers = Spectral::generateTapers(iNSamples, "hanning");
            MatrixXd matTaps = tapers.first;
            VectorXd vecTapWeights = tapers.second;

            // Compute Spectrum
            QVector<MatrixXcd> matTaperedSpectrum;
            matTaperedSpectrum = Spectral::computeTaperedSpectraMatrix(t_NSampleMat_noav, matTaps, iNfft, true);

            matSpectrum = Spectral::psdFromTaperedSpectra(matTaperedSpectrum, vecTapWeights, iNfft, dSampFreq, true);

        }

        qDebug() << "matSpectrum" << matSpectrum.length() << " " << matSpectrum.at(0).size();

        //Compute Mean Spectrum
        MatrixXd matSpectrumMean = MatrixXd::Zero(matSpectrum.at(0).size(),1);

/*        for (int i = 0; i < 60; ++i) {
            matSpectrumMean.array() += matSpectrum.at(i).array().abs().log();
        }

        matSpectrumMean /= 60.0;*/

        matSpectrumMean = matSpectrum.at(31);

        matSpectrumMean.array() = matSpectrumMean.array() + 1;

        matSpectrumMean.array() = matSpectrumMean.array().log();

        int maxIndex, wasteIndex;

        matSpectrumMean.maxCoeff(&maxIndex, &wasteIndex);

        qDebug() << FFTFreqs(1) - FFTFreqs(0) << FFTFreqs(FFTFreqs.size() - 1) << FFTFreqs(maxIndex) << " Hz";

        //MatrixXd matSpectrumOut = MatrixXd::Zero(matTapSpectrumSeed.size(),matTapSpectrumSeed.at(0).cols());

        //for (int i = 0; i < matTapSpectrumSeed.size(); ++i) {
        //    matSpectrumOut.row(i) = matTapSpectrumSeed.at(i).array().abs();
        //}

        qDebug() << "matSpectrumMean" << matSpectrumMean.rows() << " " << matSpectrumMean.cols();

        //matSpectrumMean.transposeInPlace();

        //Send the data to the connected plugins and the online display
        m_pSpectrumOutput->data()->setValue(matSpectrumMean);

        //move matrix entries one block to the left
        if(m_iIntervallLengthFactor > 1)
            for(int i=0; i<m_iIntervallLengthFactor-1; ++i)
                t_NSampleMat.block(0,i*m_iNTimeSteps,m_iNChannels,m_iNTimeSteps) = t_NSampleMat.block(0,(i+1)*m_iNTimeSteps,m_iNChannels,m_iNTimeSteps);

    }
}
