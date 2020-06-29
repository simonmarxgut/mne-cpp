//=============================================================================================================
/**
 * @file     spectral.cpp
 * @author   Daniel Strohmeier <Daniel.Strohmeier@tu-ilmenau.de>;
 *           Lorenz Esch <lesch@mgh.harvard.edu>
 * @since    0.1.0
 * @date     March, 2018
 *
 * @section  LICENSE
 *
 * Copyright (C) 2018, Daniel Strohmeier, Lorenz Esch. All rights reserved.
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
 * @note Notes:
 * - Some of this code was adapted from mne-python (https://martinos.org/mne) with permission from Alexandre Gramfort.
 * - This code is prepared for adding spectral estimation with multitapers, which is, however not yet supported.
 * - This code only allows FFT based spectral estimation. Time-frequency transforms are not yet supported.
 *
 * @brief    Declaration of Spectral class.
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "spectral.h"
#include "math.h"

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <unsupported/Eigen/FFT>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QtMath>
#include <QtConcurrent>
#include <QVector>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace UTILSLIB;
using namespace Eigen;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

MatrixXcd Spectral::computeTaperedSpectraRow(const RowVectorXd &vecData,
                                             const MatrixXd &matTaper,
                                             int iNfft)
{
    //qDebug() << "Spectral::computeTaperedSpectra Matrixwise";
    FFT<double> fft;
    fft.SetFlag(fft.HalfSpectrum);

    //Check inputs
    if (vecData.cols() != matTaper.cols() || iNfft < vecData.cols()) {
        return MatrixXcd();
    }

    //FFT for freq domain returning the half spectrum
    RowVectorXd vecInputFFT;
    RowVectorXcd vecTmpFreq;
    MatrixXcd matTapSpectrum(matTaper.rows(), int(floor(iNfft / 2.0)) + 1);
    for (int i=0; i < matTaper.rows(); i++) {
        vecInputFFT = vecData.cwiseProduct(matTaper.row(i));
        fft.fwd(vecTmpFreq, vecInputFFT, iNfft);
        matTapSpectrum.row(i) = vecTmpFreq;
    }

    return matTapSpectrum;
}

//=============================================================================================================

QVector<MatrixXcd> Spectral::computeTaperedSpectraMatrix(const MatrixXd &matData,
                                                         const MatrixXd &matTaper,
                                                         int iNfft,
                                                         bool bUseThreads)
{
    #ifdef EIGEN_FFTW_DEFAULT
        fftw_make_planner_thread_safe();
    #endif

    QVector<MatrixXcd> finalResult;

    if(!bUseThreads) {
        // Sequential
//        QElapsedTimer timer;
//        int iTime = 0;
//        int iTimeAll = 0;
//        timer.start();

        FFT<double> fft;
        fft.SetFlag(fft.HalfSpectrum);

        RowVectorXd vecInputFFT, rowData;
        RowVectorXcd vecTmpFreq;

        MatrixXcd matTapSpectrum(matTaper.rows(), int(floor(iNfft / 2.0)) + 1);
        int j;
        for (int i = 0; i < matData.rows(); ++i) {
            rowData = matData.row(i);

            //FFT for freq domain returning the half spectrum
            for (j = 0; j < matTaper.rows(); j++) {
                vecInputFFT = rowData.cwiseProduct(matTaper.row(j));
                fft.fwd(vecTmpFreq, vecInputFFT, iNfft);
                matTapSpectrum.row(j) = vecTmpFreq;
            }

            finalResult.append(matTapSpectrum);

//            iTime = timer.elapsed();
//            qDebug() << QThread::currentThreadId() << "Spectral::computeTaperedSpectraMatrix - Row-wise computation:" << iTime;
//            iTimeAll += iTime;
//            timer.restart();
        }

//        qDebug() << QThread::currentThreadId() << "Spectral::computeTaperedSpectraMatrix - Complete computation:" << iTimeAll;
    } else {
        // Parallel
        QList<TaperedSpectraInputData> lData;
        TaperedSpectraInputData dataTemp;
        dataTemp.matTaper = matTaper;
        dataTemp.iNfft = iNfft;

        for (int i = 0; i < matData.rows(); ++i) {
            dataTemp.vecData = matData.row(i);

            lData.append(dataTemp);
        }

        QFuture<QVector<MatrixXcd> > result = QtConcurrent::mappedReduced(lData,
                                                                          compute,
                                                                          reduce,
                                                                          QtConcurrent::OrderedReduce);
        result.waitForFinished();
        finalResult = result.result();
    }

    return finalResult;
}

//=============================================================================================================

MatrixXcd Spectral::compute(const TaperedSpectraInputData& inputData)
{
    //qDebug() << "Spectral::compute";
    return computeTaperedSpectraRow(inputData.vecData,
                                    inputData.matTaper,
                                    inputData.iNfft);
}

//=============================================================================================================

void Spectral::reduce(QVector<MatrixXcd>& finalData,
                      const MatrixXcd& resultData)
{
    //qDebug() << "Spectral::reduce";
    finalData.append(resultData);
}

//=============================================================================================================

Eigen::RowVectorXd Spectral::psdFromTaperedSpectra(const Eigen::MatrixXcd &matTapSpectrum,
                                                   const Eigen::VectorXd &vecTapWeights,
                                                   int iNfft,
                                                   double dSampFreq)
{
    //Check inputs
    if (matTapSpectrum.rows() != vecTapWeights.rows()) {
        return Eigen::RowVectorXd();
    }

    //Compute PSD (average over tapers if necessary)
    //Normalization via sFreq
    //multiply by 2 due to half spectrum
    double denom = vecTapWeights.cwiseAbs2().sum() * dSampFreq;
    Eigen::RowVectorXd vecPsd = 2.0 * (vecTapWeights.asDiagonal() * matTapSpectrum).cwiseAbs2().colwise().sum() / denom;

    vecPsd(0) /= 2.0;
    if (iNfft % 2 == 0){
        vecPsd.tail(1) /= 2.0;
    }

    return vecPsd;
}

//=============================================================================================================

QVector<Eigen::VectorXd> Spectral::psdFromTaperedSpectra(const QVector<Eigen::MatrixXcd> &matTapSpectrum,
                                                         const Eigen::VectorXd &vecTapWeights,
                                                         int iNfft,
                                                         double dSampFreq,
                                                         bool bUseMultithread)
{
#ifdef EIGEN_FFTW_DEFAULT
    fftw_make_planner_thread_safe();
#endif

    QVector<Eigen::VectorXd> finalResult;
    Eigen::RowVectorXd tmpResult;

    if(!bUseMultithread)
        for(int i=0; i< matTapSpectrum.size(); i++)
        {
            tmpResult = psdFromTaperedSpectra(matTapSpectrum.at(i),vecTapWeights, iNfft, dSampFreq);
            finalResult.append(tmpResult.transpose());
        }
    else
    {
        std::function<Eigen::VectorXd(const Eigen::MatrixXcd &)> loccompute = [&vecTapWeights, iNfft, dSampFreq](const Eigen::MatrixXcd & locTapSpectrum) {
            return psdFromTaperedSpectra(locTapSpectrum,vecTapWeights,iNfft,dSampFreq).transpose();
        };
        QFuture<Eigen::VectorXd> result = QtConcurrent::mapped(matTapSpectrum, loccompute);
        result.waitForFinished();
        finalResult = result.results().toVector();
    }

    return finalResult;
}


//=============================================================================================================

Eigen::RowVectorXcd Spectral::csdFromTaperedSpectra(const Eigen::MatrixXcd &vecTapSpectrumSeed,
                                                    const Eigen::MatrixXcd &vecTapSpectrumTarget,
                                                    const Eigen::VectorXd &vecTapWeightsSeed,
                                                    const Eigen::VectorXd &vecTapWeightsTarget,
                                                    int iNfft,
                                                    double dSampFreq)
{
//    QElapsedTimer timer;
//    int iTime = 0;
//    timer.start();

    //Check inputs
    if (vecTapSpectrumSeed.rows() != vecTapSpectrumTarget.rows()) {
        return Eigen::MatrixXcd();
    }
    if (vecTapSpectrumSeed.cols() != vecTapSpectrumTarget.cols()) {
        return Eigen::MatrixXcd();
    }
    if (vecTapSpectrumSeed.rows() != vecTapWeightsSeed.rows()) {
        return Eigen::MatrixXcd();
    }
    if (vecTapSpectrumTarget.rows() != vecTapWeightsTarget.rows()) {
        return Eigen::MatrixXcd();
    }

//    iTime = timer.elapsed();
//    qDebug() << QThread::currentThreadId() << "Spectral::csdFromTaperedSpectra timer - Prepare:" << iTime;
//    timer.restart();

    // Compute PSD (average over tapers if necessary)
    // Multiply by 2 due to half spectrum
    // Normalize via sFreq
    double denom = sqrt(vecTapWeightsSeed.cwiseAbs2().sum()) * sqrt(vecTapWeightsTarget.cwiseAbs2().sum()) * dSampFreq;
    Eigen::RowVectorXcd vecCsd = 2.0 * (vecTapWeightsSeed.asDiagonal() * vecTapSpectrumSeed).cwiseProduct((vecTapWeightsTarget.asDiagonal() * vecTapSpectrumTarget).conjugate()).colwise().sum() / denom;

//    iTime = timer.elapsed();
//    qDebug() << QThread::currentThreadId() << "Spectral::csdFromTaperedSpectra timer - compute PSD:" << iTime;
//    timer.restart();

    //multiply first and last element by 2 due to half spectrum
    vecCsd(0) /= 2.0;
    if (iNfft % 2 == 0){
        vecCsd.tail(1) /= 2.0;
    }

//    iTime = timer.elapsed();
//    qDebug() << QThread::currentThreadId() << "Spectral::csdFromTaperedSpectra timer - half spectrum:" << iTime;
//    timer.restart();

    return vecCsd;
}

//=============================================================================================================

VectorXd Spectral::calculateFFTFreqs(int iNfft, double dSampFreq)
{
    //Compute FFT frequencies
    RowVectorXd vecFFTFreqs;
    if (iNfft % 2 == 0){
        vecFFTFreqs = (dSampFreq / iNfft) * RowVectorXd::LinSpaced(iNfft / 2.0 + 1, 0.0, iNfft / 2.0);
    } else {
        vecFFTFreqs = (dSampFreq / iNfft) * RowVectorXd::LinSpaced((iNfft - 1) / 2.0 + 1, 0.0, (iNfft - 1) / 2.0);
    }
    return vecFFTFreqs;
}

//=============================================================================================================

QPair<MatrixXd, VectorXd> Spectral::generateTapers(int iSignalLength, const QString &sWindowType)
{
    QPair<MatrixXd, VectorXd> pairOut;
    if (sWindowType == "hanning") {
        pairOut.first = hanningWindow(iSignalLength);
        pairOut.second = VectorXd::Ones(1);
    } else if (sWindowType == "ones") {
        pairOut.first = MatrixXd::Ones(1, iSignalLength) / double(iSignalLength);
        pairOut.second = VectorXd::Ones(1);
    } else {
        pairOut.first = hanningWindow(iSignalLength);
        pairOut.second = VectorXd::Ones(1);
    }
    return pairOut;
}

//=============================================================================================================

QPair<Eigen::VectorXd, double> Spectral::calculateARWeightsMEMRow(const Eigen::RowVectorXd &vecData, int iOrder)
{
    int n = vecData.size();

    if (n < iOrder)
        qCritical() << "Error Occured in Spectral::calculateARweightssMEM: Number of data samples smaller than order of regression! ";

    Eigen::VectorXd coeffs = VectorXd::Zero(iOrder + 1);
    Eigen::VectorXd vecf = vecData;
    Eigen::VectorXd vecb = vecData;
    Eigen::VectorXd vecwkm = VectorXd::Zero(coeffs.size());

    double dE = vecf.dot(vecf);
    double dDen = 2.0*dE;

    dE /= static_cast<double>(n);

    double dNum = 0.0;
    double dq = 1.0;
    coeffs[0] = 1.0;

    for (int k=1; k <= iOrder; k++)
    {
        dNum = 0.0;

        for (int j=0; j < (n-k); j++) // mit dot product ersetzen
            dNum += vecf[j+1]*vecb[j];
        //dNum = vecf.segment(1,n-k).dot(vecb.segment(0,n-k));

        dDen = dq*dDen - vecf[0]*vecf[0] - vecb[n-k]*vecb[n-k]; // is this correct?

        coeffs[k] = -2.0*dNum/dDen;
        dq = 1.0 - coeffs[k]*coeffs[k];
        dE *= dq;

        for (int i=1; i<k; i++)
            coeffs[i] = vecwkm[i] + coeffs[k]*vecwkm[k-i]; // check +-
        //coeffs.segment(1,k-1) = vecwkm.segment(1,k-1) + coeffs[k]*vecwkm.segment(1,k-1).reverse();

        if (k < iOrder)
        {
            for (int i=1; i<=k; i++)
                vecwkm[i] = coeffs[i];
            //vecwkm.segment(1,k) = coeffs.segment(1,k);

            for (int j=0; j<n-k; j++) // evtl. mit eigen, braucht aber zwischenspeicher für vecf
            {
                vecf[j] = vecf[j+1] + vecwkm[k]*vecb[j]; // check +-
                vecb[j] = vecb[j] + vecwkm[k]*vecf[j+1]; // check +-
            }
        }

        // check indices!!!
    }

    //for (int i=1; i<=iOrder; i++)
    //    coeffs[i] *= -1; // check +-

    QPair<Eigen::VectorXd, double> out;

    out.first = coeffs;
    out.second = dE;

    return out;
}

//=============================================================================================================

QVector<QPair<Eigen::VectorXd, double>> Spectral::calculateARWeightsMEMMatrix(const Eigen::MatrixXd &matData, int iOrder, bool bUseMultithread)
{
#ifdef EIGEN_FFTW_DEFAULT
    fftw_make_planner_thread_safe();
#endif

    QVector<QPair<Eigen::VectorXd, double>> finalResult;
    QPair<Eigen::VectorXd, double> tmpResult;

    if(!bUseMultithread)
    {
        for (int i=0; i<matData.rows(); i++)
        {
            tmpResult = calculateARWeightsMEMRow(matData.row(i),iOrder);
            finalResult.append(tmpResult);
        }
    }
    else
    {
        QVector<ARWeightsMEMInputData> lData;
        ARWeightsMEMInputData tempData;
        tempData.iOrder = iOrder;

        for (int i=0; i<matData.rows(); i++)
        {
            tempData.vecData = matData.row(i);
            lData.append(tempData);
        }

        QFuture<QPair<Eigen::VectorXd, double>> result = QtConcurrent::mapped(lData, computeAR);
        result.waitForFinished();
        finalResult = result.results().toVector();
    }
    return finalResult;
}

//=============================================================================================================

QPair<Eigen::VectorXd, double> Spectral::computeAR(const ARWeightsMEMInputData& inputData)
{
    //qDebug() << "Spectral::compute";
    return calculateARWeightsMEMRow(inputData.vecData,
                                    inputData.iOrder);
}

//=============================================================================================================

/*void Spectral::reduceAR(QVector<QPair<Eigen::VectorXd, double>>& finalData,
                      const QPair<Eigen::VectorXd, double>& resultData)
{
    //qDebug() << "Spectral::reduce";
    finalData.append(resultData);
}*/

//=============================================================================================================


Eigen::MatrixXcd Spectral::generateARSpectraWeights(double dBottomFreq, double dTopFreq, int iBins, int iEvalPerBin, bool bCentered)
{
    if (dTopFreq < dBottomFreq)
        qCritical() << "Error Occured in Spectral::generateARSpectraWeights: Bottom frequency has to be smaller than top frequency! ";

    double twopi = 2*M_PI;
    double delta, tempfreq;

    if (bCentered)
    {
        delta = (dTopFreq - dBottomFreq) / (iBins * iEvalPerBin);
        tempfreq = twopi * (dBottomFreq + delta/2);
    }
    else
    {
        delta = (dTopFreq - dBottomFreq) / ((iBins - 1) * iEvalPerBin);
        tempfreq = twopi * (dBottomFreq - delta/2);
    }

    double deltapi = twopi * delta;

    Eigen::MatrixXcd mat_Theta = MatrixXcd::Zero(iBins,iEvalPerBin);

    for (int i=0; i<iBins; i++)
        for (int j=0; j<iEvalPerBin; j++)
        {
            mat_Theta(i,j) = std::polar<double>(1.0,tempfreq);
            tempfreq += deltapi;
        }

    return mat_Theta;
}

//=============================================================================================================

Eigen::VectorXd Spectral::psdFromARSpectra(const QPair<Eigen::VectorXd, double> &ARWeights,
                                           const Eigen::MatrixXcd &matSpectraWeights,  double dSampFreq)
{
    Eigen::VectorXd coeffs = ARWeights.first;
    double dNum = ARWeights.second;

    Eigen::MatrixXcd matSpectraWeightsWrk = matSpectraWeights;

    Eigen::MatrixXcd dDen = MatrixXcd::Ones(matSpectraWeights.rows(),matSpectraWeights.cols());

    for (int i=1; i<coeffs.size(); i++) // coeff[0] is one and already in dDen
    {
        dDen += coeffs(i)*matSpectraWeightsWrk;
        if (i<coeffs.size()-1)
            matSpectraWeightsWrk = matSpectraWeightsWrk.cwiseProduct(matSpectraWeights);
    }

    Eigen::MatrixXd SpectraBinned;

    SpectraBinned.array() = (2.0 * dNum) / (dDen.real().array().pow(2) + dDen.imag().array().pow(2)); // Multiply by 2 to account for symmetry of spectrum

    Eigen::VectorXd vecSpectraOut = SpectraBinned.rowwise().mean() / dSampFreq; // Take mean over evaluations in bin

    return vecSpectraOut;
}

//=============================================================================================================

QVector<Eigen::VectorXd> Spectral::psdFromARSpectra(const QVector<QPair<Eigen::VectorXd, double>> &ARWeights,
                                           const Eigen::MatrixXcd &matSpectraWeights,  double dSampFreq,
                                          bool bUseMultithread)
{
#ifdef EIGEN_FFTW_DEFAULT
    fftw_make_planner_thread_safe();
#endif

    QVector<Eigen::VectorXd> finalResult;
    Eigen::VectorXd tmpResult;

    if(!bUseMultithread)
        for(int i=0; i< ARWeights.size(); i++)
        {
            tmpResult = psdFromARSpectra(ARWeights.at(i),matSpectraWeights, dSampFreq);
            finalResult.append(tmpResult);
        }
    else
    {
        std::function<Eigen::VectorXd(const QPair<Eigen::VectorXd,double> &)> loccompute = [&matSpectraWeights, dSampFreq](const QPair<Eigen::VectorXd,double> & locARWeights) {
            return psdFromARSpectra(locARWeights,matSpectraWeights,dSampFreq);
        };
        QFuture<Eigen::VectorXd> result = QtConcurrent::mapped(ARWeights, loccompute);
        result.waitForFinished();
        finalResult = result.results().toVector();
    }

    return finalResult;
}

//=============================================================================================================

MatrixXd Spectral::hanningWindow(int iSignalLength)
{
    MatrixXd matHann = MatrixXd::Zero(1, iSignalLength);

    //Main step of building the hanning window
    for (int n = 0; n < iSignalLength; n++) {
        matHann(0, n) = 0.5 - 0.5 * cos(2.0 * M_PI * n / (iSignalLength - 1.0));
    }
    matHann.array() /= matHann.row(0).norm();

    return matHann;
}
