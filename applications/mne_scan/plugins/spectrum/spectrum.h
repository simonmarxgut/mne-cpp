//=============================================================================================================
/**
 * @file     spectrum.h
 * @author   Christoph Dinh <chdinh@nmr.mgh.harvard.edu>;
 *           Lorenz Esch <lesch@mgh.harvard.edu>
 * @version  dev
 * @date     April, 2020
 *
 * @section  LICENSE
 *
 * Copyright (C) 2020, Christoph Dinh, Lorenz Esch, Johannes Vorwerk. All rights reserved.
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
 * @brief    Contains the declaration of the Spectrum class.
 *
 */

#ifndef SPECTRUM_H
#define SPECTRUM_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "spectrum_global.h"

#include <utils/generics/circularbuffer.h>
#include <scMeas/realtimemultisamplearray.h>
#include <scMeas/realtimespectrum.h>

#include <utils/spectral.h>

#include <scShared/Interfaces/IAlgorithm.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QPointer>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <Eigen/SparseCore>

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

namespace SCMEASLIB {
    class RealTimeMultiSampleArray;
}


namespace DISPLIB {
    class ChannelSelectionView;
    class ChannelInfoModel;
}

//=============================================================================================================
// DEFINE NAMESPACE SPECTRUMPLUGIN
//=============================================================================================================

namespace SPECTRUMPLUGIN
{

//=============================================================================================================
// SPECTRUMPLUGIN FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
/**
 * DECLARE CLASS Spectrum
 *
 * @brief The Spectrum class provides a dummy algorithm structure.
 */
class SPECTRUMSHARED_EXPORT Spectrum : public SCSHAREDLIB::IAlgorithm
{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "scsharedlib/1.0" FILE "spectrum.json") //New Qt5 Plugin system replaces Q_EXPORT_PLUGIN2 macro
    // Use the Q_INTERFACES() macro to tell Qt's meta-object system about the interfaces
    Q_INTERFACES(SCSHAREDLIB::IAlgorithm)

public:
    //=========================================================================================================
    /**
     * Constructs a Spectrum.
     */
    Spectrum();

    //=========================================================================================================
    /**
     * Destroys the Spectrum.
     */
    ~Spectrum();

    //=========================================================================================================
    /**
     * IAlgorithm functions
     */
    virtual QSharedPointer<SCSHAREDLIB::IPlugin> clone() const;
    virtual void init();
    virtual void unload();
    virtual bool start();
    virtual bool stop();
    virtual SCSHAREDLIB::IPlugin::PluginType getType() const;
    virtual QString getName() const;
    virtual QWidget* setupWidget();

    //=========================================================================================================
    /**
     * Udates the pugin with new (incoming) data.
     *
     * @param[in] pMeasurement    The incoming data in form of a generalized Measurement.
     */
    void update(SCMEASLIB::Measurement::SPtr pMeasurement);

    //=========================================================================================================
    /**
     * Inits widgets which are used to control this plugin, then emits them in form of a QList.
     */
    virtual void initPluginControlWidgets();

    //=========================================================================================================
    /**
             * Changes the bottom and top frequency of the power spectrum.
             *
             * @param[in] frequency  Desired bottom frequency.
             */
    void changeMinMaxFrequency(double minFrequency, double maxFrequency);

    //=========================================================================================================
    /**
             * Changes the method to calculate the power spectrum.
             *
             * @param[in] method  Method to calculate the power spectrum.
             */
    void changeSpectrumMethod(QString method);

    //=========================================================================================================
    /**
     * Changes the length of the evaluated intervall (in multiples of block length).
     *
     * @param[in] lengthfactor  Factor for intervall length.
     */
    void changeIntervallLengthFactor(quint32 lengthfactor);

    //=========================================================================================================
    /**
     * Changes the resolution of the power spectrum.
     *
     * @param[in] resolutionfactor  Resolution of the power spectrum gets reduced by factor.
     */
    void changeResolutionFactor(quint32 resolutionfactor);

    //=========================================================================================================
    /**
             * Changes the detrend method for the data.
             *
             * @param[in] method  Detrend method (0: none, 1: remove mean, 2: remove linear trend).
             */
    void changeDetrendMethod(quint32 method);

    //=========================================================================================================
    /**
             * Changes the number of channels for which bandpower is analyzed.
             *
             * @param[in] channels  Number of channels.
             */
    void changeBandPowerChannels(quint32 channels);

    //=========================================================================================================
    /**
             * Changes the number of bins in which bandpower is calculated for each channel.
             *
             * @param[in] bins  Number of bins.
             */
    void changeSpectrumBins(quint32 bins);

    //=========================================================================================================
    /**
             * Changes the order of the AR for spectrum estimation.
             *
             * @param[in] order  Order of the AR method.
             */
    void changeAROrder(quint32 order);

    //=========================================================================================================
    /**
             * Changes the number of evaluatin points for the AR for spectrum estimation.
             *
             * @param[in] evaluationpoints  Number of evaluation points.
             */
    void changeARNumEvaluationPoints(quint32 evaluationpoints);

    //=========================================================================================================
    /**
     * Only evaluates the channels defined in the QStringList selectedChannels
     *
     * @param [in] selectedChannels list of all channel names which are currently selected in the selection manager.
     */
    void evaluateSelectedChannelsOnly(const QStringList& selectedChannels);

    //=========================================================================================================
    /**
     * Calculate Frequency bins of Spectrum
     *
     */
    void calcFrequencies();

protected:
    //=========================================================================================================
    /**
     * IAlgorithm function
     */
    virtual void run();

    //void showYourWidget();

private:

    //=========================================================================================================
    /**
     * Shows sensor selection widget
     */
    void showSensorSelectionWidget();

    QMutex                          m_qMutex;                                    /**< The threads mutex.*/

    bool        m_bProcessData;                     /**< If data should be received for processing */

    quint32 m_iIntervallLengthFactor;
    quint32 m_iResolutionFactorFFT;
    quint32 m_iOrderAR;
    quint32 m_iEvalsAR;
    QString m_sSpectrumMethod;
    quint32 m_iDetrendMethod;
    qint32 m_iSpectrumBins;

    double m_dDataSampFreq;

    double m_dFreqMin;
    double m_dFreqMax;
    Eigen::VectorXd m_pFrequencies;

    FIFFLIB::FiffInfo::SPtr                         m_pFiffInfo_orig;       /**< Input Fiff measurement info.*/
    FIFFLIB::FiffInfo::SPtr                         m_pFiffInfo;            /**< Fiff measurement info.*/
    int                                             m_iNTimeSteps;
    int                                             m_iNChannels;
    //QSharedPointer<DummyYourWidget>                 m_pYourWidget;          /**< flag whether thread is running.*/
    QAction*                                        m_pActionShowYourWidget;/**< flag whether thread is running.*/
    QPointer<QAction>                               m_pActionSelectSensors;         /**< show roi select widget */

    QSharedPointer<DISPLIB::ChannelInfoModel>               m_pChannelInfoModel;            /**< channel info model. */
    QSharedPointer<DISPLIB::ChannelSelectionView>           m_pChannelSelectionView;        /**< ChannelSelectionView. */
    QList<int>                                    m_pSelectedChannels;

    IOBUFFER::CircularBuffer<Eigen::MatrixXd>::SPtr    m_pSpectrumBuffer;         /**< Holds incoming data.*/

    SCSHAREDLIB::PluginInputData<SCMEASLIB::RealTimeMultiSampleArray>::SPtr      m_pSpectrumInput;      /**< The RealTimeMultiSampleArray of the DummyToolbox input.*/
    SCSHAREDLIB::PluginOutputData<SCMEASLIB::RealTimeSpectrum>::SPtr     m_pSpectrumOutput;     /**< The RealTimeMultiSampleArray of the DummyToolbox output.*/
    // SCSHAREDLIB::PluginOutputData<SCMEASLIB::RealTimeMultiSampleArray>::SPtr     m_pSpectrumOutput;     /**< The RealTimeMultiSampleArray of the DummyToolbox output.*/


signals:
    //=========================================================================================================
    /**
     * Emitted when fiffInfo is available
     */
    void fiffInfoAvailable();
};
} // NAMESPACE

#endif // SPECTRUM_H
