//=============================================================================================================
/**
 * @file     bandpowersettingsview.h
 * @author   Johannes Vorwerk <johannes.vorwerk@umti-tirol.at>
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
 * @brief    Declaration of the BandPowerSettingsView Class.
 *
 */

#ifndef BANDPOWERSETTINGSVIEW_H
#define BANDPOWERSETTINGSVIEW_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../disp_global.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QWidget>
#include <QPointer>

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

class QSlider;
class QDoubleSpinBox;

namespace Ui {
class BandPowerSettingsViewWidget;
}

//=============================================================================================================
// DEFINE NAMESPACE DISPLIB
//=============================================================================================================

namespace DISPLIB
{

//=============================================================================================================
// DISPLIB FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
/**
 * DECLARE CLASS BandpowerSettingsView
 *
 * @brief The SpectrumSettingsView class provides settings for the spectrum estimation
 */
class DISPSHARED_EXPORT BandPowerSettingsView : public QWidget
{
    Q_OBJECT

public:
    typedef QSharedPointer<BandPowerSettingsView> SPtr;              /**< Shared pointer type for SpectrumSettingsView. */
    typedef QSharedPointer<const BandPowerSettingsView> ConstSPtr;   /**< Const shared pointer type for SpectrumSettingsView. */

    //=========================================================================================================
    /**
     * Constructs a BandPowerSettingsView which is a child of parent.
     *
     * @param [in] parent    parent of widget
     */
    BandPowerSettingsView(const QString& sSettingsPath = "",
                          double dSampleFreq = 512,
                          double dMin = 10,
                          double dMax = 20,
                          const QString& sSpectrumMethod = QString::fromStdString("AR"),
                          int iIntervallLength = 1,
                          int iChannels = 1,
                          int iBins = 1,
                          int iDetrend = 0,
                          bool bIsRunning = false,
                          QWidget *parent = 0,
                          Qt::WindowFlags f = Qt::Widget);

    //=========================================================================================================
    /**
     * Destroys the ScalingView.
     */
    ~BandPowerSettingsView();

    //=========================================================================================================
    /**
     * Get min value
     */
    double getMin();

    //=========================================================================================================
    /**
     * Get max slider value
     */
    double getMax();

    //=========================================================================================================
    /**
     * Get max slider value
     */
    QString getMethod();

    //=========================================================================================================
    /**
     * Emits all settings once.
     *
     */
    void emitSignals();

protected:
    //=========================================================================================================
    /**
     * Saves all important settings of this view via QSettings.
     *
     * @param[in] settingsPath        the path to store the settings to.
     */
    void saveSettings(const QString& settingsPath);

    //=========================================================================================================
    /**
     * Loads and inits all important settings of this view via QSettings.
     *
     * @param[in] settingsPath        the path to load the settings from.
     */
    void loadSettings(const QString& settingsPath);

    //=========================================================================================================
    /**
     * Slot called when scaling spin boxes change
     */
    void onUpdateSpinBoxMinScaling();

    //=========================================================================================================
    /**
     * Slot called when scaling spin boxes change
     */
    void onUpdateSpinBoxMaxScaling();

    //=========================================================================================================
    /**
     * Slot called when slider scaling change
     */
    void onUpdateSliderMinScaling();

    //=========================================================================================================
    /**
     * Slot called when slider scaling change
     */
    void onUpdateSliderMaxScaling();

    //=========================================================================================================
    /**
     * Slot called when radio button method change
     */
    void onClickedButtonMethod(int value);

    //=========================================================================================================
    /**
     * Slot called when scaling spin box changes
     */
    void onUpdateSpinBoxIntervallLength();

    //=========================================================================================================
    /**
     * Slot called when number of channels spin box changes
     */
    void onUpdateSpinBoxChannels(int value);

    //=========================================================================================================
    /**
     * Slot called when number of bins spin box changes
     */
    void onUpdateSpinBoxBins(int value);

    //=========================================================================================================
    /**
     * Slot called when radio button method change
     */
    void onClickedButtonDetrend(int value);

    QString     m_sSettingsPath;
    double      m_dSampFreq;
    double      m_dMin;
    double      m_dMax;
    QString     m_sSpectrumMethod;
    int         m_iIntervallLength;
    int         m_iChannels;
    int         m_iBins;
    int         m_iDetrend;
    bool        m_bIsRunning;

    Ui::BandPowerSettingsViewWidget *ui;

signals:
    //=========================================================================================================
    /**
     * Emitted whenever the settings changed and are ready to be retreived.
     */
    void changeMinMax(double min, double max);
    void changeMethod(const QString& sSpectrumMethod);
    void changeIntervallLength(int value);
    void changeDetrend(int value);
    void changeChannels(int value);
    void changeBins(int value);

};
} // NAMESPACE

#endif // BANDPOWERSETTINGSVIEW_H
