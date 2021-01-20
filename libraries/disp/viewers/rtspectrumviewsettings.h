//=============================================================================================================
/**
 * @file     rtspectrumviewsettings.h
 * @author   Lorenz Esch <lesch@mgh.harvard.edu>
 * @since    0.1.0
 * @date     July, 2018
 *
 * @section  LICENSE
 *
 * Copyright (C) 2018, Lorenz Esch. All rights reserved.
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
 * @brief    Declaration of the RtSpectrumViewSettings Class.
 *
 */

#ifndef RTSPECTRUMVIEWSETTINGS_H
#define RTSPECTRUMVIEWSETTINGS_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../disp_global.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QWidget>
#include <QStringList>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

namespace Ui {
    class RtSpectrumViewSettingsWidget;
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
 * DECLARE CLASS RtSpectrumViewSettings
 *
 * @brief The RtSpectrumViewSettings class provides a view to select different channel data view dependent settings
 */
class DISPSHARED_EXPORT RtSpectrumViewSettings : public QWidget
{
    Q_OBJECT

public:    
    typedef QSharedPointer<RtSpectrumViewSettings> SPtr;              /**< Shared pointer type for FiffRawViewSettings. */
    typedef QSharedPointer<const RtSpectrumViewSettings> ConstSPtr;   /**< Const shared pointer type for FiffRawViewSettings. */

    //=========================================================================================================
    /**
     * Constructs a FiffRawViewSettings which is a child of parent.
     *
     * @param [in] parent        parent of widget
     */
    RtSpectrumViewSettings(const QString& sSettingsPath = "", double dWindowLength = 10.0, double dWindowStepsize = 1.0,
                            QWidget *parent = 0,
                            Qt::WindowFlags f = Qt::Widget);

    //=========================================================================================================
    /**
     * Destroys the FiffRawViewSettings.
     */
    ~RtSpectrumViewSettings();

    //=========================================================================================================
    /**
     * Sets the values of the windowSize spin box
     *
     * @param [in] windowSize    new window size value
     */
    void setWindowSize(int windowSize);

    //=========================================================================================================
    /**
     * Set current distance time spacer combo box.
     *
     * @param [in] value     the new value of the combo box
     */
    void setColormapMax(double value);

    //=========================================================================================================
    /**
     * Set current distance time spacer combo box.
     *
     * @param [in] value     the new value of the combo box
     */
    void setColormapMin(double value);

    //=========================================================================================================
    /**
     * Returns the current window size.
     *
     * @return The current window size.
     */
    int getWindowSize();

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
     * Slot called when time window size changes
     */
    void onTimeWindowChanged(double value);

    //=========================================================================================================
    /**
     * Slot called when zoome changes
     */
    void onZoomChanged(double value);

    //=========================================================================================================
    /**
     * Call this slot whenever you want to make a screenshot of the butterfly or layout view.
     */
    void onMakeScreenshot();

    Ui::RtSpectrumViewSettingsWidget* ui;

    QString     m_sSettingsPath;                /**< The settings path to store the GUI settings to. */

    //=========================================================================================================
    /**
     * Slot called when zoome changes
     */
    void onFixColormapChanged(bool value);

    //=========================================================================================================
    /**
     * Slot called when zoome changes
     */
    void onColormapMaxChanged(double value);

    //=========================================================================================================
    /**
     * Slot called when zoome changes
     */
    void onColormapMinChanged(double value);

signals:
    //=========================================================================================================
    /**
     * Emit this signal whenever the user changes the window size.
     */
    void timeWindowChanged(double value);

    //=========================================================================================================
    /**
     * Emit this signal whenever the user wants to make a screenshot.
     *
     * @param[out] imageType     The current image type: png, svg.
     */
    void makeScreenshot(const QString& imageType);

    //=========================================================================================================
    /**
     * Emit this signal whenever the user changes the row height (zoom) of the channels.
     */
    void fixColormapChanged(bool value);

    //=========================================================================================================
    /**
     * Emit this signal whenever the user changes the row height (zoom) of the channels.
     */
    void colormapMaxChanged(double value);

    //=========================================================================================================
    /**
     * Emit this signal whenever the user changes the row height (zoom) of the channels.
     */
    void colormapMinChanged(double value);
};
} // NAMESPACE

#endif // RTSPECTRUMVIEWSETTINGS_H
