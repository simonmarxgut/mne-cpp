//=============================================================================================================
/**
 * @file     arsettingsview.h
 * @author   Johannes Vorwerk <johannes.vorwerk@umit-tirol.at>
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
 * @brief    Definition of the ARSettingsView Class.
 *
 */

#ifndef ARSETTINGSVIEW_H
#define ARSETTINGSVIEW_H

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

class QSpinBox;

namespace Ui {
class ARSettingsViewWidget;
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
class DISPSHARED_EXPORT ARSettingsView : public QWidget
{
    Q_OBJECT

public:
    typedef QSharedPointer<ARSettingsView> SPtr;              /**< Shared pointer type for SpectrumSettingsView. */
    typedef QSharedPointer<const ARSettingsView> ConstSPtr;   /**< Const shared pointer type for SpectrumSettingsView. */

    //=========================================================================================================
    /**
     * Constructs a SpectrumSettingsView which is a child of parent.
     *
     * @param [in] parent    parent of widget
     */
    ARSettingsView(const QString& sSettingsPath = "",
                   bool bChangeSamplingPoints = false,
                   QWidget *parent = 0,
                   Qt::WindowFlags f = Qt::Widget);

    //=========================================================================================================
    /**
     * Destroys the ScalingView.
     */
    ~ARSettingsView();

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
    void onUpdateSpinBoxAROrder();

    //=========================================================================================================
    /**
     * Slot called when scaling spin boxes change
     */
    void onUpdateSpinBoxNumEvaluationPoints();

    //=========================================================================================================
    /**
     * Slot called when scaling spin boxes change
     */
    void onUpdateSpinBoxSamplingPoints();

    //=========================================================================================================

    QString     m_sSettingsPath;
    bool        m_bChangeSamplingPoints;

    int         m_iOrder;
    int         m_iNumEvaluationPoints;
    int         m_iSamplingPoints;

    QSpinBox*   m_qSpinBoxSamplingPoints;

    Ui::ARSettingsViewWidget *ui;

signals:
    //=========================================================================================================
    /**
     * Emitted whenever the settings changed and are ready to be retreived.
     */
    void changeAROrder(int value);
    void changeARNumEvaluationPoints(int value);
};
} // NAMESPACE

#endif // ARSETTINGSVIEW_H
