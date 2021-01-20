//=============================================================================================================
/**
 * @file     normalizewidget.h
 * @author   Christoph Dinh <chdinh@nmr.mgh.harvard.edu>;
 *           Viktor Klueber <Viktor.Klueber@tu-ilmenau.de>;
 *           Lorenz Esch <lesch@mgh.harvard.edu>
 * @since    0.1.0
 * @date     January, 2016
 *
 * @section  LICENSE
 *
 * Copyright (C) 2016, Christoph Dinh, Viktor Klueber, Lorenz Esch. All rights reserved.
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
 * @brief    Contains the declaration of the NormalizeWidget class.
 *
 */

#ifndef NORMALIZEWIDGET_H
#define NORMALIZEWIDGET_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QWidget>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

namespace Ui{
    class NormalizeWidget;
}

//=============================================================================================================
// DEFINE NAMESPACE NORMALIZEPLUGIN
//=============================================================================================================

namespace NORMALIZEPLUGIN
{

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
/**
 * DECLARE CLASS NormalizeWidget
 *
 * @brief The NormalizeWidget class provides a normalize widget.
 */
class NormalizeWidget : public QWidget
{
    Q_OBJECT

public:
    typedef QSharedPointer<NormalizeWidget> SPtr;         /**< Shared pointer type for RereferenceWidget. */
    typedef QSharedPointer<NormalizeWidget> ConstSPtr;    /**< Const shared pointer type for RereferenceWidget. */

    //=========================================================================================================
    /**
     * Constructs a DummyToolbox.
     */
    explicit NormalizeWidget(const QString& sSettingsPath = "",
                             double dDataSampFreq = -1,
                             int iDataBlockLength = -1,
                             bool bEnabled = true,
                             int iMethod = 0,
                             int iModality = 0,
                             double dIntervallLength = 1.0,
                             QWidget *parent = 0);

    //=========================================================================================================
    /**
     * Destroys the DummyToolbox.
     */
    ~NormalizeWidget();

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
     * Slot called when radio button modality change
     */
    void onClickedButtonModality(int value);

    //=========================================================================================================
    /**
     * Slot called when radio button method change
     */
    void onClickedButtonMethod(int value);

    //=========================================================================================================
    /**
     * Slot called when checkbox enabled change
     */
    void onClickedCheckboxEnabled(bool value);

    //=========================================================================================================
    /**
     * Slot called when intervall length changed using data blocks spinbox
     */
    void onChangeIntervallLengthDataBlocks(int value);

    //=========================================================================================================
    /**
     * Slot called when intervall length changed using seconds spinbox
     */
    void onChangeIntervallLengthSeconds(double value);


    QString                     m_sSettingsPath;    /**< The settings path to store the GUI settings to. */
    double m_dDataSampFreq;
    int m_iDataBlockLength;
    bool m_bEnabled;
    int m_iModality;
    int m_iMethod;
    double m_dIntervallLength;

    Ui::NormalizeWidget*     ui;              /**< The UI class specified in the designer. */

signals:
    //=========================================================================================================
    /**
     * Emitted whenever the settings changed and are ready to be retreived.
     */
    void changeEnabled(bool value);
    void changeModality(int value);
    void changeMethod(int value);
    void changeIntervallLength(double value);
};
}   //namespace

#endif // NORMALIZEWIDGET_H
