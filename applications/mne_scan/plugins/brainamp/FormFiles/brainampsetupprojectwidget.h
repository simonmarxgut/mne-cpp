//=============================================================================================================
/**
 * @file     brainampsetupprojectwidget.h
 * @author   Lorenz Esch <lesch@mgh.harvard.edu>;
 *           Viktor Klueber <Viktor.Klueber@tu-ilmenau.de>
 * @since    0.1.0
 * @date     October, 2016
 *
 * @section  LICENSE
 *
 * Copyright (C) 2016, Lorenz Esch, Viktor Klueber. All rights reserved.
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
 * @brief    Contains the declaration of the BrainAMPSetupProjectWidget class.
 *
 */

#ifndef BRAINAMPSETUPPROJECTWIDGET_H
#define BRAINAMPSETUPPROJECTWIDGET_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QWidget>

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

namespace Ui {
    class BrainAMPSetupProjectWidget;
}

//=============================================================================================================
// DEFINE NAMESPACE BRAINAMPPLUGIN
//=============================================================================================================

namespace BRAINAMPPLUGIN
{

//=============================================================================================================
// BRAINAMPPLUGIN FORWARD DECLARATIONS
//=============================================================================================================

class BrainAMP;

//=============================================================================================================
/**
 * DECLARE CLASS BrainAMPSetupProjectWidget
 *
 * @brief The BrainAMPSetupProjectWidget class provides the BrainAMPSetupProjectWidget configuration window.
 */
class BrainAMPSetupProjectWidget : public QWidget
{
    Q_OBJECT

public:
    //=========================================================================================================
    /**
     * Constructs a BrainAMPSetupProjectWidget which is a child of parent.
     *
     * @param[in] parent pointer to parent widget; If parent is 0, the new BrainAMPSetupProjectWidget becomes a window. If parent is another widget, BrainAMPSetupWidget becomes a child window inside parent. BrainAMPSetupWidget is deleted when its parent is deleted.
     * @param[in] pBrainAMP a pointer to the corresponding ECGSimulator.
     */
    explicit BrainAMPSetupProjectWidget(BrainAMP* pBrainAMP,
                                        QWidget *parent = 0);

    //=========================================================================================================
    /**
     * Destructs a BrainAMPSetupProjectWidget which is a child of parent.
     *
     */
    ~BrainAMPSetupProjectWidget();

    //=========================================================================================================
    /**
     * Inits the GUI
     *
     */
    void initGui();

private:
    BrainAMP*                           m_pBrainAMP;        /**< a pointer to corresponding BrainAMP.*/
    Ui::BrainAMPSetupProjectWidget*     ui;                 /**< the user interface for the BrainAMPSetupWidget.*/

    //=========================================================================================================
    /**
     * Sets the project dir
     */
    void changeCardinalMode(const QString &text);

    //=========================================================================================================
    /**
     * Sets the project dir
     */
    void onCardinalComboBoxChanged();

    //=========================================================================================================
    /**
     * Sets the project dir
     *
     * @param[in] sPath  The file path of the 3D layout file.
     */
    void updateCardinalComboBoxes(const QString& sPath);

    //=========================================================================================================
    /**
     * Sets the dir where the eeg cap file is located
     */
    void changeCap();

    //=========================================================================================================
    /**
     * Sets the dir where the cardinal file is located
     */
    void changeCardinalFile();

    //=========================================================================================================
    /**
     * Changes the EEG cap and file path variables in the BrainAMP class
     */
    void changeQLineEdits();

signals:    
    //=========================================================================================================
    /**
     * Emit this signal whenever the cardinal points changed.
     *
     * @param[in] sLPA       The channel name to take as the LPA.
     * @param[in] dLPA       The amount (in m) to translate the LPA channel position on the z axis.
     * @param[in] sRPA       The channel name to take as the RPA.
     * @param[in] dRPA       The amount (in m) to translate the RPA channel position on the z axis.
     * @param[in] sNasion    The channel name to take as the Nasion.
     * @param[in] dNasion    The amount (in m) to translate the Nasion channel position on the z axis.
     */
    void cardinalPointsChanged(const QString& sLPA,
                               double dLPA,
                               const QString& sRPA,
                               double dRPA,
                               const QString& sNasion,
                               double dNasion);
};
} // NAMESPACE

#endif // BRAINAMPSETUPPROJECTWIDGET_H
