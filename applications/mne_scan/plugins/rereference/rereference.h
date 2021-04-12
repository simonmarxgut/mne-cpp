//=============================================================================================================
/**
 * @file     rereference.h
 * @author   Christoph Dinh <chdinh@nmr.mgh.harvard.edu>;
 *           Lorenz Esch <lesch@mgh.harvard.edu>;
 *           Viktor Klueber <Viktor.Klueber@tu-ilmenau.de>
 * @since    0.1.0
 * @date     February, 2013
 *
 * @section  LICENSE
 *
 * Copyright (C) 2013, Christoph Dinh, Lorenz Esch, Viktor Klueber. All rights reserved.
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
 * @brief    Contains the declaration of the Rereference class.
 *
 */

#ifndef REREFERENCE_H
#define REREFERENCE_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "rereference_global.h"

#include <scShared/Plugins/abstractalgorithm.h>
#include <utils/generics/circularbuffer.h>
#include <scMeas/realtimemultisamplearray.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QtWidgets>
#include <QtCore/QtPlugin>
#include <QDebug>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

namespace DISPLIB {
    class ChannelSelectionView;
    class ChannelInfoModel;
}

namespace SCMEASLIB {
    class RealTimeMultiSampleArray;
}

//=============================================================================================================
// DEFINE NAMESPACE REREFERENCEPLUGIN
//=============================================================================================================

namespace REREFERENCEPLUGIN
{

//=============================================================================================================
// REREFERENCEPLUGIN FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
/**
 * DECLARE CLASS Rereference
 *
 * @brief The Rereference class provides a dummy algorithm structure.
 */

class REREFERENCESHARED_EXPORT Rereference : public SCSHAREDLIB::AbstractAlgorithm
{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "scsharedlib/1.0" FILE "rereference.json") //New Qt5 Plugin system replaces Q_EXPORT_PLUGIN2 macro
    // Use the Q_INTERFACES() macro to tell Qt's meta-object system about the interfaces
    Q_INTERFACES(SCSHAREDLIB::AbstractAlgorithm)

public:
    //=========================================================================================================
    /**
     * Constructs a rereference.
     */
    Rereference();

    //=========================================================================================================
    /**
     * Destroys the rereference.
     */
    ~Rereference();

    //=========================================================================================================
    /**
     * IAlgorithm functions
     */
    virtual QSharedPointer<SCSHAREDLIB::AbstractPlugin> clone() const;
    virtual void init();
    virtual void unload();
    virtual bool start();
    virtual bool stop();
    virtual SCSHAREDLIB::AbstractPlugin::PluginType getType() const;
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
     * Only evaluates the channels defined in the QStringList selectedChannels
     *
     * @param [in] selectedChannels list of all channel names which are currently selected in the selection manager.
     */
    void evaluateSelectedChannelsOnly(const QStringList& selectedChannels);

    //=========================================================================================================
    /**
     * Update the matrix for rereferencing
     */
    void updateRereferenceMatrix();

    //=========================================================================================================
    /**
     * Change modality for which rereferencing is performed
     *
     * @param [in] modality (0 = EEG, 1 = MEG, 2 = EEG + MEG)
     */
    void changeModality(int modality);

    //=========================================================================================================
    /**
     * Change rereferencing method
     *
     * @param [in] rerefrencing method (0 = CAR, 1 = Selection, 2 = File)
     */
    void changeMethod(int method);

    //=========================================================================================================
    /**
     * Change enabled status
     *
     * @param [in] set rereferencing enabled
     */
    void changeEnabled(bool enabled);

    //=========================================================================================================
    /**
     * Change matrix filename
     *
     * @param [in] rereferencing matrix filename
     */
    void changeMatrixFilename(QString sMatrixFilename);

protected:
    //=========================================================================================================
    /**
     * Inits widgets which are used to control this plugin, then emits them in form of a QList.
     */
    virtual void initPluginControlWidgets();

    //=========================================================================================================
    /**
     * IAlgorithm function
     */
    virtual void run();

private:

    //=========================================================================================================
    /**
     * Shows sensor selection widget
     */
    void showSensorSelectionWidget();

    QMutex                          m_qMutex;                                    /**< The threads mutex.*/

    FIFFLIB::FiffInfo::SPtr                         m_pFiffInfo;                /**< Fiff measurement info.*/
    int                                             m_iNChannels;               /** number of incoming channels */
    bool                                            m_bRereferenceEnabled;
    int                                             m_iRereferenceMethod;
    int                                             m_iRereferenceModality;
    Eigen::RowVectorXi                              m_pEEGindex;                /** indices of EEG channels */
    Eigen::RowVectorXi                              m_pMEGindex;                /** indices of MEG channels */
    Eigen::RowVectorXi                              m_pMEEGindex;                /** indices of MEEG channels */
    Eigen::MatrixXd                                 m_pRereferenceMatrix;             /** matrix for rereferencing */
    QString                                         m_sMatrixFilename;

    QSharedPointer<Rereference>                     m_pYourWidget;              /**< The widget used to control this plugin by the user.*/
    QPointer<QAction>                               m_pActionSelectSensors;         /**< show roi select widget */

    QSharedPointer<DISPLIB::ChannelInfoModel>       m_pChannelInfoModel;            /**< channel info model. */
    QSharedPointer<DISPLIB::ChannelSelectionView>   m_pChannelSelectionView;        /**< ChannelSelectionView. */
    QList<int>                                      m_pSelectedChannels;

    UTILSLIB::CircularBuffer<Eigen::MatrixXd>::SPtr m_pRereferenceBuffer;          /**< Holds incoming data.*/

    SCSHAREDLIB::PluginInputData<SCMEASLIB::RealTimeMultiSampleArray>::SPtr      m_pRereferenceInput;      /**< The incoming data.*/
    SCSHAREDLIB::PluginOutputData<SCMEASLIB::RealTimeMultiSampleArray>::SPtr     m_pRereferenceOutput;     /**< The outgoing data.*/

signals:
    //=========================================================================================================
    /**
     * Emitted when fiffInfo is available
     */
    void fiffInfoAvailable();
};
} // NAMESPACE

#endif // REREFERENCE_H
