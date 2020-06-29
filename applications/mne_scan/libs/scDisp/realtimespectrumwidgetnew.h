//=============================================================================================================
/**
 * @file     realtimespectrumwidgetnew.h
 * @author   Johannes Vorwerk <johannes.vorwerk@umit.at>
 * @version  dev
 * @date     April, 2020
 *
 * @section  LICENSE
 *
 * Copyright (C) 2020, Johannes Vorwerk. All rights reserved.
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
 * @brief     RealTimeSpectrumWidgetNew class declaration.
 *
 */

#ifndef SCDISPLIB_REALTIMESPECTRUMWIDGETNEW_H
#define SCDISPLIB_REALTIMESPECTRUMWIDGETNEW_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "scdisp_global.h"
#include "measurementwidget.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QSharedPointer>
#include <QPointer>
#include <QMap>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <Eigen/Core>

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

class QTime;
class QVBoxLayout;
class QLabel;

namespace FIFFLIB {
    class FiffInfo;
}

namespace SCMEASLIB {
    class RealTimeSpectrum;
}

namespace DISPLIB {
    class ModalitySelectionView;
    class ImageSc;
}

//=============================================================================================================
// DEFINE NAMESPACE SCDISPLIB
//=============================================================================================================

namespace SCDISPLIB {


//=============================================================================================================
// SCDISPLIB FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
/**
 * Description of what this class is intended to do (in detail).
 *
 * @brief Brief description of this class.
 */
class RealTimeSpectrumWidgetNew : public MeasurementWidget
{
    Q_OBJECT

public:
    //=========================================================================================================
    /**
     * Constructs a RealTimeSpectrumWidget which is a child of parent.
     *
     * @param [in] pNE           pointer to noise estimation measurement.
     * @param [in] pTime         pointer to application time.
     * @param [in] parent        pointer to parent widget; If parent is 0, the new NumericWidget becomes a window. If parent is another widget, NumericWidget becomes a child window inside parent. NumericWidget is deleted when its parent is deleted.
     */
    RealTimeSpectrumWidgetNew(QSharedPointer<SCMEASLIB::RealTimeSpectrum> pNE,
                            QSharedPointer<QTime> &pTime,
                            QWidget* parent = 0);

    typedef QSharedPointer<RealTimeSpectrumWidgetNew> SPtr;            /**< Shared pointer type for RealTimeSpectrumWidgetNew. */
    typedef QSharedPointer<const RealTimeSpectrumWidgetNew> ConstSPtr; /**< Const shared pointer type for RealTimeSpectrumWidgetNew. */

    //=========================================================================================================
    /**
     * Destroys the RealTimeSpectrumWidget.
     */
    ~RealTimeSpectrumWidgetNew();

    //=========================================================================================================
    /**
     * Is called when new data are available.
     *
     * @param [in] pMeasurement  pointer to measurement -> not used because its direct attached to the measurement.
     */
    virtual void update(SCMEASLIB::Measurement::SPtr pMeasurement);

    //=========================================================================================================
    /**
     * Is called when new data are available.
     */
    virtual void getData();

    //=========================================================================================================
    /**
     * Initialise the RealTimeSpectrumWidget.
     */
    virtual void init();

    //=========================================================================================================
    /**
     * Initialise the SettingsWidget.
     */
    void initSettingsWidget();

    //bool eventFilter(QObject *object, QEvent *event);

protected:
    QPointer<QAction>                       m_pActionSelectModality;                /**< Modality selection action */
    QPointer<QVBoxLayout>                   m_pFSLayout;                           /**< Widget layout */
    QPointer<QLabel>                        m_pLabelInit;                           /**< Initialization label */

    bool                                    m_bInitialized;                         /**< Is Initialized */

    QSharedPointer<FIFFLIB::FiffInfo>       m_pFiffInfo;                            /**< The Fiff Info. */

    QMap<QString, bool>                     m_modalityMap;                          /**< Map of different modalities. */

    QPointer<DISPLIB::ImageSc>              m_pImageSc;                             /**< The covariance colormap */

    QSharedPointer<SCMEASLIB::RealTimeSpectrum>                 m_pFS;                              /**< The frequency spectrum measurement. */

    Eigen::MatrixXd m_matSpectrumData;
    float m_fLowerFrqBound;         /**< Lower frequency bound */
    float m_fUpperFrqBound;         /**< Upper frequency bound */
    int m_iDisplayLength;
    int m_iDisplayIndex;
};

//=============================================================================================================
// INLINE DEFINITIONS
//=============================================================================================================


} // namespace SCDISPLIB

#endif // SCDISPLIB_REALTIMESPECTRUMWIDGETNEW_H
