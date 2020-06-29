//=============================================================================================================
/**
 * @file     realtimespectrumwidgetnew.cpp
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
 * @brief    RealTimeSpectrumWidgetNew class definition.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "realtimespectrumwidgetnew.h"

#include <scMeas/realtimespectrum.h>

#include <disp/viewers/modalityselectionview.h>
#include <disp/plots/imagesc.h>

#include <fiff/fiff_info.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QLabel>
#include <QFont>
#include <QDebug>
#include <QVBoxLayout>
#include <QSharedPointer>
#include <QAction>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <Eigen/Core>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace SCDISPLIB;
using namespace SCMEASLIB;
using namespace DISPLIB;
using namespace FIFFLIB;
using namespace Eigen;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

RealTimeSpectrumWidgetNew::RealTimeSpectrumWidgetNew(QSharedPointer<RealTimeSpectrum> pFS,
                                                 QSharedPointer<QTime> &pTime,
                                                 QWidget* parent)
: MeasurementWidget(parent)
, m_pFS(pFS)
, m_fLowerFrqBound(0)
, m_fUpperFrqBound(300)
, m_iDisplayLength(100)
, m_iDisplayIndex(0)
, m_bInitialized(false)
{
    Q_UNUSED(pTime)

    /*m_pActionFrequencySettings = new QAction(QIcon(":/images/frqResolution.png"), tr("Shows the frequency spectrum settings widget (F12)"),this);
    m_pActionFrequencySettings->setShortcut(tr("F12"));
    m_pActionFrequencySettings->setStatusTip(tr("Shows the frequency spectrum settings widget (F12)"));
    connect(m_pActionFrequencySettings.data(), &QAction::triggered,
            this, &RealTimeSpectrumWidget::showSpectrumSettingsView);
    addDisplayAction(m_pActionFrequencySettings);

    m_pActionFrequencySettings->setVisible(false);

    m_pSpectrumView = new SpectrumView(this, Qt::Window);*/

    //set vertical layout
    m_pFSLayout = new QVBoxLayout(this);

    m_pLabelInit= new QLabel;
    m_pLabelInit->setText("Acquiring Data");
    m_pLabelInit->setAlignment(Qt::AlignCenter);
    QFont font;font.setBold(true);font.setPointSize(20);
    m_pLabelInit->setFont(font);
    m_pFSLayout->addWidget(m_pLabelInit);

    m_pImageSc = new ImageSc;
    m_pFSLayout->addWidget(m_pImageSc);

    //set layouts
    this->setLayout(m_pFSLayout);

    getData();
}

//=============================================================================================================

RealTimeSpectrumWidgetNew::~RealTimeSpectrumWidgetNew()
{
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::update(SCMEASLIB::Measurement::SPtr)
{
    getData();
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::getData()
{
    if(!m_bInitialized)
    {
        if(m_pFS->isInit())
        {
            init();

            m_matSpectrumData = MatrixXd::Zero(m_pFS->getValue().rows(),m_iDisplayLength);

            if (!(m_iDisplayIndex<m_iDisplayLength))
                    m_iDisplayIndex = 0;

            m_matSpectrumData.col(m_iDisplayIndex) = m_pFS->getValue().col(0);
            m_iDisplayIndex++;

            m_pImageSc->updateData(m_matSpectrumData);
            //initSettingsWidget();
        }
    } else {
        if (!(m_iDisplayIndex<m_iDisplayLength))
                m_iDisplayIndex = 0;

        m_matSpectrumData.col(m_iDisplayIndex) = m_pFS->getValue().col(0);
        m_iDisplayIndex++;

        qDebug() << "m_pFS" << m_pFS->getValue().rows() << " " << m_pFS->getValue().cols();
        qDebug() << "m_matSpectrumData" << m_matSpectrumData.rows() << " " << m_matSpectrumData.cols();

        m_pImageSc->updateData(m_matSpectrumData);
    }
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::init()
{
    if(m_pFS->getValue().size() > 0)
    {
        m_pFiffInfo = m_pFS->getFiffInfo();

        m_pFSLayout->removeWidget(m_pLabelInit);
        m_pLabelInit->hide();

        m_pImageSc->setTitle(m_pFS->getName());

        //onNewModalitySelection(m_modalityMap);

        m_bInitialized = true;
    }
}

//=============================================================================================================
