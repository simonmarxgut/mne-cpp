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

#include <disp/viewers/rtspectrumviewsettings.h>
#include <disp/plots/imagesc.h>

#include <fiff/fiff_info.h>

#include <disp/plots/helpers/colormap.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QLabel>
#include <QFont>
#include <QDebug>
#include <QGridLayout>
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

RealTimeSpectrumWidgetNew::RealTimeSpectrumWidgetNew(QSharedPointer<QTime> &pTime,
                                                 QWidget* parent)
: MeasurementWidget(parent)
, m_sColorMap("Jet")
, m_dMinValue(0)
, m_dMaxValue(256)
, m_dMinLimit(0)
, m_dMaxLimit(256)
, m_dXSampleRate(2)
, m_dYSampleRate(2)
, m_dMaxFreq(30)
, m_dMinFreq(4)
, m_iXLimit(100)
, m_iDisplayIndex(0)
, m_bColormapLimits(false)
, m_bInitialized(false)
{
    Q_UNUSED(pTime)

    //set layouts
    this->setMinimumSize(300,50);
}

//=============================================================================================================

RealTimeSpectrumWidgetNew::~RealTimeSpectrumWidgetNew()
{
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::update(SCMEASLIB::Measurement::SPtr pMeasurement)
{
    if(!m_pFS) {
        m_pFS = qSharedPointerDynamicCast<RealTimeSpectrum>(pMeasurement);
    }

    if(m_pFS->isInit() && !m_pFiffInfo)
    {
        m_pFiffInfo = m_pFS->getFiffInfo();

        m_dXSampleRate = m_pFiffInfo->sfreq;

        // set min/max frequency
        m_dMinFreq = m_pFS->getValue().col(1).minCoeff();
        m_dMaxFreq = m_pFS->getValue().col(1).maxCoeff();

        // set data matrix
        m_matSpectrumData = MatrixXd::Zero(m_pFS->getValue().rows(),m_iXLimit);

        qDebug() << "RealTimeSpectrumWidgetNew::update m_matSpectrumData" << m_pFS->getValue().rows() << m_iXLimit;

        m_dYSampleRate = (m_dMaxFreq - m_dMinFreq)/m_matSpectrumData.rows();

        if(!m_bInitialized) {
            //m_pLayout->removeWidget(m_pLabelInit);
            //m_pLabelInit->hide();
            initDisplayWidget();
            initSettingsWidget();
        }
    } else if(!(m_pFS->getValue().cols() == 0))
    {
        getData();
        updateDataPixmap();
    }
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::getData()
{
    // check if frequencies changed
    if ((m_dMinFreq != m_pFS->getValue().col(1).minCoeff())||(m_dMaxFreq != m_pFS->getValue().col(1).maxCoeff()))
    {
        m_dMinFreq = m_pFS->getValue().col(1).minCoeff();
        m_dMaxFreq = m_pFS->getValue().col(1).maxCoeff();
        updateCoeffsPixmap();

        m_matSpectrumData = MatrixXd::Zero(m_pFS->getValue().rows(),m_iXLimit);
        m_iDisplayIndex = 0;
    }

    if (!(m_iDisplayIndex<m_iXLimit))
        m_iDisplayIndex = 0;

    if (m_pFS->getValue().rows() == m_matSpectrumData.rows()) {
        m_matSpectrumData.col(m_iDisplayIndex) = m_pFS->getValue().col(0);
        updateData();
        m_iDisplayIndex++;
    } else {
        m_matSpectrumData = MatrixXd::Zero(m_pFS->getValue().rows(),m_iXLimit);
        m_matSpectrumData.col(0) = m_pFS->getValue().col(0);
        updateData();
        m_iDisplayIndex = 0;
    }
    //updateCoeffsPixmap();
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::init()
{
    //Colormap
    pColorMapper = ColorMap::valueToColor;

}

//=============================================================================================================
void RealTimeSpectrumWidgetNew::initDisplayWidget()
{
    //setup image
    QImage * image_to_tf_plot = new QImage(m_matSpectrumData.cols(), m_matSpectrumData.rows(), QImage::Format_RGB32);

    //setup pixelcolors in image
    QColor color;
    for ( qint32 y = 0; y < m_matSpectrumData.rows(); y++ ) {
        for ( qint32 x = 0; x < m_matSpectrumData.cols(); x++ ) {
            color.setRgb(0,0,0);
            image_to_tf_plot->setPixel(x, m_matSpectrumData.rows() - 1 -  y,  color.rgb());
        }
    }

    *image_to_tf_plot = image_to_tf_plot->scaled(m_matSpectrumData.cols(), m_matSpectrumData.cols()/2, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
    *image_to_tf_plot = image_to_tf_plot->scaledToWidth(this->width()-200, Qt::SmoothTransformation);
    //image to pixmap
    m_pFSPixmap = new QGraphicsPixmapItem(QPixmap::fromImage(*image_to_tf_plot));
    //m_pFSPixmap->setScale(100);
    m_pFSScene = new QGraphicsScene();
    m_pFSScene->addItem(m_pFSPixmap);

    m_pLayout = new QGridLayout();
    m_pLayout->setContentsMargins(0,0,0,0);
    m_pView = new QGraphicsView();
    m_pView->setObjectName("tf_view");
    m_pView->setScene(m_pFSScene);

    m_pFSScene->setBackgroundBrush(Qt::white);

    QImage * coeffs_image = new QImage(5, m_matSpectrumData.rows(), QImage::Format_RGB32);
    //qreal norm = m_matSpectrumData.maxCoeff();
    //setup pixelcolors in image
    for(qint32 it = 0; it < m_matSpectrumData.rows(); it++) {
        for ( qint32 x = 0; x < 5; x++ ) {
            // color.setRgb(pColorMapper(it*norm/m_matSpectrumData.rows(),m_sColorMap));
            color.setRgb(pColorMapper(1.0 - double(it)/double(m_matSpectrumData.rows()),m_sColorMap));
            coeffs_image->setPixel(x, it,  color.rgb());
        }
    }

    *coeffs_image = coeffs_image->scaled(3, m_matSpectrumData.cols()/2, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
    *coeffs_image = coeffs_image->scaledToHeight(image_to_tf_plot->height(), Qt::SmoothTransformation);

    //setup x-axis
    m_x_axis_name = new QGraphicsTextItem("time [sec]", m_pFSPixmap);
    m_x_axis_name->setFont(QFont("arial", 14));

    qreal scaleXText = (m_matSpectrumData.cols() - 1) /  m_dXSampleRate / 20.0;  // divide signallength

    for(qint32 j = 0; j < 21; j++) {
        QGraphicsTextItem* text_item = new QGraphicsTextItem(QString::number(j * scaleXText, 'f', 2), m_pFSPixmap);
        text_item->setFont(QFont("arial", 10));
        m_x_axis_values.append(text_item);    // scalevalue as string
        QGraphicsLineItem* x_line_item = new QGraphicsLineItem(m_pFSPixmap);
        x_line_item->setLine(0,-3,0,3);
        x_line_item->setPen(QPen(Qt::darkGray, 2, Qt::SolidLine, Qt::SquareCap, Qt::MiterJoin));
        m_x_axis_lines.append(x_line_item);                    // scalelines
    }

    m_x_axis_name->setPos(m_pFSPixmap->boundingRect().width()/2 - m_x_axis_name->boundingRect().width()/2,
                          m_pFSPixmap->boundingRect().height() + 0.8 * m_x_axis_values.at(0)->boundingRect().height());

    qreal scale_x = qreal(m_pFSPixmap->boundingRect().width()) / qreal(m_x_axis_values.length()-1);

    for(qint32 i = 0; i < m_x_axis_values.length(); i++) {
        m_x_axis_values.at(i)->setPos(qreal(i)*scale_x - m_x_axis_values.at(0)->boundingRect().width()/2,
                                      m_pFSPixmap->boundingRect().height());
        m_x_axis_lines.at(i)->setPos(qreal(i)*scale_x,
                                     m_pFSPixmap->boundingRect().height());

    }//end x axis

    //y-axis
    m_y_axis_name = new QGraphicsTextItem("frequency [Hz]", m_pFSPixmap);
    m_y_axis_name->setFont(QFont("arial", 14));

    qreal scale_y_text = 0;

    if(m_dMinFreq == 0  && m_dMaxFreq == 0) {
        scale_y_text = 0.5* m_dYSampleRate / 10.0;                       // divide signallength
    } else {
        scale_y_text = (m_dMaxFreq - m_dMinFreq) / 10.0;
    }

    for(qint32 j = 0; j < 11; j++) {
        QGraphicsTextItem* text_item = new QGraphicsTextItem(QString::number(m_dMinFreq + j*scale_y_text,//pow(10, j)/pow(10, 11) * max_frequency,//scale_y_text,
                                                                             'f', 1), m_pFSPixmap);
        text_item->setFont(QFont("arial", 10));
        m_y_axis_values.append(text_item);    // scalevalue as string
        QGraphicsLineItem* y_line_item = new QGraphicsLineItem(m_pFSPixmap);
        y_line_item->setLine(-3,0,3,0);
        y_line_item->setPen(QPen(Qt::darkGray, 2, Qt::SolidLine, Qt::SquareCap, Qt::MiterJoin));
        m_y_axis_lines.append(y_line_item);                    // scalelines
    }

    m_y_axis_name->setPos(- m_x_axis_values.at(0)->boundingRect().width() - m_y_axis_name->boundingRect().height(),
                          m_pFSPixmap->boundingRect().height()/2 + m_y_axis_name->boundingRect().height()/2);
    m_y_axis_name->setRotation(-90);

    qreal scale_y = qreal(m_pFSPixmap->boundingRect().height()) / qreal(m_y_axis_values.length()-1);

    for(qint32 i = 0; i < m_y_axis_values.length(); i++) {
        m_y_axis_values.at(i)->setPos( -m_y_axis_values.last()->boundingRect().width()
                                       -0.5*m_y_axis_lines.last()->boundingRect().width()
                                       -1
                                       ,m_pFSPixmap->boundingRect().height()
                                       - m_y_axis_values.last()->boundingRect().height()/2
                                       - qreal(i)*scale_y);
        m_y_axis_lines.at(i)->setPos(0, qreal(i)*scale_y);
    }
    //end y axis

    m_pItemPixmap = m_pFSScene->addPixmap(QPixmap::fromImage(*coeffs_image));//addItem();
    m_pItemPixmap->setParentItem(m_pFSPixmap);
    m_pItemPixmap->setPos(m_pFSPixmap->boundingRect().width() +5, 0);

    QGraphicsSimpleTextItem *axis_name_item = new QGraphicsSimpleTextItem("Colorbar", m_pItemPixmap);
    m_axis_zero_item = new QGraphicsSimpleTextItem("0", m_pItemPixmap);
    m_axis_one_item = new QGraphicsSimpleTextItem("1", m_pItemPixmap);
    axis_name_item->setFont(QFont("arial", 14));
    m_axis_zero_item->setFont(QFont("arial", 10));
    m_axis_one_item->setFont(QFont("arial", 10));

    axis_name_item->setPos(m_pItemPixmap->boundingRect().width() + 1,
                           m_pItemPixmap->boundingRect().height()/2 + axis_name_item->boundingRect().height()/2);
    axis_name_item->setRotation(-90);
    m_axis_zero_item->setPos( 1 + m_pItemPixmap->boundingRect().width(),
                            m_pItemPixmap->boundingRect().height()- m_axis_zero_item->boundingRect().height());
    m_axis_one_item->setPos( 1 + m_pItemPixmap->boundingRect().width(), 0);
    //end coeffs picture

    m_pView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    m_pView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    m_pView->fitInView(m_pFSScene->sceneRect(),Qt::KeepAspectRatio);
    m_pLayout->addWidget(m_pView);
    this->setLayout(m_pLayout);
    //m_pLayout->update();

    m_bInitialized = true;
}

//=============================================================================================================
void RealTimeSpectrumWidgetNew::initSettingsWidget()
{
    QString sFSName = m_pFS->getName();
    //Init control widgets
    QList<QWidget*> lControlWidgets;

    RtSpectrumViewSettings* pRtSpectrumViewSettings = new RtSpectrumViewSettings(QString("FS/%1").arg(m_pFS->getName()), (m_iXLimit-1)/m_dXSampleRate, 1.0/m_dXSampleRate);

    pRtSpectrumViewSettings->setObjectName("group_tab_View_Spectrum");
    lControlWidgets.append(pRtSpectrumViewSettings);

    connect(pRtSpectrumViewSettings, &RtSpectrumViewSettings::timeWindowChanged,this, &RealTimeSpectrumWidgetNew::setXLimitSeconds);
    connect(pRtSpectrumViewSettings, &RtSpectrumViewSettings::fixColormapChanged,this, &RealTimeSpectrumWidgetNew::setColormapLimits);
    connect(pRtSpectrumViewSettings, &RtSpectrumViewSettings::colormapMaxChanged,this, &RealTimeSpectrumWidgetNew::setColormapMax);
    connect(pRtSpectrumViewSettings, &RtSpectrumViewSettings::colormapMinChanged,this, &RealTimeSpectrumWidgetNew::setColormapMin);

    connect(this, &RealTimeSpectrumWidgetNew::colormapMaxChanged, pRtSpectrumViewSettings, &RtSpectrumViewSettings::setColormapMax);
    connect(this, &RealTimeSpectrumWidgetNew::colormapMinChanged, pRtSpectrumViewSettings, &RtSpectrumViewSettings::setColormapMin);

    emit displayControlWidgetsChanged(lControlWidgets, sFSName);
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::setXLimit(int value)
{
    QMutexLocker locker(&m_qMutex);

    Eigen::MatrixXd matSpectrumDataTemp = Eigen::MatrixXd::Zero(m_matSpectrumData.rows(),value);
    if(value > m_iXLimit)
    {
        matSpectrumDataTemp.leftCols(m_iXLimit) = m_matSpectrumData;
    } else if(m_iDisplayIndex < value)
    {
        matSpectrumDataTemp = m_matSpectrumData.leftCols(value);
    } else
    {
        matSpectrumDataTemp = m_matSpectrumData.block(0,m_matSpectrumData.cols()-value-1,m_matSpectrumData.rows(),value);
        m_iDisplayIndex -= (m_matSpectrumData.cols()-value);
    }
    m_matSpectrumData = matSpectrumDataTemp;
    m_iXLimit = value;

    if (m_iDisplayIndex > m_iXLimit - 1)
        m_iDisplayIndex = 0;

    updateData();
    updateCoeffsPixmap();
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::setXLimitSeconds(double value)
{
    int nSamples = std::max(static_cast<int>(static_cast<double>(m_dXSampleRate)*value) + 1,1);
    setXLimit(nSamples);
}


//=============================================================================================================

void RealTimeSpectrumWidgetNew::setColorMap(const QString &p_sColorMap)
{
    m_sColorMap = p_sColorMap;

    updateCoeffsPixmap();
    updateDataPixmap();
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::setColormapLimits(bool bColormapLimits)
{
    m_bColormapLimits = bColormapLimits;

    updateColorbarLimitsPixmap();
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::setColormapMax(double dColormapMax)
{
    QMutexLocker locker(&m_qMutex);

    m_dMaxLimit = dColormapMax;

    if(m_bColormapLimits)
    {
        updateData();
        updateCoeffsPixmap();
        updateDataPixmap();
        updateColorbarLimitsPixmap();
    }
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::setColormapMin(double dColormapMin)
{
    QMutexLocker locker(&m_qMutex);

    m_dMinLimit = dColormapMin;

    if(m_bColormapLimits)
    {
        updateData();
        updateCoeffsPixmap();
        updateDataPixmap();
        updateColorbarLimitsPixmap();
    }
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::updateDataPixmap()
{
    QImage * image_to_tf_plot = new QImage(m_matCentNormData.cols(), m_matCentNormData.rows(), QImage::Format_RGB32);

    //setup pixelcolors in image
    QColor color;
    for ( qint32 y = 0; y < m_matCentNormData.rows(); y++ ) {
        for ( qint32 x = 0; x < m_matCentNormData.cols(); x++ ) {
            color.setRgb(pColorMapper(m_matCentNormData(y, x),m_sColorMap));
            image_to_tf_plot->setPixel(x, m_matCentNormData.rows() - 1 - y,  color.rgb());
        }
    }

    *image_to_tf_plot = image_to_tf_plot->scaled(m_pFSPixmap->boundingRect().width(), m_pFSPixmap->boundingRect().height(), Qt::IgnoreAspectRatio, Qt::FastTransformation);

    m_pFSPixmap->setPixmap(QPixmap::fromImage(*image_to_tf_plot));
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::updateCoeffsPixmap()
{
    qreal scaleXText = (m_matSpectrumData.cols() - 1) /  (m_dXSampleRate * 20.0);  // divide signallength

    for(qint32 j = 0; j < 21; j++) {
        m_x_axis_values[j]->setPlainText(QString::number(j * scaleXText, 'f', 2));    // scalevalue as string
    }

    qreal scale_y_text = 0;

    //y-axis
    if(m_dMinFreq == 0  && m_dMaxFreq == 0) {
        scale_y_text = 0.5* m_dYSampleRate / 10.0;                       // divide signallength
    } else {
        scale_y_text = (m_dMaxFreq - m_dMinFreq) / 10.0;
    }

    for(qint32 j = 0; j < 11; j++) {
        m_y_axis_values[j]->setPlainText(QString::number(m_dMinFreq + j*scale_y_text,//pow(10, j)/pow(10, 11) * max_frequency,//scale_y_text,
                                                         'f', 1));    // scalevalue as string
    }
    //end y axis
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::updateColorbarLimitsPixmap()
{
    m_axis_zero_item->setText(QString::number(m_dMinLimit, 'e', 2));
    m_axis_one_item->setText(QString::number(m_dMaxLimit, 'e', 2));
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::updateData()
{
    if(m_matSpectrumData.col(m_iDisplayIndex).array().isNaN().any())
    {
        qDebug() << "RealTimeSpectrumWidgetNew::updateData(): NaN";
        m_matSpectrumData.col(m_iDisplayIndex) = Eigen::RowVectorXd::Zero(m_matSpectrumData.rows());
    }
    m_dMinValue = m_matSpectrumData.minCoeff();
    m_dMaxValue = m_matSpectrumData.maxCoeff();

    m_matCentNormData = m_matSpectrumData;

    if(m_bColormapLimits)
    {
        if(m_dMaxValue > m_dMaxLimit)
            m_matCentNormData = m_matCentNormData.cwiseMin(m_dMaxLimit);
        if(m_dMinValue < m_dMinLimit)
            m_matCentNormData = m_matCentNormData.cwiseMax(m_dMinLimit);
        double scale = m_dMaxLimit - m_dMinLimit;
        m_matCentNormData.array() -= m_dMinLimit;
        m_matCentNormData *= 1.0/scale;
    } else
    {
        //double scale = m_dMaxValue - m_dMinValue;
        m_matCentNormData.array() -= m_dMinValue;
        m_matCentNormData *= 1.0/(m_dMaxValue - m_dMinValue);

        m_dMinLimit = m_dMinValue;
        m_dMaxLimit = m_dMaxValue;
        updateColorbarLimitsPixmap();

        emit colormapMaxChanged(m_dMaxLimit);
        emit colormapMinChanged(m_dMinLimit);
    }
}

//=============================================================================================================

void RealTimeSpectrumWidgetNew::resizeEvent(QResizeEvent* event)
{
    if(m_pView && m_pLayout)
        m_pView->fitInView(m_pFSScene->sceneRect(),Qt::KeepAspectRatio);
}
