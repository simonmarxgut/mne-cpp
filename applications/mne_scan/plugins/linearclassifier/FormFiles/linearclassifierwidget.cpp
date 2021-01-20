//=============================================================================================================
/**
 * @file     linearclassifierwidget.cpp
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
 * @brief    Definition of the LinearClassifierWidget class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "linearclassifierwidget.h"
#include "../ui_linearclassifierwidget.h"

#include "../linearclassifier.h"

#include <utils/ioutils.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QSettings>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <Eigen/Core>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace LINEARCLASSIFIERPLUGIN;
using namespace UTILSLIB;
using namespace Eigen;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

LinearClassifierWidget::LinearClassifierWidget(LinearClassifier* pLinearClassifier, const QString& sSettingsPath, bool bIsRunning, QWidget *parent)
: QWidget(parent)
, ui(new Ui::LinearClassifierWidget)
, m_pLinearClassifier(pLinearClassifier)
, m_sSettingsPath(sSettingsPath)
, m_iInputs(1)
, m_iOutputs(1)
, m_bIsRunning(bIsRunning)
{
    ui->setupUi(this);

    ui->tableWidget_Manual->setItemDelegate(new SpinBoxDelegate());

    ui->tableWidget_Manual->setRowCount(m_iInputs);
    ui->tableWidget_Manual->setColumnCount(m_iOutputs);

    QString header = QString("Output 0");
    ui->tableWidget_Manual->setVerticalHeaderItem(0, new QTableWidgetItem(header));
    header = QString("Input 0");
    ui->tableWidget_Manual->setHorizontalHeaderItem(0, new QTableWidgetItem(header));

    QTableWidgetItem* tempItem = new QTableWidgetItem(tr("0"));
    tempItem->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
    ui->tableWidget_Manual->setItem(0,0,tempItem);

    ui->spinBox_Inputs->setValue(m_iInputs);
    ui->spinBox_Outputs->setValue(m_iOutputs);

    if(!m_bIsRunning){
        connect(ui->spinBox_Inputs,&QSpinBox::editingFinished,this,&LinearClassifierWidget::onChangeInputs);
        connect(ui->spinBox_Outputs,&QSpinBox::editingFinished,this,&LinearClassifierWidget::onChangeOutputs);
    } else {
        ui->spinBox_Inputs->setEnabled(false);
        ui->spinBox_Outputs->setEnabled(false);
    }

    connect(ui->pushButton_InputFile,&QPushButton::released,this,&LinearClassifierWidget::onChangeInputFile);
    connect(ui->pushButton_OutputFile,&QPushButton::released,this,&LinearClassifierWidget::onChangeOutputFile);
    connect(ui->pushButton_SaveOutputFile,&QPushButton::released,this,&LinearClassifierWidget::saveOutputMatrix);

    connect(ui->tableWidget_Manual,&QTableWidget::cellChanged,this,&LinearClassifierWidget::onChangeMatrixWidget);

    loadSettings(m_sSettingsPath);
}

//=============================================================================================================

LinearClassifierWidget::~LinearClassifierWidget()
{
    saveSettings(m_sSettingsPath);

    delete ui;
}

//=============================================================================================================

void LinearClassifierWidget::saveSettings(const QString& settingsPath)
{
    if(settingsPath.isEmpty()) {
        return;
    }

    QSettings settings;

    //settings.setValue(settingsPath + QString("/modality"), m_iModality);
    //settings.setValue(settingsPath + QString("/method"), m_iMethod);
}

//=============================================================================================================

void LinearClassifierWidget::loadSettings(const QString& settingsPath)
{
    if(settingsPath.isEmpty()) {
        return;
    }

    QSettings settings;

    //m_iModality = settings.value(settingsPath + QString("/modality"), m_iModality).toInt();
    //m_iMethod = settings.value(settingsPath + QString("/method"), m_iMethod).toInt();
}

//=============================================================================================================

void LinearClassifierWidget::onChangeInputFile()
{
    QString path = QFileDialog::getOpenFileName(this,
                                                "Select linear classifier matrix file",
                                                "resources/mne_scan/plugins/rereferencing/loc_files",
                                                 tr("Simple text files (*.*)"));

    if(path==NULL){
        path = ui->lineEdit_InputFile->text();
    }

    ui->lineEdit_InputFile->setText(path);

    m_sInputMatrixFilename = path;

    if(readInputMatrix())
        onChangeMatrixExternal();
}

//=============================================================================================================

void LinearClassifierWidget::onChangeOutputFile()
{
    QString path = QFileDialog::getSaveFileName(this,
                                                "Select file to store Linear Classifier matrix",
                                                "resources/mne_scan/plugins/linearclassifier/loc_files",
                                                 tr("Simple text files (*.*)"));

    if(path==NULL){
        path = ui->lineEdit_OutputFile->text();
    }

    ui->lineEdit_OutputFile->setText(path);

    m_sOutputMatrixFilename = path;

    saveOutputMatrix();
}

//=============================================================================================================

void LinearClassifierWidget::onChangeMatrixWidget()
{
    m_pClassifierMatrix = Eigen::MatrixXd::Zero(m_iOutputs,m_iInputs);

    for(int i=0; i<m_iOutputs; ++i)
        for(int j=0; j<m_iInputs; ++j)
            if(ui->tableWidget_Manual->item(i,j))
                m_pClassifierMatrix(i,j) = ui->tableWidget_Manual->item(i,j)->text().toDouble();
            else
                m_pClassifierMatrix(i,j) = 0.0;

    m_pLinearClassifier->setLinearClassifierMatrix(m_pClassifierMatrix);
}

//=============================================================================================================

void LinearClassifierWidget::onChangeMatrixExternal()
{
    m_pClassifierMatrix = Eigen::MatrixXd::Zero(m_iOutputs,m_iInputs);

    for(int i=0; i<m_iOutputs; ++i)
        for(int j=0; j<m_iInputs; ++j)
            ui->tableWidget_Manual->item(i,j)->setText(QString("%1").arg(m_pClassifierMatrix(i,j)));

    m_pLinearClassifier->setLinearClassifierMatrix(m_pClassifierMatrix);
}

//=============================================================================================================

void LinearClassifierWidget::onChangeInputs()
{
    m_iInputs = ui->spinBox_Inputs->value();

    //emit changeClassifierInputs(m_iInputs);
    m_pLinearClassifier->setNInputChannels(m_iInputs);

    changeMatrixColumns(m_iInputs);
}

//=============================================================================================================

void LinearClassifierWidget::onChangeOutputs()
{
    m_iOutputs = ui->spinBox_Outputs->value();

    //emit changeClassifierOutputs(m_iOutputs);
    m_pLinearClassifier->setNOutputChannels(m_iOutputs);

    changeMatrixRows(m_iOutputs);
}

//=============================================================================================================

void LinearClassifierWidget::changeMatrixRows(int value)
{
    int iRowsPre = ui->tableWidget_Manual->rowCount();

    {
        const QSignalBlocker blocker(ui->tableWidget_Manual);

        ui->tableWidget_Manual->setRowCount(value);

        if(value > iRowsPre)
        {
            for(int i=std::max(iRowsPre,0); i<value; i++)
            {
                QString header = QString("Output %1").arg(i);
                ui->tableWidget_Manual->setVerticalHeaderItem(i, new QTableWidgetItem(header));

                for(int j=0; j<ui->tableWidget_Manual->columnCount(); ++j){
                    QTableWidgetItem* tempItem =  new QTableWidgetItem(tr("0"));
                    tempItem->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
                    ui->tableWidget_Manual->setItem(i,j,tempItem);
                }
            }
        }
    }

    onChangeMatrixWidget();
}

//=============================================================================================================

void LinearClassifierWidget::changeMatrixColumns(int value)
{
    int iColumnsPre = ui->tableWidget_Manual->columnCount();

    {
        const QSignalBlocker blocker(ui->tableWidget_Manual);

        ui->tableWidget_Manual->setColumnCount(value);

        if(value > iColumnsPre)
        {
            for(int i=std::max(iColumnsPre,0); i<value; i++)
            {
                QString header = QString("Input %1").arg(i);
                ui->tableWidget_Manual->setHorizontalHeaderItem(i, new QTableWidgetItem(header));

                for(int j=0; j<ui->tableWidget_Manual->rowCount(); ++j){
                    QTableWidgetItem* tempItem =  new QTableWidgetItem(tr("0"));
                    tempItem->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
                    ui->tableWidget_Manual->setItem(j,i,tempItem);
                }
            }
        }
    }

    onChangeMatrixWidget();
}

//=============================================================================================================

Eigen::MatrixXd LinearClassifierWidget::getClassifierMatrix()
{
    return m_pClassifierMatrix;
}

//=============================================================================================================

bool LinearClassifierWidget::readInputMatrix()
{
    Eigen::MatrixXd tempLinearClassifier;

    if(!IOUtils::read_eigen_matrix(tempLinearClassifier,m_sInputMatrixFilename))
    {
        qWarning() << "[LinearClassifier::readInputMatrix] Could not read linear classifier matrix file " << m_sInputMatrixFilename;
        return false;
    } else if(tempLinearClassifier.rows()!=m_iOutputs || tempLinearClassifier.cols()!=m_iInputs)
    {
        qWarning() << "[LinearClassifier::readInputMatrix] Linear classifier matrix does not match number of input or output channels";
        return false;
    }

    m_pClassifierMatrix = tempLinearClassifier;

    return true;
}

//=============================================================================================================

bool LinearClassifierWidget::saveOutputMatrix()
{
    return !IOUtils::write_eigen_matrix(m_pClassifierMatrix,m_sOutputMatrixFilename);
}

//=============================================================================================================

SpinBoxDelegate::SpinBoxDelegate(QObject *parent)
    : QStyledItemDelegate(parent)
{
}

//=============================================================================================================

QWidget *SpinBoxDelegate::createEditor(QWidget *parent,
                                       const QStyleOptionViewItem &/* option */,
                                       const QModelIndex &/* index */) const
{
    QDoubleSpinBox *editor = new QDoubleSpinBox(parent);
    editor->setFrame(false);
    editor->setButtonSymbols(QAbstractSpinBox::NoButtons);
    editor->setDecimals(4);
    //editor->setValue(0.0);
    editor->setSingleStep(std::numeric_limits<double>::epsilon());
    editor->setMinimum(std::numeric_limits<double>::min());
    editor->setMaximum(std::numeric_limits<double>::max());

    return editor;
}

//=============================================================================================================

void SpinBoxDelegate::setEditorData(QWidget *editor,
                                    const QModelIndex &index) const
{
    double value = index.model()->data(index, Qt::EditRole).toDouble();

    QDoubleSpinBox *spinBox = static_cast<QDoubleSpinBox*>(editor);
    spinBox->setValue(value);
}

//=============================================================================================================

void SpinBoxDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                   const QModelIndex &index) const
{
    QDoubleSpinBox *spinBox = static_cast<QDoubleSpinBox*>(editor);
    spinBox->interpretText();
    double value = spinBox->value();

    model->setData(index, value, Qt::EditRole);
}

//=============================================================================================================

void SpinBoxDelegate::updateEditorGeometry(QWidget *editor,
                                           const QStyleOptionViewItem &option,
                                           const QModelIndex &/* index */) const
{
    editor->setGeometry(option.rect);
}
