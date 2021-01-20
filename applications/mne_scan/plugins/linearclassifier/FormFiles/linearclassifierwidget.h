//=============================================================================================================
/**
 * @file     linearclassifierwidget.h
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
 * @brief    Contains the declaration of the LinearClassifierWidget class.
 *
 */

#ifndef LINEARCLASSIFIERWIDGET_H
#define LINEARCLASSIFIERWIDGET_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../linearclassifier.h"

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
    class LinearClassifierWidget;
}

//=============================================================================================================
// DEFINE NAMESPACE LINEARCLASSIFIERPLUGIN
//=============================================================================================================

namespace LINEARCLASSIFIERPLUGIN
{

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
/**
 * DECLARE CLASS LinearClassifierWidget
 *
 * @brief The LinearClassifier class provides a rereference widget.
 */
class LinearClassifierWidget : public QWidget
{
    Q_OBJECT

public:
    typedef QSharedPointer<LinearClassifierWidget> SPtr;         /**< Shared pointer type for RereferenceWidget. */
    typedef QSharedPointer<LinearClassifierWidget> ConstSPtr;    /**< Const shared pointer type for RereferenceWidget. */

    //=========================================================================================================
    /**
     * Constructs a LinearClassifierWidget.
     */
    explicit LinearClassifierWidget(LinearClassifier* pLinearClassifier, const QString& sSettingsPath = "", bool bIsRunning = false,
                             QWidget *parent = 0);

    //=========================================================================================================
    /**
     * Destroys the LinearClassifierWidget.
     */
    ~LinearClassifierWidget();

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
     * Slot called when input classifier matrix file change
     */
    void onChangeInputFile();

    //=========================================================================================================
    /**
     * Slot called when input classifier matrix file change
     */
    void onChangeOutputFile();

    //=========================================================================================================
    /**
     * Slot called when input classifier matrix changed in widget
     */
    void onChangeMatrixWidget();

    //=========================================================================================================
    /**
     * Slot called when input classifier matrix changed externally
     */
    void onChangeMatrixExternal();
    //=========================================================================================================
    /**
     * Slot called when number of inputs changes
     */
    void onChangeInputs();

    //=========================================================================================================
    /**
     * Slot called when number of outputs changes
     */
    void onChangeOutputs();

    //=========================================================================================================
    /**
     * Slot called when number of matrix rows changes
     */
    void changeMatrixRows(int value);

    //=========================================================================================================
    /**
     * Slot called when number of matrix columns changes
     */
    void changeMatrixColumns(int value);

    //=========================================================================================================
    /**
     * Returns the classifier matrix
     */
    Eigen::MatrixXd getClassifierMatrix();

    //=========================================================================================================
    /**
     * Read input matrix
     *
     */
    bool readInputMatrix();

    //=========================================================================================================
    /**
     * Save output matrix
     *
     */
    bool saveOutputMatrix();

    LinearClassifier*               m_pLinearClassifier;	/**< Holds a pointer to corresponding Rereference-plugin.*/

    QString                     m_sSettingsPath;    /**< The settings path to store the GUI settings to. */
    int                         m_iInputs;
    int                         m_iOutputs;
    bool                        m_bIsRunning;
    Eigen::MatrixXd             m_pClassifierMatrix;
    QString                                         m_sOutputMatrixFilename;
    QString                                         m_sInputMatrixFilename;

    Ui::LinearClassifierWidget*     ui;              /**< The UI class specified in the designer. */

signals:
    //=========================================================================================================
    /**
     * Emitted whenever the settings changed and are ready to be retreived.
     */
    void changeClassifierMatrix();
    void changeClassifierInputs(int value);
    void changeClassifierOutputs(int value);
};

class SpinBoxDelegate : public QStyledItemDelegate
{
    Q_OBJECT

public:
    SpinBoxDelegate(QObject *parent = nullptr);

    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                          const QModelIndex &index) const override;

    void setEditorData(QWidget *editor, const QModelIndex &index) const override;
    void setModelData(QWidget *editor, QAbstractItemModel *model,
                      const QModelIndex &index) const override;

    void updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option,
                              const QModelIndex &index) const override;
};
}   //namespace

#endif // LINEARCLASSIFIERWIDGET_H
