//=============================================================================================================
/**
 * @file     rereferencesetupwidget.h
 * @author   Christoph Dinh <chdinh@nmr.mgh.harvard.edu>;
 *           Viktor Klueber <Viktor.Klueber@tu-ilmenau.de>
 * @since    0.1.0
 * @date     February, 2013
 *
 * @section  LICENSE
 *
 * Copyright (C) 2013, Christoph Dinh, Viktor Klueber. All rights reserved.
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
 * @brief    Contains the declaration of the RereferenceSetupWidget class.
 *
 */

#ifndef REREFERENCESETUPWIDGET_H
#define REREFERENCESETUPWIDGET_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../rereference.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QtWidgets>

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

namespace Ui {
    class RereferenceSetupWidgetClass;
}

//=============================================================================================================
// DEFINE NAMESPACE REREFERENCETOOLBOXPLUGIN
//=============================================================================================================

namespace REREFERENCEPLUGIN
{

//=============================================================================================================
// REREFERENCETOOLBOXPLUGIN FORWARD DECLARATIONS
//=============================================================================================================

class Rereference;

//=============================================================================================================
/**
 * DECLARE CLASS RereferenceSetupWidget
 *
 * @brief The RereferenceSetupWidget class provides the Rereference-plugin configuration window.
 */
class RereferenceSetupWidget : public QWidget
{
    Q_OBJECT

public:

    //=========================================================================================================
    /**
     * Constructs a RereferenceSetupWidget which is a child of parent.
     *
     * @param [in] toolbox a pointer to the corresponding Rereference-plugin.
     * @param [in] parent pointer to parent widget; If parent is 0, the new RereferenceSetupWidget becomes a window. If parent is another widget, RereferenceSetupWidget becomes a child window inside parent. RereferenceSetupWidget is deleted when its parent is deleted.
     */
    RereferenceSetupWidget(Rereference* pRereference, QWidget *parent = 0);

    //=========================================================================================================
    /**
     * Destroys the RereferenceSetupWidget.
     * All RereferenceSetupWidget's children are deleted first. The application exits if RereferenceSetupWidget is the main widget.
     */
    ~RereferenceSetupWidget();

private:
    Rereference*               m_pRereference;	/**< Holds a pointer to corresponding Rereference-plugin.*/

    Ui::RereferenceSetupWidgetClass*  ui;              /**< Holds the user interface for the RereferenceSetupWidget.*/
};
} // NAMESPACE

#endif // REREFERENCESETUPWIDGET_H
