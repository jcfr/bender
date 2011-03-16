/*==============================================================================

  Program: 3D Slicer

  Copyright (c) 2010 Kitware Inc.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Julien Finet, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

#ifndef __qSlicerSettingsGeneralPanel_h
#define __qSlicerSettingsGeneralPanel_h

// Qt includes
#include <QWidget>

// CTK includes
#include <ctkSettingsPanel.h>

// QtGUI includes
#include "qSlicerBaseQTGUIExport.h"

class QSettings;
class qSlicerSettingsGeneralPanelPrivate;

class Q_SLICER_BASE_QTGUI_EXPORT qSlicerSettingsGeneralPanel
  : public ctkSettingsPanel
{
  Q_OBJECT
public:
  /// Superclass typedef
  typedef ctkSettingsPanel Superclass;

  /// Constructor
  explicit qSlicerSettingsGeneralPanel(QWidget* parent = 0);

  /// Destructor
  virtual ~qSlicerSettingsGeneralPanel();

protected slots:
  void onFontChanged(const QFont& font);
  void onShowToolTipsToggled(bool);
  void onShowToolButtonTextToggled(bool enable);

protected:
  QScopedPointer<qSlicerSettingsGeneralPanelPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerSettingsGeneralPanel);
  Q_DISABLE_COPY(qSlicerSettingsGeneralPanel);
};

#endif