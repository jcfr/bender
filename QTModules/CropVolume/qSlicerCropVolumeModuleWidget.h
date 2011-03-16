#ifndef __qSlicerCropVolumeModuleWidget_h
#define __qSlicerCropVolumeModuleWidget_h

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerCropVolumeModuleExport.h"

class qSlicerCropVolumeModuleWidgetPrivate;
class vtkMRMLNode;
class vtkMRMLCropVolumeParametersNode;

class Q_SLICER_QTMODULES_CROPVOLUME_EXPORT qSlicerCropVolumeModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerCropVolumeModuleWidget(QWidget *parent=0);
  virtual ~qSlicerCropVolumeModuleWidget();

public slots:

protected:
  QScopedPointer<qSlicerCropVolumeModuleWidgetPrivate> d_ptr;
  
  virtual void setup();

protected slots:
  void initializeNode(vtkMRMLNode*);
  void onInputVolumeChanged();
  void onInputROIChanged();
  void onOutputVolumeChanged();
  void onROIVisibilityChanged();
  void onInterpolationModeChanged();
  void onApply();
  void updateWidget();

private:
  Q_DECLARE_PRIVATE(qSlicerCropVolumeModuleWidget);
  Q_DISABLE_COPY(qSlicerCropVolumeModuleWidget);

  vtkMRMLCropVolumeParametersNode *parametersNode;
};

#endif