#ifndef __qSlicerAnnotationModuleWidget_h
#define __qSlicerAnnotationModuleWidget_h

// SlicerQT includes
#include "qSlicerAbstractModuleWidget.h"
#include "qSlicerAnnotationsModuleExport.h"
#include "GUI/qSlicerAnnotationModulePropertyDialog.h"

// CTK includes
#include <ctkPimpl.h>
#include <ctkVTKObject.h>


class qSlicerAnnotationModuleReportDialog;
class qSlicerAnnotationModuleSnapShotDialog;
class vtkSlicerNodeSelectorWidget;
class vtkMRMLAnnotationRulerNode;
class vtkMRMLAnnotationAngleNode;
class qSlicerAnnotationModuleWidgetPrivate;
class QTableWidgetItem;
class QLineEdit;
class QTextDocument;
class QGridLayout;
class vtkMRMLAnnotationLinesNode;
class vtkMRMLAnnotationControlPointsNode;
class vtkMRMLAnnotationNode;
class vtkMRMLAnnotationFiducialNode;
class vtkMRMLNode;
class vtkMRMLAnnotationDisplayNode;
class vtkMRMLAnnotationLineDisplayNode;
class vtkMRMLAnnotationPointDisplayNode;

class vtkMRMLInteractionNode;
class vtkObserverManager;

/// \ingroup Slicer_QtModules_Annotation
class Q_SLICER_QTMODULES_ANNOTATIONS_EXPORT qSlicerAnnotationModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT
  QVTK_OBJECT
public:
    typedef qSlicerAbstractModuleWidget Superclass;
    qSlicerAnnotationModuleWidget(QWidget *parent=0);
    ~qSlicerAnnotationModuleWidget();


  


    //BTX
    /// Different Annotation Types
    enum
      {
        TextNode = 1000,
        FiducialNode,
        AngleNode,
        StickyNode,
        SplineNode,
        RulerNode,
        BidimensionalNode,
        ROINode,
        Screenshot,
      };
    //ETX

    /// the logic observes the interaction node, update the gui to keep in
    /// step with the mouse modes tool bar. If interactionNode is null, try to
    /// get it from the scene.
    void updateWidgetFromInteractionMode(vtkMRMLInteractionNode *interactionNode);

protected:

public slots:
    void refreshTree();

    /// a public slot allowing other modules to open up the screen capture
    /// dialog
    void grabSnapShot();

    /// Update the label showing the active annotation hierarchy, triggered from
    /// the logic modifying the active hierarchy node
    void updateActiveHierarchyLabel();

protected slots:

    /// User clicked on property edit button
    void propertyEditButtonClicked(QString mrmlId);

    // Table and Property Modify
    void moveDownSelected();
    void moveUpSelected();  


  //------------------------------------------------------------------
  // Daniel's approved code starting here

  /// Add a new hierarchy.
  void onAddHierarchyButtonClicked();

  /// Jump the 2d Slices to the first control point of the selected annotation
  void onJumpSlicesButtonClicked();

  /// Select all annotations
  void selectAllButtonClicked();

  /// Unselect all annotations
  void unselectAllButtonClicked();

  /// Toggle the visibility of selected Annotations
  void visibleSelectedButtonClicked();

  /// Un-/Lock selected Annotations
  void lockSelectedButtonClicked();

  /// Delete selected Annotations
  void deleteSelectedButtonClicked();

  /// Make active hierarchy annotations visible/invisible
  void invisibleHierarchyButtonClicked();
  void visibleHierarchyButtonClicked();
  /// Un/Lock annotations in active hierarchy
  void lockHierarchyButtonClicked();
  void unlockHierarchyButtonClicked();

  // Property dialog
  void propertyRestored();
  void propertyAccepted();

  // Report dialog
  void reportDialogRejected();
  void reportDialogAccepted();
  void onReportButtonClicked();

protected:
  QScopedPointer<qSlicerAnnotationModuleWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerAnnotationModuleWidget);
  Q_DISABLE_COPY(qSlicerAnnotationModuleWidget);

  virtual void setup();

  qSlicerAnnotationModulePropertyDialog* m_PropertyDialog;
  qSlicerAnnotationModuleReportDialog* m_ReportDialog;

  /// Type of current Annotations - described by enum
  int m_CurrentAnnotationType;

};

#endif
