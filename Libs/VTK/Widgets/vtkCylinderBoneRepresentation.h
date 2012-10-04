/*=========================================================================

  Program: Bender

  Copyright (c) Kitware Inc.

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

=========================================================================*/

#ifndef __vtkCylinderBoneRepresentation_h
#define __vtkCylinderBoneRepresentation_h

// Bender includes
#include "vtkBoneRepresentation.h"

class vtkActor;
class vtkPolyDataMapper;
class vtkPolyData;
class vtkProperty;
class vtkTubeFilter;

class VTK_BENDER_WIDGETS_EXPORT vtkCylinderBoneRepresentation
  : public vtkBoneRepresentation
{
public:

  // Description:
  // Instantiate this class.
  static vtkCylinderBoneRepresentation *New();

  // Description:
  // Standard methods for the class.
  vtkTypeMacro(vtkCylinderBoneRepresentation, vtkBoneRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // These are methods that satisfy vtkWidgetRepresentation's API.
  virtual void BuildRepresentation();
  virtual double *GetBounds();

  // Description:
  // Get the radius of the cylinder. The radius is automatically computed
  // from the distance between the two endpoints.
  vtkGetMacro(Radius,double);

  // Description:
  // Set/Get if the cylinder is capped or not. Default is true.
  vtkSetMacro(Capping, int);
  vtkGetMacro(Capping, int);

  // Description:
  // Set/Get the number of sides of the cylinder. The minimum is 3 and
  // the default is 5.
  vtkSetMacro(NumberOfSides, int);
  vtkGetMacro(NumberOfSides, int);

  // Description:
  // Get the cylinder properties. The properties of the cylinder
  // can be manipulated.
  vtkGetObjectMacro(CylinderProperty,vtkProperty);

  // Description:
  // Methods supporting the rendering process.
  virtual void GetActors(vtkPropCollection *pc);
  virtual void ReleaseGraphicsResources(vtkWindow*);
  virtual int RenderOpaqueGeometry(vtkViewport*);
  virtual int RenderTranslucentPolygonalGeometry(vtkViewport*);
  virtual int HasTranslucentPolygonalGeometry();
  
protected:
  vtkCylinderBoneRepresentation();
  ~vtkCylinderBoneRepresentation();

  // the Cylinder
  vtkActor*          CylinderActor;
  vtkPolyDataMapper* CylinderMapper;
  vtkTubeFilter*     CylinderGenerator;

  // Properties used to control the appearance of selected objects and
  // the manipulator in general.
  vtkProperty* CylinderProperty;
  void         CreateDefaultProperties();

  //Cylinder properties
  double Radius;
  int    Capping;
  int    NumberOfSides;

  void RebuildCylinder();

private:
  vtkCylinderBoneRepresentation(const vtkCylinderBoneRepresentation&);  //Not implemented
  void operator=(const vtkCylinderBoneRepresentation&);  //Not implemented
};

#endif
