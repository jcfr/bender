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

//Rip off vtkLineAnnotation

#ifndef __vtkArmatureRepresentation_h
#define __vtkArmatureRepresentation_h

// VTK includes
#include "vtkWidgetRepresentation.h"

// Bender includes
#include "vtkBenderWidgetsExport.h"

class vtkPointHandleRepresentation3D;

class VTK_BENDER_WIDGETS_EXPORT vtkArmatureRepresentation
  : public vtkWidgetRepresentation
{
public:
  // Description:
  // Instantiate the class.
  static vtkArmatureRepresentation *New();

  // Description:
  // Standard methods for the class.
  vtkTypeMacro(vtkArmatureRepresentation,vtkWidgetRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void BuildRepresentation();

protected:
  vtkArmatureRepresentation();
  ~vtkArmatureRepresentation();

private:
  vtkArmatureRepresentation(const vtkArmatureRepresentation&);  //Not implemented
  void operator=(const vtkArmatureRepresentation&);  //Not implemented
};

#endif
