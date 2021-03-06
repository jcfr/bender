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

#ifndef __vtkMRMLBoneDisplayNode_h
#define __vtkMRMLBoneDisplayNode_h

// Slicer includes
#include "vtkMRMLAnnotationDisplayNode.h"

// Armatures includes
#include "vtkBenderArmaturesModuleMRMLCoreExport.h"

/// \ingroup Bender_MRML
/// \brief Annotation to represent a bone.
///
/// \sa vtkMRMLBoneNode, vtkMRMLArmatureNode
class VTK_BENDER_ARMATURES_MRML_CORE_EXPORT vtkMRMLBoneDisplayNode
  : public vtkMRMLAnnotationDisplayNode
{
public:
  //--------------------------------------------------------------------------
  // VTK methods
  //--------------------------------------------------------------------------

  static vtkMRMLBoneDisplayNode *New();
  vtkTypeMacro(vtkMRMLBoneDisplayNode,vtkMRMLAnnotationDisplayNode);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  //--------------------------------------------------------------------------
  // MRMLNode methods
  //--------------------------------------------------------------------------

  /// Instantiate a bone node.
  virtual vtkMRMLNode* CreateNodeInstance();

  /// Get node XML tag name (like Volume, Model).
  virtual const char* GetNodeTagName() {return "BoneDisplay";};

  /// Read node attributes from XML file.
  virtual void ReadXMLAttributes( const char** atts);

  /// Write this node's information to a MRML file in XML format.
  virtual void WriteXML(ostream& of, int indent);

  /// Copy the node's attributes to this object.
  virtual void Copy(vtkMRMLNode *node);

  virtual void UpdateScene(vtkMRMLScene *scene);
  virtual void ProcessMRMLEvents(vtkObject* caller,
                                 unsigned long event,
                                 void* callData);


  //--------------------------------------------------------------------------
  // Bone methods
  //--------------------------------------------------------------------------

protected:
  vtkMRMLBoneDisplayNode();
  ~vtkMRMLBoneDisplayNode();

  vtkMRMLBoneDisplayNode(const vtkMRMLBoneDisplayNode&); /// not implemented
  void operator=(const vtkMRMLBoneDisplayNode&); /// not implemented
};

#endif
