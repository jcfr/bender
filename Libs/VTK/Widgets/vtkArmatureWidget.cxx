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

// Bender includes
#include "vtkArmatureRepresentation.h"
#include "vtkArmatureWidget.h"
#include "vtkBoneRepresentation.h"
#include "vtkBoneWidget.h"
#include "vtkCylinderBoneRepresentation.h"
#include "vtkDoubleConeBoneRepresentation.h"

// VTK includes
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkCollection.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkWidgetCallbackMapper.h>
#include <vtkWidgetEvent.h>

// STD includes
#include <algorithm>

vtkStandardNewMacro(vtkArmatureWidget);

struct ArmatureTreeNode
{
  vtkBoneWidget* Bone;
  std::vector<ArmatureTreeNode*> Children;
  typedef std::vector<ArmatureTreeNode*>::iterator ChildrenNodeIteratorType;
  ArmatureTreeNode* Parent;
  bool HeadLinkedToParent;

  ArmatureTreeNode()
    {
    this->Bone = 0;
    this->Children.clear();
    this->Parent = 0;
    this->HeadLinkedToParent = false;
    }

  // Add the child to the parent's vector
  // and add the parent to the child's reference
  void AddChild(ArmatureTreeNode* child)
    {
    this->Children.push_back(child);
    child->Parent = this;
    }

  // Remove bone, obsevers and this
  void Delete()
    {
    this->Bone->RemoveAllObservers();
    this->Bone->Delete();

    delete this;
    }

  // Link child's children to parent and vise versa
  // and then remove the child.
  void RemoveChild(ArmatureTreeNode* child)
    {
    // Rebuild linkage
    for (ChildrenNodeIteratorType it = child->Children.begin();
      it != child->Children.end(); ++it)
      {
      this->Children.push_back(*it);
      (*it)->Parent = this;
      }

    // Remove child from parent
    ChildrenNodeIteratorType childPosition =
      std::find(this->Children.begin(), this->Children.end(), child);
    this->Children.erase(childPosition);
    }

  // Find a new root (if any), update children and returns the new root
  ArmatureTreeNode* RemoveRoot()
    {
    ArmatureTreeNode* newRoot = NULL;

    for (ChildrenNodeIteratorType it = this->Children.begin();
      it != this->Children.end(); ++it)
      {
      if (!newRoot)
        {
        (*it)->Parent = NULL;
        newRoot = *it;
        }
      else
        {
        (*it)->Parent = newRoot;
        newRoot->Children.push_back(*it);
        }
      }

    return newRoot;
    }
};

class ArmatureTreeNodeVectorType
  : public std::vector<ArmatureTreeNode*> {};
typedef ArmatureTreeNodeVectorType::iterator NodeIteratorType;
typedef std::vector<ArmatureTreeNode*> ChildrenVectorType;
typedef ChildrenVectorType::iterator ChildrenNodeIteratorType;

class vtkArmatureWidget::vtkArmatureWidgetCallback : public vtkCommand
{
public:
  static vtkArmatureWidgetCallback *New()
    {
    return new vtkArmatureWidgetCallback;
    }
  vtkArmatureWidgetCallback() { this->ArmatureWidget = 0; }
  virtual void Execute(vtkObject* caller, unsigned long eventId, void*)
    {
    switch (eventId)
      {
      case vtkBoneWidget::RestChangedEvent:
        {
        // Assume the cast isn't null. The only elements that should send
        // rest changed event should be vtkBoneWidget
        vtkBoneWidget* bone = vtkBoneWidget::SafeDownCast(caller);

        ArmatureTreeNode* node = ArmatureWidget->GetNode(bone);
        if (node && node->HeadLinkedToParent && node->Parent)
          {
          double parentTail[3];
          node->Parent->Bone->GetCurrentWorldTail(parentTail);
          node->Bone->SetWorldHeadRest(parentTail);
          }

        ArmatureWidget->UpdateChildrenWidgetStateToRest(node);

        break;
        }
      case vtkBoneWidget::PoseChangedEvent:
        {
        // Assume the cast isn't null. The only elements that should send
        // pose changed event should be vtkBoneWidget
        vtkBoneWidget* bone = vtkBoneWidget::SafeDownCast(caller);

        ArmatureTreeNode* node = ArmatureWidget->GetNode(bone);
        ArmatureWidget->UpdateChildrenWidgetStateToPose(node);

        break;
        }
      }
    }

  vtkArmatureWidget *ArmatureWidget;
};

//----------------------------------------------------------------------------
vtkArmatureWidget::vtkArmatureWidget()
{
  // Init callback
  this->ArmatureWidgetCallback =
    vtkArmatureWidget::vtkArmatureWidgetCallback::New();
  this->ArmatureWidgetCallback->ArmatureWidget = this;

  // Init map and root
  this->Bones = new ArmatureTreeNodeVectorType;
  // Init bones properties
  this->BonesRepresentationType = vtkArmatureWidget::None;
  this->WidgetState = vtkArmatureWidget::Rest;
  this->AxesVisibility = vtkBoneWidget::Hidden;
  this->ShowParenthood = true;
}

//----------------------------------------------------------------------------
vtkArmatureWidget::~vtkArmatureWidget()
{
  // Delete all the bones in the map
  for (NodeIteratorType it = this->Bones->begin();
    it != this->Bones->end(); ++it)
    {
    (*it)->Bone->RemoveAllObservers();
    (*it)->Bone->Delete(); // Delete bone
    }

  this->ArmatureWidgetCallback->Delete();
}

//----------------------------------------------------------------------------
void vtkArmatureWidget::CreateDefaultRepresentation()
{
  vtkArmatureRepresentation* armatureRepresentation =
    this->GetArmatureRepresentation();
  if (!armatureRepresentation)
    {
    armatureRepresentation = vtkArmatureRepresentation::New();
    this->WidgetRep = armatureRepresentation;
    }
}

//----------------------------------------------------------------------------
void vtkArmatureWidget::SetEnabled(int enabling)
{
  if (enabling)
    {
    // Set interactor and renderer to all the bones
    for (NodeIteratorType it = this->Bones->begin();
      it != this->Bones->end(); ++it)
      {
      (*it)->Bone->SetInteractor(this->Interactor);
      (*it)->Bone->SetCurrentRenderer(this->CurrentRenderer);
      }
    }

  if (this->WidgetRep)
    {
    this->WidgetRep->SetVisibility(enabling);
    }

  // Enabled all the bones
  for (NodeIteratorType it = this->Bones->begin();
    it != this->Bones->end(); ++it)
    {
    if (!(*it)->Bone->GetBoneRepresentation())
      {
      this->SetBoneRepresentation(
        (*it)->Bone, this->BonesRepresentationType);
      }
    (*it)->Bone->SetEnabled(enabling);
    }

  this->Superclass::SetEnabled(enabling);
}

//----------------------------------------------------------------------------
void vtkArmatureWidget::SetRepresentation(vtkArmatureRepresentation* r)
{
  this->Superclass::SetWidgetRepresentation(r);
}

//----------------------------------------------------------------------------
vtkArmatureRepresentation* vtkArmatureWidget::GetArmatureRepresentation()
{
  return vtkArmatureRepresentation::SafeDownCast(this->WidgetRep);
}

//----------------------------------------------------------------------------
void vtkArmatureWidget::SetProcessEvents(int pe)
{
  this->Superclass::SetProcessEvents(pe);

  for (NodeIteratorType it = this->Bones->begin();
    it != this->Bones->end(); ++it)
    {
    (*it)->Bone->SetProcessEvents(pe);
    }
}

//----------------------------------------------------------------------------
vtkBoneWidget* vtkArmatureWidget
::CreateBone(vtkBoneWidget* parent, const vtkStdString& name)
{
  vtkBoneWidget* newBone = vtkBoneWidget::New();
  newBone->SetName(name);
  this->UpdateBoneWithArmatureOptions(newBone, parent);
  return newBone;
}


//----------------------------------------------------------------------------
vtkBoneWidget* vtkArmatureWidget::CreateBone(vtkBoneWidget* parent,
                                             double tail[3],
                                             const vtkStdString& name)
{
  if (!parent)
    {
    vtkErrorMacro("The bone inserted with this method must have a parent.");
    return NULL;
    }

  vtkBoneWidget* newBone = vtkBoneWidget::New();
  newBone->SetName(name);
  this->UpdateBoneWithArmatureOptions(newBone, parent);
  newBone->SetWorldHeadRest(parent->GetWorldTailRest());
  newBone->SetWorldTailRest(tail);

  return newBone;
}

//----------------------------------------------------------------------------
vtkBoneWidget* vtkArmatureWidget::CreateBone(vtkBoneWidget* parent,
                                             double xTail,
                                             double yTail,
                                             double zTail,
                                             const vtkStdString& name)
{
  double tail[3];
  tail [0] = xTail;
  tail [1] = yTail;
  tail [2] = zTail;
  return this->CreateBone(parent, tail, name);
}

//----------------------------------------------------------------------------
void vtkArmatureWidget
::AddBone(vtkBoneWidget* bone, vtkBoneWidget* parent, bool linkedWithParent)
{
  bone->Register(this);

  // Create bone
  ArmatureTreeNode* newNode = this->CreateNewMapElement(parent);
  newNode->Bone = bone;
  if (! parent)
    {
    this->TopLevelBones.push_back(bone);
    }

  if (! parent || ! linkedWithParent)
    {
    newNode->HeadLinkedToParent = false;
    }
  else
    {
    bone->SetWorldHeadRest(parent->GetWorldTailRest());
    newNode->HeadLinkedToParent = true;
    }

  this->Bones->push_back(newNode);
  newNode->Bone->SetDebugBoneID(this->Bones->size()); // Debug
}

//----------------------------------------------------------------------------
ArmatureTreeNode* vtkArmatureWidget::CreateNewMapElement(vtkBoneWidget* parent)
{
  // Create new bone map element
  ArmatureTreeNode* newNode = new ArmatureTreeNode();

  // Take care of parent stuff
  newNode->Parent = this->GetNode(parent);

  if (newNode->Parent) // Add this to the parent children list
    {
    newNode->Parent->AddChild(newNode);
    }

  return newNode;
}

//----------------------------------------------------------------------------
void vtkArmatureWidget::SetWidgetState(int state)
{
  if (state == this->WidgetState)
    {
    return;
    }

  this->WidgetState = state;
  if (this->TopLevelBones.empty())
    {
    vtkErrorMacro("Could not find any root element !"
      " Cannot set armature state.");
    return;
    }

  switch (this->WidgetState)
    {
    case vtkArmatureWidget::Rest:
      {
      // No need for any recursivity here !
      for (NodeIteratorType it = this->Bones->begin();
        it != this->Bones->end(); ++it)
        {
        (*it)->Bone->SetWidgetStateToRest();
        }
      break;
      }
    case vtkArmatureWidget::Pose:
      {
      // could be smart and check if something has changed since the last
      // state switch and only udpate if it's the case
      for (BoneVectorIterator it = this->TopLevelBones.begin();
        it != this->TopLevelBones.end(); ++it)
        {
        this->SetBoneWorldToParentPoseTransform(*it, NULL);
        }

      for (NodeIteratorType it = this->Bones->begin();
        it != this->Bones->end(); ++it)
        {
        (*it)->Bone->SetWidgetStateToPose();
        }
      break;
      }
    default:
      {
      vtkErrorMacro("Unknown state: "
        <<this->WidgetState<<"/n  ->State unchanged");
      break;
      }
    }
}

//----------------------------------------------------------------------------
void vtkArmatureWidget::SetAxesVisibility(int axesVisibility)
{
  if (axesVisibility == this->AxesVisibility)
    {
    return;
    }

  this->AxesVisibility = axesVisibility;
  for (NodeIteratorType it = this->Bones->begin();
    it != this->Bones->end(); ++it)
    {
    (*it)->Bone->SetAxesVisibility(this->AxesVisibility);
    }
}

//----------------------------------------------------------------------------
void vtkArmatureWidget::SetShowParenthood(int parenthood)
{
  if (this->ShowParenthood == parenthood)
    {
    return;
    }

  this->ShowParenthood = parenthood;
  for (NodeIteratorType it = this->Bones->begin();
    it != this->Bones->end(); ++it)
    {
    (*it)->Bone->SetShowParenthood(this->ShowParenthood);
    }
}

//----------------------------------------------------------------------------
void vtkArmatureWidget::SetBonesRepresentation(int representationType)
{
  if (representationType == this->BonesRepresentationType)
    {
    return;
    }

  if (representationType < 0
    || representationType > vtkArmatureWidget::DoubleCone)
    {
    vtkErrorMacro("Unknown representation type: "
      <<this->BonesRepresentationType<<"/n  ->Representations unchanged");
    }

  this->BonesRepresentationType = representationType;
  for (NodeIteratorType it = this->Bones->begin();
    it != this->Bones->end(); ++it)
    {
    this->SetBoneRepresentation((*it)->Bone, this->BonesRepresentationType);
    }
}

//----------------------------------------------------------------------------
void vtkArmatureWidget
::SetBoneRepresentation(vtkBoneWidget* bone, int representationType)
{
  switch (representationType)
    {
    case vtkArmatureWidget::None:
      {
      bone->SetRepresentation(NULL);

      break;
      }
    case vtkArmatureWidget::Bone:
      {
      vtkSmartPointer<vtkBoneRepresentation> boneRep =
        vtkSmartPointer<vtkBoneRepresentation>::New();
      bone->SetRepresentation(boneRep);

      break;
      }
    case vtkArmatureWidget::Cylinder:
      {
      vtkSmartPointer<vtkCylinderBoneRepresentation> cylinderRep =
          vtkSmartPointer<vtkCylinderBoneRepresentation>::New();
      bone->SetRepresentation(cylinderRep);

      break;
      }
    case vtkArmatureWidget::DoubleCone:
      {
      vtkSmartPointer<vtkDoubleConeBoneRepresentation> simsIconRep =
          vtkSmartPointer<vtkDoubleConeBoneRepresentation>::New();
      bone->SetRepresentation(simsIconRep);

      break;
      }
    default:
      {
      vtkErrorMacro("Unknown representation type: "
        <<representationType<<"/n  ->Representation unchanged");
      break;
      }
    }
}

//----------------------------------------------------------------------------
bool vtkArmatureWidget::HasBone(vtkBoneWidget* bone)
{
  return this->GetNode(bone) != 0;
}

//----------------------------------------------------------------------------
vtkBoneWidget* vtkArmatureWidget::GetBoneParent(vtkBoneWidget* bone)
{
  ArmatureTreeNode* node = this->GetNode(bone);
  if (node && node->Parent)
    {
    return node->Parent->Bone;
    }

  return NULL;
}

//----------------------------------------------------------------------------
vtkCollection* vtkArmatureWidget::FindBoneChildren(vtkBoneWidget* parent)
{
  vtkCollection* children = vtkCollection::New();

  ArmatureTreeNode* node = this->GetNode(parent);
  if (node)
    {
    for (NodeIteratorType it = node->Children.begin();
      it != node->Children.end(); ++it)
      {
      children->AddItem((*it)->Bone);
      }
    }

  return children;
}

//----------------------------------------------------------------------------
bool vtkArmatureWidget::RemoveBone(vtkBoneWidget* bone)
{
  for (NodeIteratorType it = this->Bones->begin();
    it != this->Bones->end(); ++it)
    {
    if ((*it)->Bone == bone)
      {
      if ((*it)->Parent) // Stitch children to the father
        {
        (*it)->Parent->RemoveChild(*it);
        this->UpdateChildren((*it)->Parent);
        }
      else // It was root
        {
        ArmatureTreeNode* replacementRoot = (*it)->RemoveRoot();
        if (replacementRoot)
          {
          this->TopLevelBones.push_back(replacementRoot->Bone);
          this->UpdateChildren(replacementRoot);
          }

        BoneVectorIterator rootsIterator
          = std::find(this->TopLevelBones.begin(),
            this->TopLevelBones.end(),
            (*it)->Bone);
        if (rootsIterator != this->TopLevelBones.end())
          {
          this->TopLevelBones.erase(rootsIterator);
          }
        }

      (*it)->Delete();
      this->Bones->erase(it);
      return true;
      }
    }

  return false;
}

//----------------------------------------------------------------------------
bool vtkArmatureWidget
::SetBoneName(vtkBoneWidget* bone, const vtkStdString& name)
{
  ArmatureTreeNode* node = this->GetNode(bone);
  if (node)
    {
    node->Bone->SetName(name);
    return true;
    }

  return false;
}

//----------------------------------------------------------------------------
vtkStdString vtkArmatureWidget::GetBoneName(vtkBoneWidget* bone)
{
  ArmatureTreeNode* node = this->GetNode(bone);
  if (node)
    {
    return node->Bone->GetName();
    }

  return "";
}

//----------------------------------------------------------------------------
vtkBoneWidget* vtkArmatureWidget::GetBoneByName(const vtkStdString& name)
{
  for (NodeIteratorType it = this->Bones->begin();
    it != this->Bones->end(); ++it)
    {
    if ((*it)->Bone->GetName() == name)
      {
      return (*it)->Bone;
      }
    }
  return NULL;
}

//----------------------------------------------------------------------------
bool vtkArmatureWidget::GetBoneLinkedWithParent(vtkBoneWidget* bone)
{
  ArmatureTreeNode* node = this->GetNode(bone);
  if (node)
    {
    return node->HeadLinkedToParent;
    }
  return false;
}

//----------------------------------------------------------------------------
void vtkArmatureWidget
::SetBoneLinkedWithParent(vtkBoneWidget* bone, bool linked)
{
  ArmatureTreeNode* node = this->GetNode(bone);
  if (node && node->HeadLinkedToParent != linked)
    {
    node->HeadLinkedToParent = linked;
    if (node->Parent)
      {
      node->Bone->SetWorldHeadRest(node->Parent->Bone->GetWorldTailRest());
      }
    }
}

//----------------------------------------------------------------------------
void vtkArmatureWidget::AddBoneObservers(vtkBoneWidget* bone)
{
  bone->AddObserver(vtkBoneWidget::RestChangedEvent,
    this->ArmatureWidgetCallback, this->Priority);
  bone->AddObserver(vtkBoneWidget::PoseChangedEvent,
    this->ArmatureWidgetCallback, this->Priority);
}

//----------------------------------------------------------------------------
void vtkArmatureWidget
::UpdateBoneWithArmatureOptions(vtkBoneWidget* bone, vtkBoneWidget* parent)
{
  this->SetBoneRepresentation(bone, this->BonesRepresentationType);

  if (this->WidgetState == vtkArmatureWidget::Rest)
    {
    this->SetBoneWorldToParentRestTransform(bone, parent);
    bone->SetWidgetStateToRest();
    }
  else
    {
    this->SetBoneWorldToParentPoseTransform(bone, parent);
    bone->SetWidgetStateToPose();
    }

  this->AddBoneObservers(bone);
}

//----------------------------------------------------------------------------
void vtkArmatureWidget
::SetBoneWorldToParentRestTransform(vtkBoneWidget* bone, vtkBoneWidget* parent)
{
  double parentWorldToBoneRestRotation[4] = {1.0, 0.0, 0.0, 0.0};
  double parentWorldToBoneRestTranslation[3] = {0.0, 0.0, 0.0};
  if (parent) // For all non-root elements
    {
    parent->GetWorldToBoneRestRotation(
      parentWorldToBoneRestRotation);
    parent->GetWorldToBoneTailRestTranslation(
      parentWorldToBoneRestTranslation);
    }

  bone->SetWorldToParentRestRotationAndTranslation(
    parentWorldToBoneRestRotation, parentWorldToBoneRestTranslation);
}

//----------------------------------------------------------------------------
void vtkArmatureWidget
::SetBoneWorldToParentPoseTransform(vtkBoneWidget* bone, vtkBoneWidget* parent)
{
  double parentWorldToBonePoseRotation[4] = {1.0, 0.0, 0.0, 0.0};
  double parentWorldToBonePoseTranslation[3] = {0.0, 0.0, 0.0};
  if (parent) // For all non-root elements
    {
    parent->GetWorldToBonePoseRotation(
      parentWorldToBonePoseRotation);
    parent->GetWorldToBoneTailPoseTranslation(
      parentWorldToBonePoseTranslation);
    }

  bone->SetWorldToParentPoseRotationAndTranslation(
    parentWorldToBonePoseRotation, parentWorldToBonePoseTranslation);
}

//----------------------------------------------------------------------------
ArmatureTreeNode* vtkArmatureWidget::GetNode(vtkBoneWidget* bone)
{
  for (NodeIteratorType it = this->Bones->begin();
    it != this->Bones->end(); ++it)
    {
    if ((*it)->Bone == bone)
      {
      return *it;
      }
    }
  return NULL;
}

//----------------------------------------------------------------------------
void vtkArmatureWidget::UpdateChildren(ArmatureTreeNode* parentNode)
{
  if (!parentNode)
    {
    return;
    }

  if (parentNode->Bone->GetWidgetState() == vtkBoneWidget::Rest)
    {
    this->UpdateChildrenWidgetStateToRest(parentNode);
    }
  else if (parentNode->Bone->GetWidgetState()== vtkBoneWidget::Pose)
    {
    this->UpdateChildrenWidgetStateToPose(parentNode);
    }
}

//----------------------------------------------------------------------------
void vtkArmatureWidget
::UpdateChildrenWidgetStateToRest(ArmatureTreeNode* parentNode)
{
  if (!parentNode)
  {
  return;
  }

  for (ChildrenNodeIteratorType it = parentNode->Children.begin();
    it != parentNode->Children.end(); ++it)
    {
    this->SetBoneWorldToParentRestTransform((*it)->Bone, parentNode->Bone);

    if ((*it)->HeadLinkedToParent)
      {
      (*it)->Bone->SetWorldHeadRest(parentNode->Bone->GetWorldTailRest());
      }
    }
}

//----------------------------------------------------------------------------
void vtkArmatureWidget
::UpdateChildrenWidgetStateToPose(ArmatureTreeNode* parentNode)
{
  if (!parentNode)
  {
  return;
  }

  for (ChildrenNodeIteratorType it = parentNode->Children.begin();
    it != parentNode->Children.end(); ++it)
    {
    this->SetBoneWorldToParentPoseTransform((*it)->Bone, parentNode->Bone);
    }
}

//----------------------------------------------------------------------------
void vtkArmatureWidget::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Armature Widget " << this << "\n";
  // TO DO
}
