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

// Models includes
#include "vtkSlicerModelsLogic.h"

// Armatures includes
#include "vtkSlicerArmaturesLogic.h"

// MRML includes
#include <vtkMRMLModelNode.h>

// VTK includes
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkMatrix3x3.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLDataElement.h>
#include <vtkXMLDataParser.h>
#include <vtksys/SystemTools.hxx>

// STD includes
#include <cassert>

vtkCxxSetObjectMacro(vtkSlicerArmaturesLogic, ModelsLogic, vtkSlicerModelsLogic);

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerArmaturesLogic);

//----------------------------------------------------------------------------
vtkSlicerArmaturesLogic::vtkSlicerArmaturesLogic()
{
  this->ModelsLogic = 0;
  this->LPS2RAS = 1.;
}

//----------------------------------------------------------------------------
vtkSlicerArmaturesLogic::~vtkSlicerArmaturesLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerArmaturesLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
vtkMRMLModelNode* vtkSlicerArmaturesLogic
::AddArmatureFile(const char* fileName, bool applyPose)
{
  if (this->ModelsLogic == 0)
    {
    vtkErrorMacro( << " Models logic is NULL");
    return 0;
    }
  vtkNew<vtkXMLDataParser> armatureParser;
  armatureParser->SetFileName(fileName);
  int res = armatureParser->Parse();
  if (res == 0)
    {
    vtkErrorMacro(<<"Fail to read" << fileName << ". Not a valid armature file");
    return 0;
    }
  vtkNew<vtkPolyData> restArmature;
  this->ReadArmature(armatureParser.GetPointer(), restArmature.GetPointer(), false, !applyPose);

  vtkNew<vtkPolyData> posedArmature;
  this->ReadArmature(armatureParser.GetPointer(), posedArmature.GetPointer(), true, applyPose);

  vtkPolyData* mainArmature =
    applyPose? posedArmature.GetPointer() : restArmature.GetPointer();

  vtkPolyData* otherArmature =
    applyPose? restArmature.GetPointer() : posedArmature.GetPointer();

  vtkDataArray* otherPoints =
    applyPose? restArmature->GetPoints()->GetData() :
               posedArmature->GetPoints()->GetData();
  otherPoints->SetName(applyPose ? "RestPoints" : "PosedPoints");
  mainArmature->GetPointData()->AddArray(otherPoints);

  vtkDataArray* otherFrames =
    applyPose? restArmature->GetCellData()->GetArray("Frames") :
               posedArmature->GetCellData()->GetArray("Frames");
  otherFrames->SetName(applyPose ? "RestFrames" : "PosedFrames");
  mainArmature->GetCellData()->AddArray(otherFrames);

  if (!applyPose)
    {
    vtkNew<vtkDoubleArray> transforms;
    transforms->SetNumberOfComponents(12);
    transforms->SetName("Transforms");
    mainArmature->GetCellData()->AddArray(transforms.GetPointer());

    for (vtkIdType i = 0; i < mainArmature->GetNumberOfCells(); ++i)
      {
      vtkNew<vtkIdList> mainPts;
      mainArmature->GetCellPoints(i, mainPts.GetPointer());
      vtkNew<vtkIdList> otherPts;
      otherArmature->GetCellPoints(i, otherPts.GetPointer());
      double mainHead[3];
      double mainTail[3];
      mainArmature->GetPoint(mainPts->GetId(0), mainHead);
      mainArmature->GetPoint(mainPts->GetId(1), mainTail);
      double otherHead[3];
      double otherTail[3];
      otherArmature->GetPoint(otherPts->GetId(0), otherHead);
      otherArmature->GetPoint(otherPts->GetId(1), otherTail);
      double translation[3];
      vtkMath::Subtract(otherHead, mainHead, translation);
      vtkMath::Subtract(mainTail, mainHead, mainTail);
      vtkMath::Subtract(otherTail, otherHead, otherTail);
      double rotation[3][3] = {{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
      this->ComputeTransform(mainTail, otherTail, rotation);
      double transform[4][3];
      memcpy(transform, rotation, 3*3*sizeof(double));
      memcpy(transform[3], translation, 3*sizeof(double));
      transforms->InsertNextTuple(transform[0]);

      //vtkMath::Add(mainHead, transform[3], mainHead);
      //vtkNew<vtkMatrix3x3> matrix;
      //matrix->DeepCopy(transform[0]);
      //matrix->Transpose();
      //matrix->MultiplyPoint(mainTail, mainTail);
      //vtkMath::Add(mainTail, mainHead, mainTail);
      //mainArmature->GetPoints()->SetPoint(mainPts->GetId(0), mainHead);
      //mainArmature->GetPoints()->SetPoint(mainPts->GetId(1), mainTail);
      }
    }
  vtkMRMLModelNode* modelNode= this->ModelsLogic->AddModel(mainArmature);
  std::string modelName = vtksys::SystemTools::GetFilenameName(fileName);
  modelNode->SetName(modelName.c_str());
  return modelNode;
}

//----------------------------------------------------------------------------
void vtkSlicerArmaturesLogic::ReadArmature(vtkXMLDataParser* armatureParser,
  vtkPolyData* armature,
  bool applyPose, bool applyVertexGroup)
{
  vtkNew<vtkPoints> points;
  points->SetDataTypeToDouble();
  armature->SetPoints(points.GetPointer());
  armature->Allocate(100);

  vtkNew<vtkUnsignedCharArray> colors;
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");
  armature->GetPointData()->SetScalars(colors.GetPointer());

  vtkNew<vtkIdTypeArray> associatedCells;
  associatedCells->SetNumberOfComponents(2);
  associatedCells->SetName("SurfaceCells");
  armature->GetCellData()->AddArray(associatedCells.GetPointer());

  vtkNew<vtkDoubleArray> frames;
  frames->SetNumberOfComponents(9);
  frames->SetName("Frames");
  armature->GetCellData()->AddArray(frames.GetPointer());
  if (applyPose)
    {
    vtkNew<vtkDoubleArray> transforms;
    transforms->SetNumberOfComponents(9);
    transforms->SetName("Transforms");
    armature->GetCellData()->AddArray(transforms.GetPointer());
    }

  vtkXMLDataElement* armatureElement = armatureParser->GetRootElement();

  double origin[3] = {0., 0., 0.};
  armatureElement->GetVectorAttribute("location", 3, origin);

  double scale[3] = {1.0, 1.0, 1.0};
  armatureElement->GetVectorAttribute("scale", 3, scale);
  double scaleMat[3][3];
  vtkMath::Identity3x3(scaleMat);
  scaleMat[0][0] = scale[0];
  scaleMat[1][1] = scale[1];
  scaleMat[2][2] = scale[2];

  double orientationXYZW[4] = {0.0, 0.0, 0.0, 1.0};
  armatureElement->GetVectorAttribute("orientation", 4, orientationXYZW);
  double orientationWXYZ[4] = {orientationXYZW[3], orientationXYZW[0],
                               orientationXYZW[1], orientationXYZW[2]};
  double mat[3][3];
  vtkMath::QuaternionToMatrix3x3(orientationWXYZ, mat);

  vtkMath::Multiply3x3(mat, scaleMat, mat);

  for (int child = 0; child < armatureElement->GetNumberOfNestedElements(); ++child)
    {
    this->ReadBone(armatureElement->GetNestedElement(child), armature,
                   applyPose, applyVertexGroup,
                   origin, mat);
    }
}

//----------------------------------------------------------------------------
void vtkSlicerArmaturesLogic::ReadBone(vtkXMLDataElement* boneElement,
                                       vtkPolyData* polyData,
                                       bool applyPose, bool applyVertexGroup,
                                       const double origin[3],
                                       const double parentMatrix[3][3],
                                       const double parentLength)
{
  double parentTransMatrix[3][3];
  vtkMath::Transpose3x3(parentMatrix, parentTransMatrix);

  double localOrientationXYZW[4] = {0., 0., 0., 1.};
  boneElement->GetVectorAttribute("orientation", 4, localOrientationXYZW);
  double localOrientationWXYZ[4] = {localOrientationXYZW[3], localOrientationXYZW[0],
                                    localOrientationXYZW[1], localOrientationXYZW[2]};
  double mat[3][3];
  vtkMath::QuaternionToMatrix3x3(localOrientationWXYZ, mat);
  vtkMath::Invert3x3(mat, mat);
  vtkMath::Multiply3x3(mat, parentMatrix, mat);

  double restTail[3] = {mat[1][0], mat[1][1], mat[1][2]};

  if (applyPose)
    {
    vtkXMLDataElement * poseElement =
      boneElement->FindNestedElementWithName("pose");
    double poseRotationXYZW[4] = {0., 0., 0., 0.};
    poseElement->GetVectorAttribute("rotation", 4, poseRotationXYZW);
    double poseRotationWXYZ[4] = {poseRotationXYZW[3], poseRotationXYZW[0],
                                  poseRotationXYZW[1], poseRotationXYZW[2]};
    double poseMat[3][3];
    vtkMath::QuaternionToMatrix3x3(poseRotationWXYZ, poseMat);
    vtkMath::Invert3x3(poseMat, poseMat);
    vtkMath::Multiply3x3(poseMat, mat, mat);
    }

  double localHead[3] = {0., 0., 0.};
  boneElement->GetVectorAttribute("head", 3, localHead);

  double head[3] = {0., 0., 0.};
  head[0] = localHead[0];
  head[1] = localHead[1] + parentLength;
  head[2] = localHead[2];
  vtkMath::Multiply3x3(parentTransMatrix, head, head);
  vtkMath::Add(origin, head, head);

  double localTail[3] = {0., 0. ,0.};
  boneElement->GetVectorAttribute("tail", 3, localTail);

  double tail[3] = {0., 0. ,0};
  vtkMath::Subtract(localTail, localHead, tail);
  double length = vtkMath::Norm(tail);

  tail[0] = mat[1][0]; tail[1] = mat[1][1]; tail[2] = mat[1][2];
  vtkMath::MultiplyScalar(tail, length);
  vtkMath::MultiplyScalar(restTail, length);

  vtkDataArray* frames = polyData->GetCellData()->GetArray("Frames");
  frames->InsertNextTuple9(
    mat[0][0] * length, mat[0][1] * length, mat[0][2] * length,
    mat[1][0] * length, mat[1][1] * length, mat[1][2] * length,
    mat[2][0] * length, mat[2][1] * length, mat[2][2] * length);

  if (applyPose)
    {
    double transform[3][3] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
    this->ComputeTransform(restTail, tail, transform);
    polyData->GetCellData()->GetArray("Transforms")->InsertNextTuple9(
      transform[0][0], transform[0][1], transform[0][2],
      transform[1][0], transform[1][1], transform[1][2],
      transform[2][0], transform[2][1], transform[2][2]);
    }

  vtkMath::Add(head, tail, tail);

  vtkIdType indexes[2];
  indexes[0] = polyData->GetPoints()->InsertNextPoint(
    head[0] * this->LPS2RAS, head[1] * this->LPS2RAS, head[2]);
  indexes[1] = polyData->GetPoints()->InsertNextPoint(
    tail[0] * this->LPS2RAS, tail[1] * this->LPS2RAS, tail[2]);
  static unsigned char color[3] = {255, 255, 255};
  unsigned char offset = 20;
  polyData->GetPointData()->GetScalars("Colors")->InsertNextTuple3(
    color[0], color[1], color[2]);
  color[0] -= offset; color[1] -= offset; color[2] -= offset;
  polyData->GetPointData()->GetScalars("Colors")->InsertNextTuple3(
    color[0], color[1], color[2]);
  color[0] -= offset; color[1] -= offset; color[2] -= offset;

  int cellId = polyData->InsertNextCell(VTK_LINE, 2, indexes);

  polyData->GetCellData()->GetArray("SurfaceCells")->InsertNextTuple2(0., 0.);
  vtkIdType boneCellsId =
    polyData->GetCellData()->GetArray("SurfaceCells")->GetNumberOfTuples();
  --boneCellsId;
  vtkIdType surfaceId = 0;

  for (int child = 0; child < boneElement->GetNumberOfNestedElements(); ++child)
    {
    vtkXMLDataElement* childElement = boneElement->GetNestedElement(child);
    if (strcmp(childElement->GetName(),"bone") == 0)
      {
      this->ReadBone(childElement, polyData, applyPose, applyVertexGroup, head, mat, length);
      }
    else if (strcmp(childElement->GetName(),"pose") == 0)
      {
      // already parsed
      }
    else if (strcmp(childElement->GetName(),"vertex_group") == 0)
      {
      if (applyVertexGroup)
        {
        vtkNew<vtkIdTypeArray> associatedCells;
        associatedCells->SetNumberOfComponents(1);
        std::ostringstream associatedCellsName;
        associatedCellsName << "SurfaceCells" << cellId << '-' << surfaceId;
        associatedCells->SetName(associatedCellsName.str().c_str());

        std::stringstream vertices;
        vertices << childElement->GetCharacterData();
        while (!vertices.eof())
          {
          vtkIdType pointId;
          vertices >> pointId;
          associatedCells->InsertNextTuple1(pointId);
          }
        vtkIdType surfaceCellsArrayId =
          polyData->GetCellData()->AddArray(associatedCells.GetPointer());

        vtkIdTypeArray::SafeDownCast(polyData->GetCellData()->GetArray("SurfaceCells"))
          ->SetComponent(boneCellsId, surfaceId++, surfaceCellsArrayId);
        }
      }
    else
      {
      vtkWarningMacro( << "XML element " << childElement->GetName()
                       << "is not supported");
      }
    }
}

//----------------------------------------------------------------------------
void vtkSlicerArmaturesLogic::ReadPose(vtkXMLDataElement* poseElement,
                                       double parentOrientation[4], bool preMult)
{
  double poseRotationXYZW[4] = {0., 0., 0., 0.};
  poseElement->GetVectorAttribute("rotation", 4, poseRotationXYZW);
  double poseRotationWXYZ[4] = {poseRotationXYZW[3], poseRotationXYZW[0],
                                poseRotationXYZW[1], poseRotationXYZW[2]};
  if (preMult)
    {
    vtkMath::MultiplyQuaternion(poseRotationWXYZ, parentOrientation, parentOrientation);
    }
  else
    {
    vtkMath::MultiplyQuaternion(parentOrientation, poseRotationWXYZ, parentOrientation);
    }
}

//----------------------------------------------------------------------------
void vtkSlicerArmaturesLogic::ComputeTransform(double start[3], double end[3], double mat[3][3])
{
  double startNormalized[3] = {start[0], start[1], start[2]};
  double startNorm = vtkMath::Normalize(startNormalized);
  double endNormalized[3] = {end[0], end[1], end[2]};
  double endNorm = vtkMath::Normalize(endNormalized);

  double rotationAxis[3] = {0., 0., 0.};
  vtkMath::Cross(startNormalized, endNormalized, rotationAxis);
  vtkMath::Normalize(rotationAxis);
  if (rotationAxis[0] != 0. || rotationAxis[1] != 0. || rotationAxis[2] != 0.)
    {
    double angle = vtkSlicerArmaturesLogic::ComputeAngle(startNormalized, endNormalized);
    vtkSlicerArmaturesLogic::ComputeAxisAngleMatrix(rotationAxis, angle, mat);
    }
  else
    {
    vtkMath::Identity3x3(mat);
    }

  double scaleMatrix[3][3] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};
  scaleMatrix[0][0] = scaleMatrix[1][1] = scaleMatrix[2][2] = endNorm / startNorm;
  vtkMath::Multiply3x3(mat, scaleMatrix, mat);
  //vtkMath::Multiply3x3(scaleMatrix, transform, transform);
}

//-----------------------------------------------------------------------------
double vtkSlicerArmaturesLogic::ComputeAngle(double v1[3], double v2[3])
{
  double dot = vtkMath::Dot(v1, v2);
  double angle = acos(dot);
  return angle;
}

//-----------------------------------------------------------------------------
void vtkSlicerArmaturesLogic::ComputeAxisAngleMatrix(double axis[3], double angle, double mat[3][3])
{
  /* rotation of angle radials around axis */
  double vx, vx2, vy, vy2, vz, vz2, co, si;

  vx = axis[0];
  vy = axis[1];
  vz = axis[2];
  vx2 = vx * vx;
  vy2 = vy * vy;
  vz2 = vz * vz;
  co = cos(angle);
  si = sin(angle);

  mat[0][0] = vx2 + co * (1.0f - vx2);
  mat[0][1] = vx * vy * (1.0f - co) + vz * si;
  mat[0][2] = vz * vx * (1.0f - co) - vy * si;
  mat[1][0] = vx * vy * (1.0f - co) - vz * si;
  mat[1][1] = vy2 + co * (1.0f - vy2);
  mat[1][2] = vy * vz * (1.0f - co) + vx * si;
  mat[2][0] = vz * vx * (1.0f - co) + vy * si;
  mat[2][1] = vy * vz * (1.0f - co) - vx * si;
  mat[2][2] = vz2 + co * (1.0f - vz2);
}
