/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Julien Finet, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// Armatures includes
#include <vtkSlicerArmaturesLogic.h>

// VTK includes
#include <vtkMath.h>
#include <vtkMatrix3x3.h>
#include <vtkNew.h>

bool TestComputeTransform();
bool TestComputeTransform(double startX, double startY, double startZ, double endX, double endY, double endZ);

//----------------------------------------------------------------------------
int vtkSlicerArmaturesLogicTest4(int argc, char* argv[])
{
  bool res = true;
  res = TestComputeTransform() && res;
  return res ? EXIT_SUCCESS : EXIT_FAILURE;
}

//----------------------------------------------------------------------------
bool TestComputeTransform()
{
  bool res = true;
  double scales[3] = {1., 10., 0.10};
  for (int i = 0; i < 3; ++i)
    {
    double s = scales[i];
    // X
    res = TestComputeTransform(1.,0.,0., s,0.,0.) && res;
    res = TestComputeTransform(1.,0.,0., -s,0.,0.) && res;
    res = TestComputeTransform(1.,0.,0., 0.,s,0.) && res;
    res = TestComputeTransform(1.,0.,0., 0.,-s,0.) && res;
    res = TestComputeTransform(1.,0.,0., 0.,0.,s) && res;
    res = TestComputeTransform(1.,0.,0., 0.,0.,-s) && res;

    res = TestComputeTransform(-1.,0.,0., s,0.,0.) && res;
    res = TestComputeTransform(-1.,0.,0., -s,0.,0.) && res;
    res = TestComputeTransform(-1.,0.,0., 0.,s,0.) && res;
    res = TestComputeTransform(-1.,0.,0., 0.,-s,0.) && res;
    res = TestComputeTransform(-1.,0.,0., 0.,0.,s) && res;
    res = TestComputeTransform(-1.,0.,0., 0.,0.,-s) && res;

    // Y
    res = TestComputeTransform(0.,1.,0., s,0.,0.) && res;
    res = TestComputeTransform(0.,1.,0., -s,0.,0.) && res;
    res = TestComputeTransform(0.,1.,0., 0.,s,0.) && res;
    res = TestComputeTransform(0.,1.,0., 0.,-s,0.) && res;
    res = TestComputeTransform(0.,1.,0., 0.,0.,s) && res;
    res = TestComputeTransform(0.,1.,0., 0.,0.,-s) && res;

    res = TestComputeTransform(0.,-1.,0., s,0.,0.) && res;
    res = TestComputeTransform(0.,-1.,0., -s,0.,0.) && res;
    res = TestComputeTransform(0.,-1.,0., 0.,s,0.) && res;
    res = TestComputeTransform(0.,-1.,0., 0.,-s,0.) && res;
    res = TestComputeTransform(0.,-1.,0., 0.,0.,s) && res;
    res = TestComputeTransform(0.,-1.,0., 0.,0.,-s) && res;

    // Z
    res = TestComputeTransform(0.,0.,1., s,0.,0.) && res;
    res = TestComputeTransform(0.,0.,1., -s,0.,0.) && res;
    res = TestComputeTransform(0.,0.,1., 0.,s,0.) && res;
    res = TestComputeTransform(0.,0.,1., 0.,-s,0.) && res;
    res = TestComputeTransform(0.,0.,1., 0.,0.,s) && res;
    res = TestComputeTransform(0.,0.,1., 0.,0.,-s) && res;

    res = TestComputeTransform(0.,0.,-1., s,0.,0.) && res;
    res = TestComputeTransform(0.,0.,-1., -s,0.,0.) && res;
    res = TestComputeTransform(0.,0.,-1., 0.,s,0.) && res;
    res = TestComputeTransform(0.,0.,-1., 0.,-s,0.) && res;
    res = TestComputeTransform(0.,0.,-1., 0.,0.,s) && res;
    res = TestComputeTransform(0.,0.,-1., 0.,0.,-s) && res;
    }
  return res;
}

//----------------------------------------------------------------------------
bool TestComputeTransform(double startX, double startY, double startZ, double endX, double endY, double endZ)
{
  double start[3] = {startX, startY, startZ};
  double end[3] = {endX, endY, endZ};
  double rotation[3][3] = {{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
  vtkSlicerArmaturesLogic::ComputeTransform(start, end, rotation);

  vtkNew<vtkMatrix3x3> matrix;
  matrix->DeepCopy(rotation[0]);
  matrix->Transpose();
  matrix->MultiplyPoint(start, start);
  double diff[3];
  vtkMath::Subtract(start, end, diff);
  return vtkMath::Norm(diff) < 0.000001;
}

