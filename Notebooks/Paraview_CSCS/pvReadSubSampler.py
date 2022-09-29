# state file generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = GetRenderView()
renderView1.Background = [0.3686274509803922]*3

fname = '/scratch/snx3000/farbodh/Lightcone_data/LCDM_lightcones/lcdm_lightcone1_0104_cdm'

reader = LightConeSeriesReader(registrationName='lcdm_lightcone1', FileNames=[fname])
reader.PointArrayStatus = ['velocity', 'id']
reader.DistributedSnapshot = 0
reader.UpdatePipeline()

programmableFilter1 = ProgrammableFilter(registrationName='ProgrammableFilter1', Input=reader)
programmableFilter1.OutputDataSetType = 'Same as Input'
programmableFilter1.Script = """
from vtk.vtkCommonDataModel import vtkPolyData
from vtk.vtkCommonCore import vtkIdList
from vtk.numpy_interface import dataset_adapter as dsa
from vtk import VTK_POLY_VERTEX
import numpy as np

Ratio = 1000

def SubSampler(input, output):
  my_input = dsa.WrapDataObject(input)
  nbpts = input.GetNumberOfPoints()
  indices = np.arange(0, nbpts, Ratio)
  nnodes = indices.size
  print("found ", nbpts, ". Extracting ", nnodes)
  ptIds = vtkIdList()
  ptIds.SetNumberOfIds(nnodes)

  my_output = dsa.WrapDataObject(output)
  for a in range(nnodes):
    ptIds.SetId(a , a)

  output.Allocate(1)
  output.InsertNextCell(VTK_POLY_VERTEX , ptIds)
  my_output.Points = my_input.Points[indices]
  my_output.PointData.append(my_input.PointData["velocity"][indices], "velocity")
  my_output.PointData.append(my_input.PointData["id"][indices], "id")
  
if inputs[0].IsA("vtkMultiBlockDataSet"):
    output.CopyStructure(inputs[0].VTKObject)
    iter = inputs[0].NewIterator()
    iter.UnRegister(None)
    iter.InitTraversal()
    while not iter.IsDoneWithTraversal():
        curInput = iter.GetCurrentDataObject()
        curOutput = curInput.NewInstance()
        curOutput.UnRegister(None)
        output.SetDataSet(iter, curOutput)
        SubSampler(curInput, curOutput)
        iter.GoToNextItem();
else:
  SubSampler(input, output)"""

rep1 = Show(programmableFilter1, renderView1, 'GeometryRepresentation')
rep1.Representation = 'Points'
ColorBy(rep1, ('POINTS', 'velocity', 'Magnitude'))

myselect2 = QuerySelect(QueryString="mag(velocity) >= 1000.", Source=reader)
sel2     = ExtractSelection(Input=reader)

rep2 = Show(sel2, renderView1, 'GeometryRepresentation')
rep2.Representation = 'Points'
ColorBy(rep2, ('POINTS', 'velocity', 'Magnitude'))

ResetCamera()
# use the current default object just defined above
SaveData("/scratch/snx3000/jfavre/clip3.csv")

