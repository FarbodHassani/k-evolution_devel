# state file generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = GetRenderView()
renderView1.ViewSize=[1024,1024]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [599995.65625, 599999.046875, 599999.640625]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [1056871.9570939988, 1205136.5295483228, 1702873.4923122516]
renderView1.CameraFocalPoint = [599995.65625, 599999.046875, 599999.640625]
renderView1.CameraViewUp = [0.15166002726454209, 0.8387167447855731, -0.5230233820264734]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 346397.9469475421
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1


# create a new 'LightCone Series Reader'
reader = LightConeSeriesReader(registrationName='lcdm_Mn0d13_gevolution_snap000_cdm_redist.0', FileNames=['/scratch/snx3000/farbodh/Gadget2_snapshot/output_redist/lcdm_Mn0d13_gevolution_snap000_cdm_redist.0'])
reader.PointArrayStatus = ['velocity']

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=reader)
clip1.ClipType = 'Sphere'
clip1.Scalars = ['POINTS', '']
clip1.ClipType.Center = [600000.0, 600000.0, 600000.0]
clip1.ClipType.Radius = 600000.0


# show data from clip1
clip1Display = Show(clip1)

# get color transfer function/color map for 'velocity'
velocityLUT = GetColorTransferFunction('velocity')
velocityLUT.RGBPoints = [0.06293770433971857, 0.231373, 0.298039, 0.752941, 3093.7891192153024, 0.865003, 0.865003, 0.865003, 6187.515300726265, 0.705882, 0.0156863, 0.14902]
velocityLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'velocity'
velocityPWF = GetOpacityTransferFunction('velocity')
velocityPWF.Points = [0.06293770433971857, 0.0, 0.5, 0.0, 6187.515300726265, 1.0, 0.5, 0.0]
velocityPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
clip1Display.Representation = 'Points'
clip1Display.ColorArrayName = ['POINTS', 'velocity']
clip1Display.LookupTable = velocityLUT


# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [-3088.29833984375, 0.0, 0.5, 0.0, 1970.0611572265625, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [-3088.29833984375, 0.0, 0.5, 0.0, 1970.0611572265625, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for velocityLUT in view renderView1
velocityLUTColorBar = GetScalarBar(velocityLUT, renderView1)
velocityLUTColorBar.Title = 'velocity'
velocityLUTColorBar.ComponentTitle = 'Magnitude'

# set color bar visibility
velocityLUTColorBar.Visibility = 1

# show color legend
clip1Display.SetScalarBarVisibility(renderView1, True)


SetActiveSource(clip1)
ResetCamera()
SaveScreenshot("/scratch/snx3000/jfavre/clip.png")
