#!/usr/bin/env python
# coding: utf-8
# ParaView LightCone Reader Minimal Test
# Created by Jean M. Favre, August 16, 2021
# tested with ParaView-EGL v 5.9
# In[2]:


from paraview.simple import *
from paraview.selection import *


# Using a terminal on this node, we request a new allocation:
# 
# salloc -N 64 -n 64-A csstaff -C gpu --time=00:15:00
# 
# salloc: Granted job allocation 30807331 salloc: Waiting for resource configuration salloc: Nodes nid0[3508-3511] are ready for job
# 
# We *first* run pvserver on the newly allocated nodes
# 
# srun pvserver
# 
# The jupyter client must connect to the first MPI task on the new allocation, i.e. connect to the first node

# In[3]:


Connect("nid01980")


# In[4]:


LoadPlugin("/users/jfavre/Projects/Adamek/ParaViewLightConePlugin/build59/lib64/paraview-5.9/plugins/pvLightConeReader/pvLightConeReader.so", ns=globals())


# In[5]:


renderView1 = GetRenderView()

reader = LightConeSeriesReader(FileNames=['/scratch/snx3000/farbodh/Gadget2_snapshot/output_redist/lcdm_Mn0d13_gevolution_snap000_cdm_redist.0'])
reader.PointArrayStatus = ['velocity']
reader.UpdatePipelineInformation()


# In[6]:


reader.DistributedSnapshot = 0 # to read a single file (you should be able to run on the single node already allocated)
reader.DistributedSnapshot = 1 # to read all files (you must do a second node allocation


# In[7]:


reader.UpdatePipeline()


# In[8]:


from ipyparaview.widgets import PVDisplay
disp = PVDisplay(GetActiveView())
w = display(disp)


# In[9]:


clip1 = Clip(Input=reader)
clip1.ClipType = 'Sphere'
clip1.ClipType.Center = [600000.0, 600000.0, 600000.0]
clip1.ClipType.Radius = 200000.0


# In[10]:


repr2 = Show(clip1)

repr2.Representation = 'Points'
ColorBy(repr2, ('POINTS', 'velocity', 'Magnitude'))
ResetCamera()


# In[11]:


Disconnect()


# In[ ]:




