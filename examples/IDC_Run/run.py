#!/bin/env python
import os
import subprocess
import time


#--------------------------------------
# !!! user input !!!
#--------------------------------------


# model file, UNV or TIF3D format
#source_file ='aa.unv'

# track of particle
#track_file='t.dat'

# position of inject particle
#particle_file='p.dat'

# weight of particle
#weight=1.0

# background mesh resolution, cm
#resolution = 0.05

# genius install dir
GENIUS_PATH = '/opt/cogenda/current/'


def run(source_file, track_file, particle_file, weight, resolution):
  source_file = os.path.realpath(str(source_file))
  file_dir = os.path.dirname(source_file)
  file_name, file_ext = os.path.splitext(os.path.basename(source_file))

  format = ''

  if file_ext == '.unv' :
    format = 'UNVFILE'

  if file_ext == '.tif3d' :
    format = 'TIF3DFILE'

  result_file = os.path.join(file_dir, file_name+'.vtu')

  params = {
   'source_file':  source_file,
   'FORMAT'     :  format,
   'track_file' :  track_file,
   'particle_file' : particle_file,
   'weight'        : weight,
   'result_file'   : result_file,
   'resolution'    : resolution
  }



  tmpl_cmd=r'''
  GLOBAL T=300 cm=1e5 second=1e-5
  IMPORT {FORMAT:s}="{source_file:s}"

  CONTACT Property=ElectrodeRegion Type=FixedPotential

  METHOD Type=RIC NS=Basic LS=gmres pc=lu Damping=no toler.relax=1e5 maxit=25
  HOOK   id=d1 Load=particle_capture_data  string<track.data>="{track_file:s}" string<particle.data>="{particle_file:s}"  \
         real<weight>={weight:g}  real<resolution>={resolution:g}  real<fraction>=1

  SOLVE Type=steadystate

  # export result
  EXPORT   VTKFILE="{result_file:s}"
  '''

  inp_file = os.path.join(file_dir, 'run.inp')

  with open(inp_file, 'w') as fout:
      fout.write(tmpl_cmd.format(**params))


  p = subprocess.Popen( [GENIUS_PATH + '/genius/bin/genius', '-i', inp_file])
  out, err = p.communicate()
