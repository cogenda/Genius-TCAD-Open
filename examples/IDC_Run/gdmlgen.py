#!/bin/env python
import os
import subprocess
import time


# genius install dir
GENIUS_PATH = '/opt/cogenda/current/'


def gdmlgen(source_file):
  source_file = os.path.realpath(str(source_file))
  file_dir = os.path.dirname(source_file)
  file_name, file_ext = os.path.splitext(os.path.basename(source_file))

  format = ''

  if file_ext == '.unv' :
    format = 'UNVFILE'

  if file_ext == '.tif3d' :
    format = 'TIF3DFILE'


  params = {
   'source_file':  source_file,
   'gdml_file'  :  os.path.join(file_dir, file_name+'.gdml'),
   'FORMAT'     :  format
  }


  tmpl_gdml=r'''
  GLOBAL T=300 cm=1e5 second=1e-5
  IMPORT {FORMAT:s}="{source_file:s}"
  EXPORT GDML.SURFACE="{gdml_file:s}"
  '''

  inp_file = os.path.join(file_dir, 'gdml.inp')

  with open(inp_file, 'w') as fout:
      fout.write(tmpl_gdml.format(**params))

  p = subprocess.Popen( [GENIUS_PATH + '/genius/bin/genius', '-i', inp_file])
  out, err = p.communicate()


#--------------------------------------
# !!! user input !!!
#--------------------------------------
# source structure
#source_file ='aa.unv'
#gdmlgen(source_file)




