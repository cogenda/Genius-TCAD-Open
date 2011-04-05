#! /usr/bin/env python
# encoding: utf-8
# WARNING! All changes made to this file will be lost!

from waflib import Task, TaskGen
import os.path

@TaskGen.extension('.l')
def process(self, node):
  flags = getattr(self, 'flexflags', None)
  if flags==None:
    flags = self.env.FLEXFLAGS

  tsk = self.create_task('flex', node, node.change_ext('.yy.c'))
  tsk.flexflags = flags

class flex(Task.Task):
  color='BLUE'
  def run(self):
    src = self.inputs[0].abspath()
    tgt = self.outputs[0].abspath()
    cmd = [self.env.FLEX, '-o%s' % tgt]
    cmd.extend(self.flexflags)
    cmd.append(src)

    return self.exec_command(cmd)

def configure(conf):
	if conf.options.FLEX:
		path,exe = os.path.split(conf.options.FLEX)
		conf.find_program(exe,path_list=[path], var='FLEX')
	else:
		conf.find_program('flex',var='FLEX')

def options(opt):
	opt.add_option('--with-flex', action='store', default=None, dest='FLEX', help='flex binary to use [flex]')
