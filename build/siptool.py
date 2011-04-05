import os
from waflib import Task, Utils, Context
from waflib.Utils import subprocess
from waflib.TaskGen import extension

@extension('.sip')
def process_sip(self, node):
  self.create_task('sip', node)

class sip(Task.Task):
  color='PINK'
  ext_out = ['.h']
  
  def run(self):
    src = self.inputs[0]
    sbf = self.inputs[0].change_ext('.sbf').get_bld()
    sbf.parent.mkdir()
    cmd = '%s -b %s -I %s -c ./ %s' % (self.env.SIP_BIN, sbf.abspath(), src.parent.abspath(), src.abspath())
    cwd = src.parent.get_bld().abspath()
    self.generator.bld.exec_command(cmd, cwd=cwd)

    sbf_data = Utils.str_to_dict(sbf.read())
    headers = Utils.to_list(sbf_data['headers'])
    sources = Utils.to_list(sbf_data['sources'])

    src_dir = src.parent.get_bld()
    self.outputs = [src_dir.find_or_declare(x) for x in sources]
    self.generator.bld.raw_deps[self.uid()] = [self.signature()] + self.outputs
    self.add_c_tasks(self.outputs)

  def add_c_tasks(self, lst):
    self.more_tasks = []
    for node in lst:
      tsk = self.generator.create_compiled_task('cxx', node)
      self.more_tasks.append(tsk)

      tsk.env.append_value('INCLUDES', [node.parent.abspath()])

      if getattr(self.generator, 'link_task', None):
        self.generator.link_task.set_run_after(tsk)
        self.generator.link_task.inputs.append(tsk.outputs[0])

  def runnable_status(self):
    ret = super(sip, self).runnable_status()
    if ret == Task.SKIP_ME:
      lst = self.generator.bld.raw_deps[self.uid()]
      if not lst[0]==self.signature():
        return Task.RUN_ME

      nodes = lst[1:]
      for x in nodes:
        try:
          os.stat(x.abspath())
        except:
          return Task.RUN_ME

      nodes = lst[1:]
      self.set_outputs(nodes)
      self.add_c_tasks(nodes)

    return ret

