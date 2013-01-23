APPNAME='genius'

top='.'
out='build/default'

import sys
import os
import re

import waflib.Configure
from waflib.Task import Task
from waflib import Utils
import tempfile, string


waflib.Configure.autoconfig = True

def config_guess():
  import platform
  guess={}
  guess['machine']  = platform.machine()
  guess['platform'] = platform.system()   #FIXME for jython, has to dig further about the underlying OS
  return guess

def options(opt):
  opt.load('compiler_c')
  opt.load('compiler_cxx')
  opt.load('compiler_fc')
  opt.load('myflex mybison', tooldir='./build')

  guess = config_guess()
  if guess['platform']=='Windows':
    opt.add_option('--target', action='store', default=None, dest='target_arch', help='target system architecture (ia32 intel64) [default: native]')
    opt.add_option('--crt', action='store', default='MT', dest='crt_version', help='Win32 only: CRT library version (MT MD) [default: MT]')

  opt.add_option('--version', action='store', default=None, dest='version_str', help='Version string for the build')
  opt.add_option('--with-git', action='store', default=None, dest='GIT', help='git binary [git]')
  opt.add_option('--cc-opt', action='store', default=None, dest='cc_opt', help='CC optimization options. [default: autodetect]')
  opt.add_option('--debug', action='store_true', default=False, dest='debug', help='Enable debug')
  opt.add_option('--with-netgen-dir', action='store', default=None, dest='netgen_dir', help='Directory to Netgen.')
  opt.add_option('--with-cgns-dir', action='store', default=None, dest='cgns_dir', help='Directory to CGNS.')
  opt.add_option('--with-vtk-dir', action='store', default=None, dest='vtk_dir', help='Directory to VTK.')
  opt.add_option('--with-vtk-ver', action='store', default='vtk-5.4', dest='vtk_ver', help='Version of VTK [vtk-5.4]')
  opt.add_option('--with-petsc-dir',  action='store', default='/usr/local/petsc', dest='petsc_dir', help='Directory to Petsc.')
  opt.add_option('--with-petsc-arch', action='store', default='linux-intel-cc', dest='petsc_arch', help='Petsc Arch.')
  opt.add_option('--with-ams', action='store_true', default=False, dest='ams_enabled', help='Build with AMS')
  opt.add_option('--with-ams-dir',  action='store', default='/usr/local/ams', dest='ams_dir', help='Directory to AMS.')
  opt.add_option('--with-slepc', action='store_true', default=False, dest='slepc_enabled', help='Build with Slepc')
  opt.add_option('--with-slepc-dir',  action='store', default='/usr/local/slepc', dest='slepc_dir', help='Directory to Slepc.')

def configure(conf):
  guess = config_guess()
  platform = guess['platform']

  if platform=='Windows':
    target = conf.options.target_arch
    if target==None:
      if   guess['machine']=='AMD64': target='intel64'
      elif guess['machine']=='x86':   target='ia32'
    conf.env.TARGET_ARCH = target
    conf.env['MSVC_VERSIONS'] = ['intel 12', 'intel 11']
    conf.env['MSVC_TARGETS']  = [target]

  if not platform=='Windows':
    conf.options.check_c_compiler = 'icc ' + conf.options.check_c_compiler
    conf.options.check_cxx_compiler = 'icpc ' + conf.options.check_cxx_compiler
    conf.options.check_fc = 'ifort ' + conf.options.check_fc

  conf.load('compiler_c compiler_cxx compiler_fc')

  conf.env.PLATFORM = platform # Windows Linux Darwin ...
  if   platform=='Linux':    conf.define('LINUX', 1); conf.define('DLLHOOK',1)
  elif platform=='Windows':  conf.define('WINDOWS', 1); conf.define('WIN32',1)
  elif platform=='Darwin':   conf.define('DARWIN', 1)
  elif platform=='AIX':      conf.define('AIX', 1)

  # try to locate git
  if conf.options.GIT:
    path,exe = os.path.split(conf.options.GIT)
    conf.find_program(exe,path_list=[path],var='GIT')
  else:
    conf.find_program('git',var='GIT', mandatory=False)

  # version
  version_str = conf.options.version_str
  if version_str==None:
    if conf.env.GIT:
      try:
        out = conf.cmd_and_log([conf.env.GIT, 'rev-parse', '--short', 'HEAD'],
                               output=waflib.Context.STDOUT)
        version_str = 'git-%s' % string.strip(out)
      except: pass
  if version_str==None:
    import datetime
    version_str = datetime.date.today().strftime('%Y%m%d')
  conf.msg('Setting version string', version_str)
  conf.env.append_value('DEFINES_VERSION', 'PACKAGE_VERSION="%s"' % version_str)

  # general C flags
  ccflags_common, cxxflags_common, fcflags_common = ([],[],[])
  ccflags_optimize, cxxflags_optimize, fcflags_optimize=([],[],[])
  ldflags_common=[]
  ldflags_shlib=[]
  if platform=='Linux':
    if conf.env['COMPILER_CC'] in ['gcc', 'icc']:
      ccflags_common.extend(['-fPIC'])
      cxxflags_common.extend(['-fPIC'])
      fcflags_common.extend(['-fPIC'])

    ldflags_common.extend(['-ldl', '-Wl,--export-dynamic'])
    if conf.env['COMPILER_CC'] in ['icc']:
      ldflags_common.extend(['-static-intel'])

    conf.env['cxxshlib_PATTERN'] = '%s.so'

  elif platform=='Windows':
    if conf.env['COMPILER_CC'] in ['msvc', 'icc']:
      cxxflags_common.extend(['/EHsc'])
      if conf.options.crt_version=='MT':
        ccflags_common.extend(conf.env.CFLAGS_CRT_MULTITHREADED)
        cxxflags_common.extend(conf.env.CFLAGS_CRT_MULTITHREADED)
      else:
        ccflags_common.extend(conf.env.CFLAGS_CRT_MULTITHREADED_DLL)
        cxxflags_common.extend(conf.env.CFLAGS_CRT_MULTITHREADED_DLL)
        ldflags_common.append('/NODEFAULTLIB:LIBCMT')

  elif platform=='AIX':
    # AIX only support 64bit
    if conf.env['COMPILER_CC'] in ['xlc']:
      ccflags_common.extend(['-q64','-qpic'])
      cxxflags_common.extend(['-q64','-qpic','-qrtti=all'])

    if conf.env['COMPILER_FC'] in ['xlf90']:
      fcflags_common.extend(['-q64','-qpic'])

    ldflags_common.extend(['-q64'])
    ldflags_common.extend(['-bexpall'])
    conf.env['cxxshlib_PATTERN'] = '%s.so'


  conf.env.append_value('CFLAGS', ccflags_common)
  conf.env.append_value('CXXFLAGS', cxxflags_common)
  conf.env.append_value('FCFLAGS', fcflags_common)

  conf.env.append_value('LINKFLAGS', ldflags_common)
  conf.env.append_value('LINKFLAGS_cshlib', ldflags_shlib)
  conf.env.append_value('LINKFLAGS_cxxshlib', ldflags_shlib)


  if platform=='Windows':
    conf.check_libs_msvc('RPCRT4 User32')

  def _split_win_paths(p):
    return p.split(';')
  if platform=='Windows':
    for i in ['LIBPATH']:
      plist = conf.env[i]
      newlist = []
      for p in plist: newlist.extend(_split_win_paths(p))
      conf.env[i] = newlist

  # {{{ test_opt()
  def test_opt(opt, lang='c'):
    tmpdir = tempfile.mkdtemp()

    if lang=='c':
      fname='test.c'
      cc=conf.env['CC']
    elif lang=='cxx':
      fname='test.cc'
      cc=conf.env['CXX']

    f=open(os.path.join(tmpdir,fname), 'w')
    f.write('int main(void){return 0;}')
    f.close()

    cmd = []
    if isinstance(cc, str): cmd.append(cc)
    else:                   cmd.extend(cc)

    if isinstance(opt, list): cmd.extend(opt)
    else:                     cmd.extend(str(opt).split())
    cmd.append(fname)

    try:
      out, err = conf.cmd_and_log(cmd, cwd=tmpdir, output=waflib.Context.BOTH)

      if conf.env['COMPILER_CC'] in ['icc', 'msvc']:
        m = re.search('command line warning', err)
        return m==None
      if conf.env['COMPILER_CC'] in ['gcc']:
        m = re.search('unrecognized command line option', err)
        return m==None
      return True
    except Exception, e:
      return False
  # }}}

  # {{{ simplifyLinkerOptions()
  def simplifyLinkerOptions(opts):
    # append -Bdynamic at the end to reset to default behaviour
    opts = list(opts)
    opts.append('-Wl,-Bdynamic')

    # {{{ 1: -Bstatic/dynamic immediately following each other
    dups = ['-Wl,-Bstatic', '-Wl,-Bdynamic']
    nopts = [ opts[0] ]
    for o in opts[1:]:
      if o in dups and nopts[-1] in dups:
        nopts.pop()
      nopts.append(o)
    opts = nopts
    # }}}

    # {{{ 2: --start-group immediately followed by --end-group
    nopts =[ opts[0] ]
    for o in opts[1:]:
      if o=='-Wl,--end-group' and nopts[-1]=='-Wl,--start-group':
        nopts.pop()
      else:
        nopts.append(o)
    opts = nopts
    # }}}

    return opts
  # }}}

  # {{{ _fromCygpath()
  def _fromCygpath(cygpath):
    if not cygpath.startswith('/cygdrive'):
      return cygpath
    parts=cygpath.split('/')
    newparts=['%s:' % parts[2], '\\']
    newparts.extend(parts[3:])
    return os.path.join(*newparts)
  # }}}

  # {{{ parse_lib_str
  def parse_lib_str(str):
    toks = re.split('(?<!\\\\)\s+', str)

    flags = []
    for tok in toks:
      # {{{ -l<lib>
      res = re.match('-l(.*)', tok)
      if res:
        flags.append(conf.env.LIB_ST % res.group(1))
        # TODO: on windows, stlib and shlib has
        # different name convension!
        continue
      # }}}

      # {{{ -L<LIBPATH>
      res = re.match('-L(.*)', tok)
      if res:
        if platform=='Windows' :
          lib_dir = _fromCygpath(res.group(1))
        else:
          lib_dir = res.group(1)
        flags.append(conf.env.LIBPATH_ST % lib_dir)
        continue
      # }}}

      # {{{ xyz.lib  or xyz.a
      res = re.match('.*\.(?:lib|a)$', tok)
      if res:
        if platform=='Windows':
          f = _fromCygpath(res.group(0))
        else:
          f = res.group(0)
        flags.append(f)
        continue
      # }}}

      # {{{ options
      allow_lnk_opts = ['-Bstatic', '-Bdynamic', '--start-group', '--end-group']
      if tok in allow_lnk_opts:
        flags.append('-Wl,%s' % tok)
        continue
      if tok in ['-Wl,%s'%x for x in allow_lnk_opts]:
        flags.append(tok)
        continue
      # }}}

    return flags
  # }}}

  # {{{ parse_inc_str
  def parse_inc_str(str):
    toks = re.split('(?<!\\\\)\s+', str)
    inc_dirs=[]
    for tok in toks:
      res = re.match('-I(.*)', tok)
      if res:
        if platform=='Windows':
          inc_dirs.append(_fromCygpath(res.group(1)))
        else:
          inc_dirs.append(res.group(1))

    return inc_dirs
  # }}}

  # {{{ test_imp_lib()
  def test_imp_lib(lang='c'):
    tmpdir = tempfile.mkdtemp()

    if lang=='c':
      fname='test.c'
      cc=conf.env['CC']
      code = 'int main(void){return 0;}'
    elif lang=='cxx':
      fname='test.cc'
      cc=conf.env['CXX']
      code = 'int main(void){return 0;}'
    elif lang=='f':
      fname='test.f90'
      cc=conf.env['FC']
      code = '''
PROGRAM  TEST
 IMPLICIT  NONE
 INTEGER            :: n
END PROGRAM  TEST
'''
    else:
      raise ValueError

    f=open(os.path.join(tmpdir,fname), 'w')
    f.write(code)
    f.close()

    cmd = []
    if isinstance(cc, str): cmd.append(cc)
    else:                   cmd.extend(cc)
    cmd.extend(['-v', fname])
    out, err = conf.cmd_and_log(cmd, cwd=tmpdir, output=waflib.Context.BOTH)

    lines = err.split(os.linesep)
    for line in lines:
      cmd = line.split()
      if cmd==None or len(cmd)==0:
        continue

      cmd = os.path.basename(cmd[0])
      if not string.lower(cmd) in ['ld', 'link', 'xlink', 'collect2']:
        continue

      # link line
      return parse_lib_str(line)
    return []

  # }}}

  # {{{ test_optimize()
  def test_optimize():
    # options to turn off optimization
    if platform=='Linux':
      if conf.env['COMPILER_CC'] in ['gcc', 'icc']:
        conf.env.append_value('CFLAGS_optoff', '-O0')
        conf.env.append_value('CXXFLAGS_optoff', '-O0')
    elif platform=='Windows':
      if conf.env['COMPILER_CC'] in ['msvc', 'icc']:
        conf.env.append_value('CFLAGS_optoff', '/Od')
        conf.env.append_value('CXXFLAGS_optoff', '/Od')

    # detect optimization option
    if platform=='Linux':
      if conf.env['COMPILER_CC'] in ['gcc', 'icc']:
        conf.start_msg('Detecting optimization options')

        oopts = ['-O2 -unroll -axSSE4.2,SSE4.1,SSSE3 -msse3',
                 '-O2 -unroll -axS -msse3',
                 '-O2 -unroll -msse3',
                 '-O2']
        if conf.options.cc_opt:
          oopts.insert(0, conf.options.cc_opt)

        for oopt in oopts:
           if test_opt(oopt): break
        conf.end_msg(oopt)
        conf.env.append_value('CFLAGS_opt', oopt.split())

        for oopt,msg in [('-fvisibility-inlines-hidden', 'Checking for visibility flags')]:
          conf.start_msg(msg)
          if test_opt(oopt, lang='cxx'):
            conf.end_msg('yes')
            conf.env.append_value('CXXFLAGS_opt', oopt.split())
          else:
            conf.end_msg('no')


        conf.env.append_value('CXXFLAGS_opt', conf.env.CFLAGS_opt)
        conf.env.append_value('FCFLAGS_opt', conf.env.CFLAGS_opt)
    elif platform=='Windows':
      if conf.env['COMPILER_CC'] in ['msvc', 'icc']:
        conf.start_msg('Detecting optimization options')
        oopts = ['/O2 /Qunroll /QaxSSE4.2,SSE4.1,SSSE3 /arch:sse3',
                 '/O2 /Qunroll /QaxS /arch:sse3',
                 '/O2 /Qunroll /arch:sse3',
                 '/O2']
        if conf.options.cc_opt:
          oopts.insert(0, conf.options.cc_opt)
        for oopt in oopts:
          if test_opt(oopt): break
        conf.end_msg(oopt)
        conf.env.append_value('CFLAGS_opt', oopt.split())
        conf.env.append_value('CXXFLAGS_opt', conf.env.CFLAGS_opt)
        conf.env.append_value('FCFLAGS_opt', conf.env.CFLAGS_opt)

    elif platform=='AIX':
      if conf.env['COMPILER_CC'] in ['xlc']:
        conf.start_msg('Detecting optimization options')
        oopts = ['-O5 -qmaxmem=-1 -qipa -qstrict -qarch=auto -qtune=auto',
                 '-O4 -qmaxmem=-1 -qipa -qstrict -qarch=auto -qtune=auto',
                 '-O3 -qmaxmem=-1 -qipa -qstrict -qarch=auto -qtune=auto',
                 '-O2 -qmaxmem=-1 -qipa -qarch=auto -qtune=auto']
        if conf.options.cc_opt:
          oopts.insert(0, conf.options.cc_opt)
        for oopt in oopts:
          if test_opt(oopt): break
        conf.end_msg(oopt)
        conf.env.append_value('CFLAGS_opt', oopt.split())
        conf.env.append_value('CXXFLAGS_opt', conf.env.CFLAGS_opt)
        conf.env.append_value('FCFLAGS_opt', conf.env.CFLAGS_opt)
  # }}}

  # {{{ test_debug()
  def test_debug():
    if platform=='Linux':
      if conf.env['COMPILER_CC'] in ['gcc', 'icc']:
        conf.check_cc(cflags='-g', msg='Checking for debugging support')
        conf.env.append_value('CFLAGS', '-g')
        conf.env.append_value('CXXFLAGS', '-g')
    elif platform=='Windows':
      if conf.env['COMPILER_CC'] in ['msvc', 'icc']:
        conf.check_cc(cflags='/Zi', msg='Checking for debugging support')
        conf.env.append_value('CFLAGS', '/Zi')
        conf.env.append_value('CXXFLAGS', '/Zi')
    elif platform=='AIX':
      if conf.env['COMPILER_CC'] in ['xlc']:
        conf.check_cc(cflags='-g', msg='Checking for debugging support')
        conf.env.append_value('CFLAGS', '-g')
        conf.env.append_value('CXXFLAGS', '-g')
  # }}}

  if conf.options.debug:
    test_debug()
  else:
    test_optimize()

  # {{{ check types
  def check_types(type,name=None):
    str='''
#include <stdio.h>
int main()
{
	printf("%%d", (int)sizeof(%s));
	return 0;
}''' % type

    if name==None: name = Utils.quote_define_name(type)
    name = 'SIZEOF_%s' % name

    size = conf.check_cc(fragment=str, msg='Checking for size of %s' % type,
                         execute=True, define_ret=True)
    conf.define(name, int(size))
  # }}}
  for t in ['double', 'float', 'int', 'long int', 'long long int', 'short int']:
    check_types(t)
  check_types('void *', 'VOID_P')

  # {{{ common headers
  for h in '''fcntl.h float.h fenv.h limits.h stddef.h stdlib.h
              string.h stdio.h assert.h sys/time.h sys/types.h
              sys/stat.h stdlib.h string.h memory.h strings.h
        		  inttypes.h stdint.h unistd.h'''.split():
    try:    conf.check(header_name=h, features='c cprogram')
    except: pass

  for h,d in [('tr1/unordered_map', 'HAVE_TR1_UNORDERED_MAP'),
              ('tr1/unordered_set', 'HAVE_TR1_UNORDERED_SET'),
              ('unordered_map', 'HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER'),
              ('unordered_set', 'HAVE_TR1_UNORDERED_SET_WITH_STD_HEADER'),
              ('limits', 'HAVE_STD_LIMITS')]:
    try:    conf.check(header_name=h, compile_mode='cxx', define_name=d)
    except: pass
  # }}}

  # {{{ namespace
  def check_namespace():
    str='''
namespace Outer { namespace Inner { int i = 0; }}
int
main ()
{
using namespace Outer::Inner; return i;
  ;
  return 0;
}'''
    conf.check_cxx(fragment=str, msg='Checking for c++ namespaces',
                     define_name='HAVE_NAMESPACES')
  # }}}
  try: check_namespace()
  except: pass


  # {{{ sstream
  def check_sstream():
    str='''
#include <sstream>
#ifdef HAVE_NAMESPACES
using namespace std;
#endif
int
main ()
{
stringstream message; message << "Hello"; return 0;
  ;
  return 0;
}'''
    conf.check_cxx(fragment=str, msg='Checking for std::sstream',
                   define_name='HAVE_SSTREAM')
  # }}}
  try: check_sstream()
  except: pass


  if not platform=='Windows':
    conf.check_cc(lib='m', uselib_store='MATH')

  conf.recurse('src/contrib/brkpnts')

  # {{{ Petsc

  def config_petsc():
    found = False
    base_dir = conf.options.petsc_dir
    arch_dir = os.path.join(base_dir, conf.options.petsc_arch)

    petsc_release = True
    try:
      petsc_ver=Utils.readf(os.path.join(base_dir, 'include', 'petscversion.h'))

      vers=[]
      for s in 'MAJOR MINOR SUBMINOR'.split():
        m = re.search('#define\s+PETSC_VERSION_'+s+'\s+(\d+)', petsc_ver)
        if m: vers.append(m.group(1))
      version = '.'.join(vers)

      m = re.search('#define\s+PETSC_VERSION_RELEASE'+'\s+(\d+)', petsc_ver)
      if m: petsc_release = (m.group(1) == 1)
    except:
      conf.fatal('Could not find petscversion.h, or it can not be parsed.')
    print 'Using Petsc version %s' % version

    petsc_vars=Utils.str_to_dict(Utils.readf(os.path.join(arch_dir, 'conf/petscvariables')))

    for v in ['PACKAGES_INCLUDES', 'PETSC_CC_INCLUDES']:
      if not petsc_vars.has_key(v): continue
      inc_dirs = parse_inc_str(petsc_vars[v])
    inc_dirs.append(os.path.join(base_dir, 'include'))
    inc_dirs.append(os.path.join(arch_dir, 'include'))

    conf.check_cxx(header_name='petscversion.h',
                   cxxflags=[conf.env.CPPPATH_ST % d for d in inc_dirs],
                   uselib_store='PETSC', define_name='HAVE_PETSC')

    lib_dirs, libs = [],[]

    linkflags = []
    for v in ['PETSC_WITH_EXTERNAL_LIB', 'PETSC_LIB_BASIC', 'PACKAGES_LIBS']:
      if not petsc_vars.has_key(v): continue
      linkflags.extend( parse_lib_str(petsc_vars[v]) )

    # MPI
    conf.start_msg('Checking for MPI')
    inc_mpi = petsc_vars.get('MPI_INCLUDE', 'mpiuni').strip()
    if inc_mpi.split(os.path.sep)[-1] == 'mpiuni':
        conf.end_msg('no')
    else:
        conf.end_msg('yes')
        conf.define('HAVE_MPI', 1)

        if platform=='Linux':
          # MPI libraries
          conf.start_msg('Checking MPI library')
          cmd = [petsc_vars['PCC'],
                 '-cc=%s'%conf.env['COMPILER_CC'],
                 '-show',
                ]
          out, err = conf.cmd_and_log(cmd, output=waflib.Context.BOTH)
          linkflags.extend( parse_lib_str(out) )
          conf.end_msg('ok')

        # AIX with PE, shoule use -binitfini:poe_remote_main as link flag
        if platform=='AIX':
          linkflags.append('-binitfini:poe_remote_main')

    if not platform=='Windows':
      # Fortran library needed when using C++ linker
      conf.start_msg('Checking Fortran library')
      cxx_flags = test_imp_lib(lang='cxx')
      f_flags = test_imp_lib(lang='f')
      flags = []
      for f in f_flags:
        if f.startswith('-L') or \
           f.startswith('-l') or \
           f.endswith('.a') or \
           f.endswith('.lib'):
          if f in cxx_flags: continue

        if platform=='Linux' and f=='-lifcore':
          # ICC 11.1 mistakenly uses libifcore.a instead of libifcore_pic.a
          flags.append('-lifcore_pic')
        else:
          flags.append(f)  # append to list only if the flag not present for CXX

      linkflags.extend(simplifyLinkerOptions(flags))
      conf.end_msg('ok')
    else:
      # some default windows libs
      flags = parse_lib_str(petsc_vars['PCC_LINKER_LIBS'])
      linkflags.extend(flags)

    conf.check_cxx(linkflags=linkflags,
                   uselib_store='PETSC', msg='Checking for library Petsc')
    conf.env.append_value('LINKFLAGS_PETSC', linkflags)

  # }}}
  config_petsc()

  # {{{ config_slepc()
  def config_slepc():
    base_dir = conf.options.slepc_dir
    arch_dir = os.path.join(base_dir, conf.options.petsc_arch)

    try:
      slepc_ver=Utils.readf(os.path.join(base_dir, 'include/slepcversion.h'))
      vers=[]
      for s in 'MAJOR MINOR SUBMINOR'.split():
        m = re.search('#define\s+SLEPC_VERSION_'+s+'\s+(\d+)', slepc_ver)
        if m: vers.append(m.group(1))
      version = '.'.join(vers)
    except:
      conf.fatal('Could not find petscversion.h, or it can not be parsed.')
    print 'Using Slepc version %s' % version

    inc_dirs = [os.path.join(base_dir, 'include')]
    inc_dirs.append(os.path.join(arch_dir, 'include'))
    conf.check_cxx(header_name='slepc.h',
                   cxxflags=[conf.env.CPPPATH_ST % d for d in inc_dirs],
                   use='PETSC',
                   uselib_store='SLEPC', define_name='HAVE_SLEPC')

    lib_dirs = [os.path.join(arch_dir, 'lib')]
    libs     = ['slepc']
    conf.check_cxx(lib=libs,
                   linkflags=[conf.env.LIBPATH_ST % d for d in lib_dirs],
                   uselib_store='SLEPC', msg='Checking for library Slepc')
  # }}}
  if conf.options.slepc_enabled:
    config_slepc()


# {{{ config_ams()
  def config_ams():
    base_dir = conf.options.ams_dir

    conf.check_cxx(header_name='ams.h',
                   cxxflags=[conf.env.CPPPATH_ST % os.path.join(base_dir, 'include')],
                   uselib_store='AMS', define_name='HAVE_AMS')

    libs     = ['amsacc']
    conf.check_cxx(lib=libs,
                   linkflags=[conf.env.LIBPATH_ST % os.path.join(base_dir, 'lib')],
                   uselib_store='AMS', msg='Checking for library AMS')

# }}}
  if conf.options.ams_enabled:
    config_ams()


  # {{{ NetGen
  def config_netgen():
    found = False
    search_dirs = [None, '/usr', '/usr/local', '/opt/netgen', '/usr/local/netgen']
    if conf.options.netgen_dir:
      search_dirs = [conf.options.netgen_dir]

    for ngdir in search_dirs:
      cxxflags, linkflags = '', ''
      if ngdir:
        cxxflags  = conf.env.CPPPATH_ST % os.path.join(ngdir,'include')
        linkflags = conf.env.LIBPATH_ST % os.path.join(ngdir,'lib')
      try:
        conf.check_cxx(header_name='nglib.h', cxxflags=cxxflags, uselib_store='NETGEN')
        conf.check_cxx(lib='nglib', linkflags=linkflags, uselib_store='NETGEN')
        found=True
        break
      except: pass
    #if not found:
    #  conf.fatal('Can not find the Netgen library!')
  # }}}
  config_netgen()

  # {{{ CGNS
  def config_cgns():
    found = False
    search_dirs = [None, '/usr', '/usr/local', '/usr/local/cgns']
    if conf.options.cgns_dir:
      search_dirs = [conf.options.cgns_dir]

    for cgnsdir in search_dirs:
      cxxflags, linkflags = '',''
      if cgnsdir:
        cxxflags  = conf.env.CPPPATH_ST % os.path.join(cgnsdir,'include')
        linkflags = conf.env.LIBPATH_ST % os.path.join(cgnsdir,'lib')
      try:
        conf.check_cxx(header_name='cgnslib.h', cxxflags=cxxflags, uselib_store='CGNS')
        conf.check_cxx(lib='cgns', linkflags=linkflags, uselib_store='CGNS')
        found=True
        break
      except: pass
    if not found:
      conf.fatal('Can not find the CGNS library!')

    conf.define('HAVE_CGNS', 1)
  # }}}
  config_cgns()

  # {{{ VTK
  def config_vtk():
    libs = '''vtkRendering vtkGraphics vtkImaging vtkIO vtkFiltering
              vtkCommon vtksys vtkDICOMParser vtkpng vtktiff vtkzlib
              vtkjpeg vtkexpat vtkftgl vtkfreetype'''.split()

    if platform=='Linux':
      libs.append('pthread')

    found = False
    search_dirs = [None, '/usr', '/usr/local']
    if conf.options.vtk_dir:
      search_dirs = [conf.options.vtk_dir]
    ver=conf.options.vtk_ver

    for vtkdir in search_dirs:
      cxxflags, linkflags = '',''
      if vtkdir:
        if len(ver)>0:
          cxxflags  = conf.env.CPPPATH_ST % os.path.join(vtkdir,'include', ver)
          linkflags = conf.env.LIBPATH_ST % os.path.join(vtkdir,'lib', ver)
        else:
          cxxflags  = conf.env.CPPPATH_ST % os.path.join(vtkdir,'include')
          linkflags = conf.env.LIBPATH_ST % os.path.join(vtkdir,'lib')
      try:
        conf.check_cxx(header_name='vtkConfigure.h', cxxflags=cxxflags, uselib_store='VTK')

        conf.check_cxx(lib=libs, linkflags=linkflags, uselib_store='VTK',
                       msg='Checking for libraries for VTK')
        found=True
        break
      except: pass
    if found:
      conf.define('HAVE_VTK', 1)

  # }}}
  config_vtk()


  # {{{ SIP
  def config_sip():
    conf.start_msg('Checking for python-sip')
    try:
      import sipconfig
      sipconf = sipconfig.Configuration()
      conf.load('siptool', tooldir='./build')
      conf.env['SIP_BIN']     = sipconf.sip_bin
      conf.env['INCLUDES_SIP'] = [sipconf.sip_inc_dir, sipconf.py_inc_dir]
      conf.env['LIBPATH_SIP']  = [sipconf.py_lib_dir]
      conf.end_msg('yes')
    except:
      conf.env['SIP_BIN'] = None
      conf.end_msg('no')
  # }}}
  config_sip()

  conf.load('myflex mybison', tooldir='./build')

  conf.write_config_header('config.h')
  #print conf.env

def build(bld):
  #print bld.env
  bld.contrib_objs =[]
  bld.recurse('src')


  bld.install_files('${PREFIX}/lib',
                    ['lib/GeniusSyntax.xml','lib/material.def',
                     'lib/solar_spectrum_outspace']
                   )

  bld.install_files('${PREFIX}/bin',
                    ['bin/GeniusCtrl.py',
                     'bin/GeniusLib.py',
                     'bin/HTTPFE.py']
                   )

  bld.install_files('${PREFIX}/bin',
                     ['bin/geniusd.py'],
                     chmod=Utils.O755
                   )

  bld.install_files('${PREFIX}', bld.path.ant_glob('examples/**'), relative_trick=True)

