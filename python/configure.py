#!/usr/bin/env python

import os
import sys
import shutil
import sipconfig

class GeniusConfiguration(sipconfig.Configuration):
    """Overriding some SIP configuration values
    """

    def __init__(self):
        win32_pkg_config = {
                'default_bin_dir':  '/cygdrive/c/usr/cogenda/Python-2.6.4',
                'default_mod_dir':  '/cygdrive/c/usr/cogenda/Python-2.6.4',
                'default_sip_dir':  '/cygdrive/c/usr/cogenda/Python-2.6.4',
                'platform':         'win32-g++',
                'py_conf_inc_dir':  '/cygdrive/c/usr/cogenda/Python-2.6.4/Include',
                'py_inc_dir':       '/cygdrive/c/usr/cogenda/Python-2.6.4/Include',
                'py_lib_dir':       '/cygdrive/c/usr/cogenda/Python-2.6.4/libs',
                'sip_bin':          'C:\\usr\\cogenda\\Python-2.6.4\\sip',
                'sip_config_args':  '',
                'sip_inc_dir':      '/cygdrive/c/usr/cogenda/Python-2.6.4/Include',
                'sip_mod_dir':      '/cygdrive/c/Python25/Lib/site-packages',
                }
        
	win32_macros = {
                'AIX_SHLIB':                '',
                'AR':                       '',
                'CC':                       'win32fe.bat icl',
                'CFLAGS':                   '-O1',
                'CFLAGS_CONSOLE':           '',
                'CFLAGS_DEBUG':             '',
                'CFLAGS_EXCEPTIONS_OFF':    '',
                'CFLAGS_EXCEPTIONS_ON':     '',
                'CFLAGS_MT':                '',
                'CFLAGS_MT_DBG':            '',
                'CFLAGS_MT_DLL':            '',
                'CFLAGS_MT_DLLDBG':         '',
                'CFLAGS_RELEASE':           '',
                'CFLAGS_RTTI_OFF':          '',
                'CFLAGS_RTTI_ON':           '',
                'CFLAGS_SHLIB':             '',
                'CFLAGS_STL_OFF':           '',
                'CFLAGS_STL_ON':            '',
                'CFLAGS_THREAD':            '',
                'CFLAGS_WARN_OFF':          '-w',
                'CFLAGS_WARN_ON':           '',
                'CHK_DIR_EXISTS':           'test -d',
                'CONFIG':                   'qt warn_on release incremental flat link_prl precompile_header autogen_precompile_source copy_dir_files debug_and_release debug_and_release_target embed_manifest_dll embed_manifest_exe',
                'COPY':                     'cp -f',
                'CXX':                      'win32fe.bat icl',
                'CXXFLAGS':                 '-O1',
                'CXXFLAGS_CONSOLE':         '',
                'CXXFLAGS_DEBUG':           '',
                'CXXFLAGS_EXCEPTIONS_OFF':  '',
                'CXXFLAGS_EXCEPTIONS_ON':   '',
                'CXXFLAGS_MT':              '',
                'CXXFLAGS_MT_DBG':          '',
                'CXXFLAGS_MT_DLL':          '',
                'CXXFLAGS_MT_DLLDBG':       '',
                'CXXFLAGS_RELEASE':         '',
                'CXXFLAGS_RTTI_OFF':        '',
                'CXXFLAGS_RTTI_ON':         '',
                'CXXFLAGS_SHLIB':           '',
                'CXXFLAGS_STL_OFF':         '',
                'CXXFLAGS_STL_ON':          '',
                'CXXFLAGS_THREAD':          '',
                'CXXFLAGS_WARN_OFF':        '-w',
                'CXXFLAGS_WARN_ON':         '',
                'DEFINES':                  'UNICODE WIN32 QT_LARGEFILE_SUPPORT',
                'DEL_FILE':                 'rm -f',
                'EXTENSION_PLUGIN':         '',
                'EXTENSION_SHLIB':          '',
                'INCDIR':                   '',
                'INCDIR_OPENGL':            '',
                'INCDIR_X11':               '',
                'LFLAGS':                   '',
                'LFLAGS_CONSOLE':           '',
                'LFLAGS_CONSOLE_DLL':       '',
                'LFLAGS_DEBUG':             '',
                'LFLAGS_DLL':               '',
                'LFLAGS_OPENGL':            '',
                'LFLAGS_PLUGIN':            '',
                'LFLAGS_RELEASE':           '',
                'LFLAGS_SHLIB':             '',
                'LFLAGS_SONAME':            '',
                'LFLAGS_THREAD':            '',
                'LFLAGS_WINDOWS':           '',
                'LFLAGS_WINDOWS_DLL':       '-LD',
                'LIB':                      '',
                'LIBDIR':                   '',
                'LIBDIR_OPENGL':            '',
                'LIBDIR_X11':               '',
                'LIBS':                     '',
                'LIBS_CONSOLE':             '',
                'LIBS_CORE':                'kernel32.lib user32.lib shell32.lib uuid.lib ole32.lib advapi32.lib ws2_32.lib',
                'LIBS_GUI':                 'gdi32.lib comdlg32.lib oleaut32.lib imm32.lib winmm.lib winspool.lib ws2_32.lib ole32.lib user32.lib advapi32.lib',
                'LIBS_NETWORK':             'ws2_32.lib',
                'LIBS_OPENGL':              'opengl32.lib glu32.lib gdi32.lib user32.lib',
                'LIBS_RT':                  '',
                'LIBS_RTMT':                '',
                'LIBS_THREAD':              '',
                'LIBS_WINDOWS':             '',
                'LIBS_X11':                 '',
                'LINK':                     'win32fe.bat icl',
                'LINK_SHLIB':               '',
                'LINK_SHLIB_CMD':           '',
                'MAKEFILE_GENERATOR':       'MINGW',
                'MKDIR':                    'mkdir -p',
                'RANLIB':                   'ranlib -s',
                'RPATH':                    '',
                'STRIP':                    ''
        }


        if sys.platform == 'win32':
            super(GeniusConfiguration, self).__init__([win32_pkg_config])
            self._macros = win32_macros
        else:
            super(GeniusConfiguration, self).__init__()

class GeniusModuleMakefile(sipconfig.SIPModuleMakefile):
    def __init__(self, configuration, build_file, install_dir=None, static=0,
                 console=0, qt=0, opengl=0, threaded=0, warnings=1, debug=0,
                 dir=None, makefile="Makefile", installs=None, strip=1,
                 export_all=1):
        sipconfig.SIPModuleMakefile.__init__(self, configuration, build_file, install_dir, static,
                 console, qt, opengl, threaded, warnings, debug,
                 dir, makefile, installs, strip,
                 export_all)
        if sys.platform == 'win32':
            self._manifest=None


    def generate_target_install(self, mfile):
        mfile.write("\ninstall: $(TARGET)\n")
        if sys.platform=='win32':
          self.install_file(mfile, "$(TARGET)", "$(GENIUS_DIR)/lib/", strip=0)
        else:
          self.install_file(mfile, "$(TARGET)", "$(GENIUS_DIR)/lib", strip=0)

    def generate_macros_and_rules(self, mfile):
        mfile.write("include ../../make.defs\n")
        sipconfig.SIPModuleMakefile.generate_macros_and_rules(self, mfile)
        mfile.write("LIBS += %s\n" % ' '.join(self.extra_extra_lflags))


tmp_dir = './build-tmp'
build_dir = './build'
if os.path.exists(tmp_dir):
    shutil.rmtree(tmp_dir)
os.mkdir(tmp_dir)

if not os.path.exists(build_dir):
  os.mkdir(build_dir)

###############
# generate in a temporary directory first
###############
# The name of the SIP build file generated by SIP and used by the build
# system.
build_file = tmp_dir+"/genius.sbf"

# Get the SIP configuration information.
config = GeniusConfiguration()

# Run SIP to generate the code.
os.system(" ".join([config.sip_bin, "-c", tmp_dir, "-b", build_file, "sip/genius.sip"]))

print build_file

# Create the Makefile.
makefile = GeniusModuleMakefile(config, build_file, debug=0, makefile=tmp_dir+'/Makefile')

# Add the library we are wrapping.  The name doesn't include any platform
# specific prefixes or extensions (e.g. the "lib" prefix on UNIX, or the
# ".dll" extension on Windows).
makefile.extra_cflags=['$(PETSC_INCLUDE)', '$(GENIUS_INCLUDE)']
makefile.extra_cxxflags=['$(PETSC_INCLUDE)', '$(GENIUS_INCLUDE)']
if sys.platform=='linux2':
    makefile.extra_lflags=['-Wl,--export-dynamic']

makefile.extra_extra_lflags=['-L$(GENIUS_DIR)/src', '-lgenius', '$(PETSC_LIB)', '$(CLIBS)']
# due to a strange problem on win32, we add -g flag to everything
#if sys.platform=='win32':
#    makefile.extra_cflags.append('-g')
#    makefile.extra_cxxflags.append('-g')
#    makefile.extra_lflags.append('-g')

# Generate the Makefile itself.
makefile.generate()



##############
# merge back to the actual build directory
##############
os.system("diff -uwr --strip-trailing-cr --unidirectional-new-file -I 'Generated by SIP' build build-tmp > build-tmp.diff")
os.chdir('build')
os.system("patch -p1 < ../build-tmp.diff")
os.chdir('..')

#if os.path.exists(tmp_dir):
#    shutil.rmtree(tmp_dir)
