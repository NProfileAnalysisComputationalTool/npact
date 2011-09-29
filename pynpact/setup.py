#!/usr/bin/env python
import os, os.path, glob

from distutils.command.build import build
from distutils.core import setup
from distutils.sysconfig import customize_compiler
import distutils.ccompiler


#need to make sure all the C files that are in src get compiled and
#put into ./pynpact/bin
src_path= os.path.join(os.path.dirname(__file__),'src')
c_files=glob.glob(os.path.join(src_path, "*.c"))

class build_pynpact(build) :
    def run(self) :
        # Setup the CCompiler object that we'll use to do all the
        # compiling and linking
        self.compiler = distutils.ccompiler.new_compiler(compiler=self.compiler,
                                                         verbose=self.verbose,
                                                         dry_run=self.dry_run,
                                                         force=self.force)
        customize_compiler(self.compiler)
        objects = self.compiler.compile(c_files,
                                        output_dir=self.build_temp,
                                        debug=self.debug)
        output_dir = os.path.join(self.build_platlib,'pynpact','bin')
        for o in objects :
            basename = os.path.splitext(os.path.basename(o))[0]
            exe_name = self.compiler.executable_filename(basename,strip_dir=True)
            self.compiler.link_executable([o],
                                          exe_name,
                                          libraries=['m'],
                                          output_dir=output_dir,debug=self.debug)


        build.run(self)
    
    # def install(self):
    #     if os.path.isdir(self.build_dir):
    #         outfiles = self.copy_tree(, self.install_dir)
    #     

setup(name='pynpact',
      version='0.2',
      description='Python N-Profile Analysis Computation Tool',
      author='Luciano Brocchieri and Nathan Bird',
      author_email='nathan@acceleration.net',
      url='http://genome.ufl.edu/spat',
      packages=['pynpact'],

      requires=["biopython(>=1.57)"],
      #data_files = ['bin',[]],
      cmdclass={'build': build_pynpact }
      
     )
