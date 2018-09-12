from distutils.core import setup

setup(name='nmrsa',
      version='',
      description='NMR processing scripts',
      author='Jacob Brady and Rui Huang',
      packages=['nmrsa'],
      scripts=['scripts/makeYamlFiles.py', 'scripts/proc_with_yaml.py'],
     )
