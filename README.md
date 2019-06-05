# nmrsa
NMR scriptaholics anonymous

Scripts for processing and fitting NMR data written by Jacob Brady and Rui Huang.

More documentation is on its way...

## Installation

Still in the stone age...

download the code

cd into top directory and run `python setup.py install`. Probably best to use a virtualenv!

This should install the package and the executable scripts `makeYamlFiles.py` and `proc_with_yaml.py`. 

## Processing bruker diffusion data (requires NMRPipe to be installed)

Mainly use for processing bruker PFG diffusion NMR data. First download your data dirs containing your bruker data that you want to batch process into your working directory (assuming each dataset is a pseudo-2D).

Next run `makeYamlFiles.py`. This should automatically generate an example NMRPipe `ft.com` script along with a `proc.yaml` file and `gradients.yaml` file containing script parameters and gradient values. To get your NMRPipe `fid.com` file you need to run the `bruker` command in NMRPipe. 

## Scripts

* makeYamlFiles.py
* proc_with_yaml.py
* spec.py - overlay NMR spectra using matplotlib and nmrglue
