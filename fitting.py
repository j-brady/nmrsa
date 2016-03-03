import yaml
import numpy as np
""" Functions for fitting nmr data """

def load_yaml(yaml_file):
    """ reads files containing YAML and converts to dictionary """
    f = open(yaml_file)
    y = f.read()
    return yaml.load(y)

class R1rho:
    def __init__(self):
        pass

    def func(self,t,I0,R1):
        return I0*np.exp(-t*R1)
          
    def getDelays(self,yaml_file="time_relax_list.yaml"):

        values = load_yaml(yaml_file)
        values = np.array(values["T values"])
        return values

class Diffusion:
    pass
