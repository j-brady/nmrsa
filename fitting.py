import yaml
import numpy as np
""" Functions for fitting nmr data """

def load_yaml(yaml_file):
    """ reads files containing YAML and converts to dictionary """
    f = open(yaml_file)
    y = f.read()
    return yaml.load(y)


class R1rho:
    """ Simple exponential function "func" = I0*np.exp(-t*R1) """   

    def __init__(self):
        pass

    def func(self,t,I0,R1):
        return I0*np.exp(-t*R1)
          
    def getDelays(self,yaml_file="time_relax_list.yaml"):
        """ reads yaml file containing list of T values and
            returns a numpy array of the T values """
        values = load_yaml(yaml_file)
        values = np.array(values["T values"])
        return values

    def calcR2(self,R1,Omega_s,nu_1s):
        #theta = arccot(Omega_s/nu_1s)
        #nu_1s = 15N spin lock field strength in Hz
        #Omega_s is the resonance offset from the spin lock carrier
        pass


class Diffusion:
    """ 
        For generating the diffusion equation
        
        parameters:

        T_diff  -- the diffusion time in seconds
        delta   -- the length of the gradient in seconds 
        Dtype   -- whether the diffusant is triple (3Q), double (2Q) or single (SQ).
                  This adds a Qfactor of 3,2 or 1 respectively to the (Qfactor*gamma*G*delta)^2 term.
        bipolar -- True/False, if bipolar gradients are used delta=2*delta

    """

    def __init__(self,T_diff,delta,Dtype="1Q",bipolar=True):
        if Dtype == "3Q":
            self.Qfactor = 3.
        elif Dtype == "2Q":
            self.Qfactor = 2.
        elif Dtype == "1Q":
            self.Qfactor = 1.
        
        self.T_diff = T_diff  # diffusion time in s

        if bipolar:
            self.delta = delta*2.    # gradient duration in s 
        else:
            self.delta = delta

    def func(self,x,I0,D):
        """ 
            x  -- is gradient strength squared!!
            I0 -- intensity when 0 gradient strength is applied
            D  -- diffusion constant

            gamma is in units of rad s^-1 G^-1

            returns:

                I0*exp(-1.*D*(Qfactor*G*gamma*delta)^2*(T_diff-delta/3.))

                where G is the gradient strength, Qfactor is the factor associated with the coherence
                level of the diffusant (i.e. 1,2 or 3), delta is the gradient length (s), T_diff is the
                diffusion time (s) and D is the diffusion constant (cm^2 s^-1).

        """

        gamma = 2.67513e4 # rads-1 G-1 (267.513e6 rads-1 T-1 * 1e-4 (Conversion to Gauss))
        return I0*np.exp(-1.*D*x*np.square(self.Qfactor*gamma*self.delta)*(self.T_diff-self.delta/3.))
