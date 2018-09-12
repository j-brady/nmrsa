from numpy import log10

def pl(Tp,B2):
    # Tp = pulse width in seconds at given power level to
    # B2 is the target field strength
    B1 = 1./(4.*Tp) # B1 field
    print("B1 field = %.3e")% B1
    power_lev = -20.*log10(B2/B1)
    print("Power level for %.3f kHz field is %.3f dB" % (B2/1e3,power_lev))
    return power_lev

if __name__ == "__main__":
    #90 pulse power
    power = 9.54 #dB 
    #power = -19.56 #dB 
    #pl(7.35e-6,10e3)
    #pl(7.35e-6,30)
    #pl(7.18e-6,30)
    #pl(7.18e-6,10e3)
    print(pl(7.125e-6,30)-9.54)
    #power2 = pl(7.22e-6,30)
    #power2 = pl(12.85e-6,1e3)
    #power2 = pl(7.1e-6,30)
    #print -power2-power
