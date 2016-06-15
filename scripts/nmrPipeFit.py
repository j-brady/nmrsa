#!/usr/bin/env python

from fuda import *

#import fudaIO, Numeric, dataIO, sys, string, os, math
import fudaIO, dataIO, sys, string, os, math
#
# To measure intensities in an arried 2D experiment,
# e.g. relaxation recovery, 15N T1 and T2 exp, CPMG, 
#
# 2002-2008
#                                      D. Flemming Hansen
#                                      flemming@pound.med.utoronto.ca
#
#########################################################################
#
# this the the default values
PARAMETERFILE=['varian','procpar']
GNUPLOT="Y"
FITLIN="N"
FITEXP="Y"
FIT180EXP="N"
FITBIEXP="N"
LOGFILENR="N"  #Non-relaxation logfile (2D NOESY/COSY, 3D NOESY, 3D J(CaCd) ..) 
#
#  To do list:
#     - Incorporate real Lorentzian instead of 'glore|c=0 '
#     - Fix problem with fitting 2D spectra (non array )
#     - Fix problem regarding including peaklists, where the dimension does not fit
#       the dimension of the spectrum, e.g. including a 3D peaklist for fitting 2D array
#
#########################################################################
#
if not len(sys.argv) == 3: 
    ProgramBaseName=(string.split(sys.argv[0],'/'))
    ProgramBaseName=ProgramBaseName[len(ProgramBaseName)-1]
    sys.stderr.write("\n Wrong number of parameters \n")
    sys.stderr.write("                                      PROGRAM ABORTED\n\n")    
    sys.stderr.write(" USAGE: %s [inputfile] [output_directory] \n" % (ProgramBaseName))
    sys.stderr.write(" The format of the inputfile can be seen elsewhere\n" )
    sys.stderr.write(" [output_directory] is the directory where the output files are stored\n")
    sys.stderr.write("                                           D. Flemming Hansen \n")
    sys.stderr.write("                                           d.hansen@ucl.ac.uk \n")
    sys.exit(2)
#
#########################################################################
# Intrinsic Functions
#########################################################################
def read_procpar_param(procparfile,param,type):
    if(type[0].lower()=='v'):
        return read_procpar_param_varian(procparfile,param)
    elif(type[0].lower()=='b'):
        return read_procpar_param_bruker(procparfile)
    else:
        sys.stderr.write(" The PARAMETERFILE type: -->%s<-- \n" % (type))
        sys.stderr.write(" is not valid\n PROGRAM ABORTED\n")
        sys.exit(1)
            
def read_procpar_param_bruker(procparfile):
    out=[]
    try:
        ifs=open(procparfile,'r')
    except (IOError):
        sys.stderr.write(" Problems while opening the varian parameter file %s\n" % (procparfile))
        sys.stderr.write("                                        PROGRAM ABORTED\n")        
        sys.exit(10)
    lines=ifs.readlines()
    ifs.close()
    for l in range(len(lines)):
        its=string.split(lines[l])
        if(len(its)!=1):
            sys.stderr.write(" Format of the Bruker parameter file is not recognized \n")
            sys.stderr.write(" format is:\n param[0]\nparam[1]\n...\n\n")
            sys.stderr.write(" PROGRAM ABORTED\n")
            sys.exit(1)
        out.append(its[0])
    return out
#    
def read_procpar_param_varian(procparfile,param):
    out=[]
    try:
        ifs=open(procparfile,'r')
    except (IOError):
        sys.stderr.write(" Problems while opening the varian parameter file %s\n" % (procparfile))
        sys.stderr.write("                                        PROGRAM ABORTED\n")        
        sys.exit(10)
    lines=ifs.readlines()
    ifs.close()
    for l in range(len(lines)):
        its=string.split(lines[l])
        if ( param == its[0] ):
            array=string.split(lines[l+1])
            for i in range(1,len(array)):
                out.append(array[i])
            #
            # Do some checking
            if not ( len(out) == int(array[0]) ):
                sys.stderr.write(" Something is wrong with the procpar file %s \n" % (procparfile))
                sys.stderr.write(" Expected to find %d parameters, but found %d \n" % (int(array[0]),len(out)))
                sys.stderr.write("                                        PROGRAM ABORTED\n")
                sys.exit(10)
            else:
                break
    return out
#
def PrintArray(ofs,array,format):
    ofs.write("{")
    for a in range(len(array)):
        try:
            len(array[a])
            ofs.write("{")
            for b in range(len(array[a])):
                ofs.write(format % (array[a][b]))
                if not b==len(array[a])-1:
                    ofs.write(",")
            ofs.write("}")
            if not a==len(array)-1:
               ofs.write(",")
        except (TypeError):
           ofs.write(format % (array[a])) 
           if not a==len(array)-1:
               ofs.write(",")
    ofs.write("}\n")
    
def read_input_line(line):
    param=string.upper(string.split(string.strip(line),'=')[0])
    #
    # Did we find an '=' sign
    if ( len(string.split(string.strip(line),'=')) < 2 ):
        sys.stderr.write(" The line:---->%s<----- is not of the correct format, i.e.,\n PARAMETER=VALUE\n" % ( line))
        sys.stderr.write("                                        PROGRAM ABORTED\n")
        sys.exit(1)
    if (len(string.split(string.strip(line),'=')) == 2):
        value=string.split(string.strip(line),'=')[1]
    if not (len(string.split(string.strip(line),'=')) == 2):
        value=string.split(string.split(string.strip(line),'(')[1],')')[0]

    return [param,value]
#
def peak_name_to_no(name):
    number=-1
    for i in range (0,len(peakline)):
        if (name==peakline[i][0]):
            number=i
    return number
#
def Dummy(i):
    return i
#
# Read the parameters and spectrum
#########################################################################

inputfile = open(sys.argv[1],'r')

#Initialize the input
par=''
value=''
specfile=''
noise=0.0
peaklistfile=''
shape='GLORE'
delayfactor=1.0
isotopeshift='n'
baseline='n'
printdata='y'
write_individual=1  #Write individual output files for each peak
verboselevel=1      #Print errors
lmtol=1e-12
lmmaxfev=250
dumpparameters=0
restoreparameters=0
parameterfile=PARAMETERFILE

def_linewidth=[]
def_radius=[]
#
#Values for the not default peaks.
discard_peaks=[]
discard_peaks_lg=0
include_peaks=[]
include_peaks_lg=0
not_def_peak=[]
overlap_peaks=[]
discard_slices=[]
discard_slices_lg=0
#
for i in range (0,3):
    def_linewidth.append(0.)
    def_radius.append(0.)
#
#Read and assign input parameters
for line in inputfile.readlines():
    if ( not line[0] == "#") and ( len(line) > 2 ) : 
        line=string.split(line,'#')[0]
        par=read_input_line(line)[0]
        value=read_input_line(line)[1]
        #print ('%s%s%s%s' %('-->',par,'<->',value))
        if (par == "SPECFILE"):
            specfile=string.strip(value)            
        elif (par == "NOISE"):
            noise=float(value)
        elif (par == "PEAKLIST"):
            peaklistfile=string.strip(value)
        elif (par == "ZCOOR"):
            zcoor=string.strip(value)
        elif (par == "PARAMETERFILE" ):
            parameterfile=[]
            value=string.split(string.split(value,'(')[1],')')[0]
            for its in string.split(value,';'):
                parameterfile.append(its)
            parameterfile[0]=parameterfile[0].lower()
        elif ( par == "FITLIN"):
            FITLIN=string.strip(value)
        elif ( par == "GNUPLOT"):
            GNUPLOT=string.strip(value)
        elif ( par == "FITEXP" ):
            FITEXP=string.strip(value)
        elif ( par == "FITBIEXP" ):
            FITBIEXP=string.strip(value)
        elif ( par == "FIT180EXP" ):
            FIT180EXP=string.strip(value)
        elif ( par == "SHAPE" ):
            shape=string.strip(value)
        elif ( par == "DELAYFACTOR" ):
            delayfactor=float(value)
        elif ( par == "ISOTOPESHIFT" ):
            if ( string.strip(value) == "N" or string.strip(value) == "n" ):
                isotopeshift='n'
            else:
                value=string.split(string.split(value,'(')[1],')')[0]
                isotopeshift=[]
                for numbers in string.split(value,','):
                    if ( numbers == "FIX" ):
                        isotopeshift.append('fix')
                    else:
                        isotopeshift.append(float(numbers))
        elif ( par == "BASELINE" ):
            baseline=string.strip(value)
        elif ( par == "PRINTDATA" ):
            if ( value[0] == "y" or value[0]=="Y" ):
                printdata='y'
            if ( value[0] == "n" or value[0] =="N" ):
                printdata='n'
        elif ( par == "VERBOSELEVEL" ):
            verboselevel=int(value)
        elif ( par == "DUMPPARAMETERS" ):
            if ( value[0] == "y" or value[0] == "Y" ):
                dumpparameters=1
            else:
                dumpparameters=0
        #
        #Read Spectral default information
        #
        elif (par == "DEF_LINEWIDTH_F1"):
            def_linewidth[0]=float(value)
        elif (par == "DEF_LINEWIDTH_F2"):
            def_linewidth[1]=float(value)
        elif (par == "DEF_LINEWIDTH_F3"):
            def_linewidth[2]=float(value)
        elif (par == "DEF_RADIUS_F1"):
            def_radius[0]=float(value)
        elif (par == "DEF_RADIUS_F2"):
            def_radius[1]=float(value)
        elif (par == "DEF_RADIUS_F3"):
            def_radius[2]=float(value)
        #
        # Read information for not default peaks
        #
        elif (par == "NOT_DEF_PEAK" ):
            #
            #start to fill the array with default values
            not_def_linewidth_f1=def_linewidth[0]
            not_def_linewidth_f2=def_linewidth[1]
            not_def_linewidth_f3=def_linewidth[2]
            not_def_radius_f1=def_radius[0]
            not_def_radius_f2=def_radius[1]
            not_def_radius_f3=def_radius[2]
            not_def_isotopeshift=isotopeshift
            not_def_shape=shape
            not_def_intensity=1.
            not_def_restoreparameters=0
            #
            found_name=0
            for not_def_line in string.split(value,';'):
                par=read_input_line(not_def_line)[0]
                val=read_input_line(not_def_line)[1]
                if (par == "NAME"):
                    not_def_name=val
                    found_name=1
                elif (par == "LINEWIDTH_F1"):
                    not_def_linewidth_f1=float(val)
                elif (par == "LINEWIDTH_F2"):
                    not_def_linewidth_f2=float(val)
                elif (par == "LINEWIDTH_F3"):
                    not_def_linewidth_f3=float(val)
                elif (par == "RADIUS_F1"):
                    not_def_radius_f1=float(val)
                elif (par == "RADIUS_F2"):
                    not_def_radius_f2=float(val)
                elif (par == "RADIUS_F3"):
                    not_def_radius_f3=float(val)
                elif (par == "SHAPE"):
                    not_def_shape=val
                elif (par == "ISOTOPESHIFT" ):
                    if ( string.strip(val) == "n" or string.strip(val) == "N" ):
                        not_def_isotopeshift='n'
                    else:
                        not_def_isotopeshift=[]
                        val=string.split(val,',')
                        for numbers in val: 
                            if ( numbers == "FIX" ):
                                not_def_isotopeshift.append('fix')
                            elif ( numbers[0:4] == "TRIP" ):
                                not_def_isotopeshift.append('triplet')
                            else:
                                not_def_isotopeshift.append(float(numbers))
                elif (par == "INTENSITY"):
                    not_def_intensity=float(val)
                elif ( par == "RESTOREPARAMETERS" ):
                    if ( ( string.strip(val) == "Y" ) or (string.strip(val) == "y")) :
                        not_def_restoreparameters=1
                    else:
                        not_def_restoreparameters=string.strip(val)
                else:
                    sys.stderr.write("The parameter \'%s\' is not valid for input in NOT_DEF_PEAK\n" % (par))
                    sys.stderr.write("                                        PROGRAM ABORTED\n")
                    sys.exit(1)   
            if not ( found_name ):
                sys.stderr.write(" The NOT_DEF_PEAK line:\n %s does not contain a peak name, i.e. NAME=...\n" % (line ))
                sys.stderr.write("                                        PROGRAM ABORTED\n")
                sys.exit(1)
            #
            #Put the not_def values in the array: "not_def_peak
            temp={}
            temp['id']=not_def_name
            temp['lw1']=not_def_linewidth_f1
            temp['lw2']=not_def_linewidth_f2
            temp['lw3']=not_def_linewidth_f3
            temp['rad1']=not_def_radius_f1
            temp['rad2']=not_def_radius_f2
            temp['rad3']=not_def_radius_f3   
            temp['isotopeshift']=not_def_isotopeshift
            temp['shape']=not_def_shape
            temp['intensity']=not_def_intensity
            temp['restoreparameters']=not_def_restoreparameters
            not_def_peak.append(temp)
            temp={}

            """
            NextSlot=len(not_def_peak)
            ##
            peak_fit_household(NextSlot,'id','save',not_def_peak,not_def_name)
            peak_fit_household(NextSlot,'lw1','save',not_def_peak,not_def_linewidth_f1)
            peak_fit_household(NextSlot,'lw2','save',not_def_peak,not_def_linewidth_f2)
            peak_fit_household(NextSlot,'lw3','save',not_def_peak,not_def_linewidth_f3)            
            peak_fit_household(NextSlot,'rad1','save',not_def_peak,not_def_radius_f1)
            peak_fit_household(NextSlot,'rad2','save',not_def_peak,not_def_radius_f2)
            peak_fit_household(NextSlot,'rad3','save',not_def_peak,not_def_radius_f3)            
            peak_fit_household(NextSlot,'isotopeshift',not_def_peak,not_def_isotopeshift)
            peak_fit_household(NextSlot,'shape',not_def_peak,not_def_shape)
            peak_fit_household(NextSlot,'intensity',not_def_peak,not_def_intensity)
            peak_fit_household(NextSlot,'restoreparameters',not_def_peak,not_def_restoreparameters)
            #
            #
            not_def_peak.append([not_def_name,not_def_linewidth_f1,not_def_linewidth_f2,not_def_radius_f1,not_def_radius_f2,not_def_isotopeshift,not_def_shape,not_def_intensity,not_def_restoreparameters])
            """

            #
            #
        elif (par == "OVERLAP_PEAKS" ):
            #Get rid of brackets
            value=string.split(string.split(value,'(')[1],')')[0]
            overlap_group=[]
            disjunkt=1   #Flag to test for disjunkt
            #
            for overlap_line in string.split(value,';'):
                #Check wether the name exists in 'overlap_peaks':
                for checkline in overlap_peaks:
                    for checkpeak in checkline:
                        if (checkpeak == overlap_line):
                            disjunkt=0
                overlap_group.append(overlap_line)
            #
            #Store the overlap_group in the array overlap_peaks
            if (disjunkt == 1):
                overlap_peaks.append(overlap_group)
            if (disjunkt == 0):
                print 'The overlap group: \n'
                print value
                print 'is not disjunkt to the previous groups'
                sys.exit(2)
        elif (par == "DISCARD_SLICES" ):
            value=string.split(string.split(value,'(')[1],')')[0]
            for slice in string.split(value,';'):
                discard_slices.append(slice)
            discard_slices_lg=1
        elif (par == "DISCARD_PEAKS" ):
            #Get rid of brackets
            value=string.split(string.split(value,'(')[1],')')[0]
            for peak in string.split(value,';'):
                discard_peaks.append(peak)
            discard_peaks_lg=1
        elif ( par == "INCLUDE_PEAKS" ):
            #
            # Check that we have brackets
            if ( len(string.split(value,'(')) == 1 ):
                sys.stderr.write(" Do not know the format %s for including peaks \n USAGE: INCLUDE_PEAKS=(NAME1;NAME2;....)\n" % (value))
                sys.stderr.write("                                        PROGRAM ABORTED\n")
                sys.exit(1)
            #
            #Get rid of brackets
            value=string.split(string.split(value,'(')[1],')')[0]
            for peak in string.split(value,';'):
                include_peaks.append(peak)
            include_peaks_lg=1
        elif ( par == "LM" ):
            for field in string.split(value,';'):
                lmparam=string.split(field,'=')[0]
                lmval=string.split(field,'=')[1]                
                if ( lmparam == "MAXFEV" ):
                    lmmaxfev=int(lmval)
                elif ( lmparam == "TOL" ):
                    lmtol=float(lmval)
                else:
                    sys.stderr.write(" The LM parameter \'%s\' is unknown\n PROGRAM ABORTED \n" % (lmparam,));
                    sys.exit(1)
        else:
            sys.stderr.write("The following inputline is not meaningful:\n")
            sys.stderr.write("%s" % (line))
            sys.stderr.write("\n PROGRAM ABORTED\n")
            sys.exit(1)
#                
#Check if all the nessecary parameters have been read:
if ( (specfile=='') or (noise==0.0) or (peaklistfile=='') or (def_linewidth[0]==0) or (def_linewidth[1]==0) ):
    sys.stderr.write('\n You must specify both "SPECFILE","NOISE", DEF_LINEWIDTH_Fi, and\n')
    sys.stderr.write('"PEAKLIST".\n')
    sys.stderr.write('                                               PROGRAM  ABORTED\n')
    sys.exit(1)
if ( include_peaks_lg and discard_peaks_lg ):
    sys.stderr.write("\n Cannot both have a DISCARD_PEAKS and a INCLUDE_PEAKS variable\n");
    sys.stderr.write("                                               PROGRAM  ABORTED\n");
    sys.exit(1)
if ( discard_slices_lg == 1 and zcoor=="2D" ):
    sys.stderr.write("\n It makes no sence to discard planes when you only have one \n")
    sys.stderr.write(" spectrum .. Stupido! \n")
    sys.stderr.write("                                              PROGRAM  ABORTED\n");
    sys.exit(1)
#
inputfile.close()
#
os.system('mkdir %s >& /dev/null ' % (sys.argv[2],))
#
path=''
specfile=path+specfile
fd=fudaIO.open(specfile)
#
#Read the Limits:
sfrq=[]
first=[]
last=[]
os.system(' showhdr -verb %s | awk \'{ if ( $3=="DIM:" || $2=="FIRST:" || $2=="LAST:" || $1=="OBS" ) { print $3,$4,$5 }} \' > templine' % (specfile,))
tempfile=open('templine','r')
#
items=string.split(tempfile.readline())
DIM=int(items[1])
print DIM

items=string.split(tempfile.readline())
if not zcoor=="3D":
    if zcoor=="2D":
        sfrq=[float(items[1]),float(items[0])]
    else:
        sfrq=[float(items[1]),float(items[0]),float(items[2])]
else:
    sfrq=[float(items[2]),float(items[1]),float(items[0])]
#       
items=string.split(tempfile.readline())
if not zcoor=="3D":
    if zcoor=="2D":
        first=[float(items[1]),float(items[0])]
    else:
        first=[float(items[1]),float(items[0]),float(items[2])]
else:
    first=[float(items[2]),float(items[1]),float(items[0])]
#
items=string.split(tempfile.readline())
if not zcoor=="3D":
    if zcoor=="2D":
        last=[float(items[1]),float(items[0])]
    else:
        last=[float(items[1]),float(items[0]),float(items[2])]
else:
    last=[float(items[2]),float(items[1]),float(items[0])]
#
tempfile.close()
#
#Clean up
os.system('rm -f templine')
#tempfile.close()
#
#If an arrayed 2D spectrum, then find the recovery delay
if ( not ( zcoor == "2D" or zcoor =="3D" ) ): 
    #find directoryname:
    if ( specfile[0] == '/' ):
        tree=string.split(string.strip(specfile[1:len(specfile)]),'/')
        dirname='/'
        for level in range(len(tree)-1):
            dirname=dirname+tree[level]+'/'
    else:
        dirname=""
        tree=string.split(string.strip(specfile),'/')
        for level in range(len(tree)-1):
            dirname=dirname+tree[level]+'/'
    #
    # Did we specify a directory name in the PARAMETERFILE?
    if (parameterfile[1].find('/')<0):
        procparfile=dirname+parameterfile[1]
    else:
        procparfile=parameterfile[1]
    #
    #                        file       name  type (bruker/varian)
    tau=[]
    line=read_procpar_param(procparfile,zcoor,parameterfile[0])
    if not ( len(line) == fudaIO.get_size(fd,2) ):        
        sys.stderr.write("\n It looks like the parameter \'%s\' does not have the same dimension\n as the spectrum.\n" % (zcoor) )
        sys.stderr.write(" Spectrum has %d planes, while \'%s\' has dimension of %d \n" % (fudaIO.get_size(fd,2),zcoor,len(line)))
        sys.stderr.write("                                             PROGRAM ABORTED\n")
        sys.exit(1)
    for i in range (0,fudaIO.get_size(fd,2)):
        tau.append(delayfactor*float(line[i]))
    #
    #Find max and min tau
    taumin=-1
    taumax=-1
    #
    #Arrange tau:
    newtau=[]
    for k in range(len(tau)):
        discard_tau=0
        for slice in discard_slices:
            if ( int(slice)==k+1 ):
                discard_tau=1
        if not ( discard_tau ):
            newtau.append(tau[k])
    tau=newtau
    #
    temp=[]
    tausort=[]
    for i in range (0,fudaIO.get_size(fd,2)-len(discard_slices)):
        temp.append(tau[i])
    #
    for k in range (0,fudaIO.get_size(fd,2)-len(discard_slices)):
        for i in range (0,fudaIO.get_size(fd,2)-len(discard_slices)):
            if ( tau[i]==max(temp) ):
                tausort.append(i)
                temp[i]=0
    #
    taumax=tausort[0]
    taumin=tausort[len(tausort)-1]
else:
    if ( verboselevel > 2 ):
        sys.stdout.write(" The fitting of exponential functions to intensities is switched off,\n")
        sys.stdout.write(" because only one intensity is available\n")
    FITLIN="N"
    FITEXP="N"
    FIT180EXP="N"
    FITBIEXP="N"
    LOGFILENR="Y"
#
# Open the fitlog files and write header
if (FITEXP == "Y"):
    fitlogfile=sys.argv[2]+"/singleexp.fit"
    fitlog = open(fitlogfile,'w')
    fitlog.write("# %-13s %10s %10s %10s %10s %10s\n" %("Peak Name","Amplitude","Esd(Amp.)","Rate","Esd(Rate)","StdErr"))
if (FIT180EXP == "Y"):
    fitlogfile=sys.argv[2]+"/single180exp.fit"
    fit180log = open(fitlogfile,'w')
    fit180log.write('# Id       Amplitude1    Amp1(esd)     Amplitude1    Amp1(esd)       R1       R1(esd)     Std \n')
if (FITBIEXP == "Y"):
    fitlogfile=sys.argv[2]+"/biexp.fit"
    fitbilog = open(fitlogfile,'w')
    fitbilog.write("# %10s%10s %10s  %10s %10s  %10s %10s  %10s %10s  %10s  %10s \n" %("Peak Name","C1","Esd(C1)","L1","Esd(L1)","C2","Esd(C2)","L2","Esd(L2)","Std","F(Single->Bi)") )

    #The shortened fit: Only use, the last 4 points.
    fitlogfile=sys.argv[2]+"/shortexp.fit"
    fitshortlog=open(fitlogfile,'w')
    fitshortlog.write('# Id       Amplitude    Amp(esd)        R1       R1(esd)     Std \n')
if ( LOGFILENR == "Y" ):
    fitlogfile=sys.argv[2]+".fit"
    fitlog = open(fitlogfile,'w')
    if ( zcoor == "2D" ):
        fitlog.write("#%13s%22s%22s%22s%22s%22s%5s\n" % (\
            "Id",\
            "---Omega1(ppm)---",\
            "---Omega2(ppm)---",\
            "----LW1(Hz)----",\
            "----LW2(Hz)----",\
            "---Intensity---",\
            "Isotopeshift"\
            ))
    elif ( zcoor == "3D" ):
        fitlog.write("#%13s%22s%22s%22s%22s%22s%22s%22s%5s\n" % (\
            "Id",\
            "---Omega1(ppm)---",\
            "---Omega2(ppm)---",\
            "---Omega3(ppm)---",\
            "----LW1(Hz)----",\
            "----LW2(Hz)----",\
            "----LW3(Hz)----",\
            "---Intensity---",\
            "Isotopeshift"\
            ))
    else:
        sys.stderr.write(" Urghhh .. please check your programming skills Hansen!\n")
        sys.stderr.write(" Internal Error #54\n")
        sys.exit(54)
#
# These two values specifies the dimensions of the data file and are convenient
# parameters to use later.
#########################################################################
NDIM=fudaIO.get_dim(fd)
FILESIZE=[]
for i in range(NDIM):
    FILESIZE.append(fudaIO.get_size(fd,i))    
if zcoor=="3D":
    FILESIZE=[FILESIZE[2],FILESIZE[1],FILESIZE[0]]
else:
    if zcoor=="2D":
        if(DIM==2):
            FILESIZE=[FILESIZE[1],FILESIZE[0]]  
        elif(DIM==3):
            FILESIZE=[FILESIZE[1],FILESIZE[0],1]
        else:
            sys.stderr.write(" Wrong number of dimensions!\n")
            sys.exit(10)
    else:
        FILESIZE=[FILESIZE[1],FILESIZE[0],FILESIZE[2]]   
#
# We number points in the file from 1 and up (Convention normally used for NMR spectroscopy)
dataIO.set_pt_offset(fd,1)
#
# Map the spectrum according to the nmrPipe header
fudaIO.map_default(fd,0,'ppm')
fudaIO.map_default(fd,1,'ppm')
#
# Map the planes from 1 --> FILESIZE[2] if arrayed 2D
if ( len(FILESIZE) > 2 ):
    if ( zcoor == "3D" ):
        fudaIO.map_default(fd,2,'ppm')
    elif(zcoor!="2D"):
        dataIO.set_xmap(fd, 2, 1, 1., FILESIZE[2], float(FILESIZE[2]))
    elif(zcoor=="2D"):
        dataIO.set_xmap(fd, 2, 1, 1., 2, 2.)
#
# Read in the peaklist into the array 'peakline'
################################################################
#
peakline=[]
#
try:
    peakfile=open(string.strip(peaklistfile),'r')
except (IOError):
    sys.stderr.write("\n Cannot open the peaklist file: %s \n" % (peaklistfile.strip()))
    sys.stderr.write(" PROGRAM ABORTED\n\n")
    sys.exit(1)
#    
errorfile=open(sys.argv[2]+"/error.log",'w')
#
for line in peakfile.readlines():
    if (len(string.split(line)) < 2 ):
        continue
    if ( string.split(line)[0] == "Assignment"):
        continue
    if ( string.split(line)[1] == "peaks" ):
        continue
    # The peakline structure is:
    # name,omega1,omega2,omega3,FIT_FLAG,OVERLAP_FLAG
    # If FIT_FLAG=0     : No fit has been performed
    # If FIT_FLAG=1     : The fit has been performed and converged
    # if FIT_FLAG=-1    : The fit has been performed and did not converge
    # if OVERLAP_FLAG=-1: No overlap occour
    # if OVERLAP_FLAG=i : The peak belong to overlap_peaks[i]
    linetosave=string.split(line)
    if ( linetosave[0][0] == "?" ):
        if ( verboselevel > 3 ):
            sys.stderr.write("\n Warning:\n")
            sys.stderr.write(" The following line in the peaklist might be corrupt\n")
            sys.stderr.write("-->%s<--\n" % ( line.strip() ))
        linetosave[0]=linetosave[0][2:len(linetosave[0])]
    linetosave.append(0)  #Have not fitted anything yet: set fit_flg=0
    overlap_flag=-1       #Default is no overlap 
    #
    #Check wether overlap, and store flag corresponding to the group:
    for i in range (0,len(overlap_peaks)):
        for k in range (0,len(overlap_peaks[i])):
            if ( linetosave[0]==overlap_peaks[i][k]):
                overlap_flag=i
    linetosave.append(overlap_flag)
    if ( len(FILESIZE) == 2 or len(FILESIZE) == 3 ):
        if ( zcoor == "3D" ):
            # Make checking for folding.
            peakline.append(linetosave)
        else:
            temp=[]
            temp.append(linetosave[0])  #Store the name
            temp.append(0.)
            #
            # Check for folding - should only occour in indirect dimension.
            #
            while ( float(linetosave[1]) > first[0] ):
                linetosave[1]=float(linetosave[1])-(first[0]-last[0])
            while ( float(linetosave[1]) < last[0] ):
                linetosave[1]=float(linetosave[1])+(first[0]-last[0])
            #
            if ( float(linetosave[2]) > first[1] or float(linetosave[2]) < last[1] ):
                sys.stderr.write("\n")
                sys.stderr.write(" Hmmm .. Maybe you should check your assignments!\n")
                sys.stderr.write(" You have assignned the peak %s to be at %.2fppm in the direct dimension,\n" % (linetosave[0],float(linetosave[2])))
                sys.stderr.write(" however, the spectrum only reach from %.2fppm to %.2fppm in this dimension\n"  % (first[1],last[1]))
                sys.stderr.write("                                                              PROGRAM ABORTED\n")
                sys.exit(20)
            #
            temp.append(linetosave[1])
            temp.append(linetosave[2])
            temp.append(linetosave[3])
            temp.append(linetosave[4])
            peakline.append(temp)
    else:
        sys.stderr.write("\n\n Stupid spectrum .. neither 2D nor 3D\n")
        sys.exit(117)
#
# Convert the 'name' in not_def_peaks to id.
for i in range(0,len(not_def_peak)):
    myint=peak_name_to_no(not_def_peak[i]['id'])
    if ( myint == -1 ):
        sys.stderr.write(" Cannot locate peak: %s in the peaklist\n" % (not_def_peak[i]['id']))
        sys.stderr.write(" Error occurs during reading the NOT_DEF_PEAK \n")
        sys.stderr.write("                              PROGRAM ABORTED \n")
        sys.exit(1)
    else:
        not_def_peak[i]['id']=myint         
#
# Define the array with chemical shifts values
omega=[0.,0.,0.]
#
#Define new function, if isotopeshift is present
# (Don't know if there is a not_def_peak line that has it, therefore, always define
# these functions - only takes a fraction of a second and little memory)
ftype_composite('gausslorej',1,1,'norm_gausslore','linear')
#
ftype_product('glore2d_Jpart',0,('gausslorej','gausslorej'))
ftype_product('glore3d_Jpart',0,('gausslorej','gausslorej','gausslorej'))
#
# Duplet
ftype_sum('glore2d_J',1,('glore2d_Jpart','glore2d_Jpart'))
ftype_sum('glore3d_J',1,('glore3d_Jpart','glore3d_Jpart'))
#
# Triplet
ftype_sum('glore2d_J-J',1,('glore2d_Jpart','glore2d_Jpart','glore2d_Jpart'))
ftype_sum('glore3d_J-J',1,('glore3d_Jpart','glore3d_Jpart','glore3d_Jpart'))
#
# 3D spectrum
if ( zcoor == "3D"):
    ftype_product('gausslore3d',1,('norm_gausslore','norm_gausslore','norm_gausslore'))
#
for j in range (0,len(peakline)):
    discard=0
    use_isotopeshift=0
    # Does the peak belong to the discard group 
    if ( discard_peaks_lg ):
        for peak in discard_peaks:
            if ( peakline[j][0] == peak ):
                discard=1
    #... or not in the include group?
    if ( include_peaks_lg ):
        discard=1
        for peak in include_peaks:
            if ( peakline[j][0] == peak ):
                discard=0
    if ( discard == 1 ):
        if ( verboselevel > 2 ):
            print 'Peak: ',peakline[j][0],' is discarded from the fit'
        continue
    if not (peakline[j][4]==0):
        if ( verboselevel > 3 ):
            print 'Peak: ',peakline[j][0],' has already been fitted'
    if (peakline[j][4] == 0):
        #Find number of peaks:
        if (peakline[j][5] == -1):
            nofpeaks=1  # Flag indicates no overlap
        else:
            nofpeaks=len(overlap_peaks[peakline[j][5]])
        if ( verboselevel > 1 ) :
            sys.stdout.write("\n\n Name of Group-of-peaks: %s \n" % (peakline[j][0],))
            sys.stdout.write(" Number of peaks in the group: %d\n" % (nofpeaks))
        #
        #Get the identity of the peaks:
        #
        #Here the array peak_fit_array is convenient to store information
        #Of the group to fit:
        #peak_fit_array=[id,lw1,lw2,lw3,rad1,rad2,rad3,isotopeshift,shape,intensity,restoreparameters]
        #Initialize the peak_fit_array[]
        #
        # As of March 19th shift to non-interger array
        # Flemming
        peak_fit_array=[]
        #
        if (nofpeaks == 1):
            if ( peak_name_to_no(peakline[j][0]) == -1 ):
                sys.stderr.write(" Cannot locate peak: %s in the peaklist\n" % (peakline[j][0]))
                sys.stderr.write(" Non-overlap fitting\n")
                sys.stderr.write(" PROGRAM ABORTED \n")
                sys.exit(1)
            #
            temp={}
            temp['id']  =peak_name_to_no(peakline[j][0])
            temp['lw1'] =def_linewidth[0]
            temp['lw2'] =def_linewidth[1]
            temp['lw3'] =def_linewidth[2]            
            temp['rad1']=def_radius[0]
            temp['rad2']=def_radius[1]
            temp['rad3']=def_radius[2]
            temp['shape']=shape
            temp['intensity']=1.
            temp['isotopeshift']=isotopeshift
            temp['restoreparameters']=0
            peak_fit_array.append(temp)
            temp={}
            
            """
            peak_fit_household(0,'id','save',peak_fit_array,peak_name_to_no(peakline[j][0]))
            peak_fit_household(0,'lw1','save',peak_fit_array,def_linewidth[0])
            peak_fit_household(0,'lw2','save',peak_fit_array,def_linewidth[1])
            peak_fit_household(0,'lw3','save',peak_fit_array,def_linewidth[2])            
            peak_fit_household(0,'rad1','save',peak_fit_array,def_radius[0])
            peak_fit_household(0,'rad2','save',peak_fit_array,def_radius[1])
            peak_fit_household(0,'rad3','save',peak_fit_array,def_radius[2])
            peak_fit_household(0,'isotopeshift','save',peak_fit_array,isotopeshift)
            peak_fit_household(0,'shape','save',peak_fit_array,shape)
            peak_fit_household(0,'intensity','save',peak_fit_array,1.)
            peak_fit_household(0,'restoreparameters','save',peak_fit_array,0)

            peak_fit_array.append([\
                peak_name_to_no(peakline[j][0]),\
                def_linewidth[0],\
                def_linewidth[1],\
                def_radius[0],\
                def_radius[1],\
                isotopeshift,\
                shape,\
                1.,\
                0])

            """
        #  
        if not (nofpeaks ==1):
            for i in range (0,nofpeaks):
                if ( peak_name_to_no(overlap_peaks[peakline[j][5]][i]) == -1 ):
                    sys.stderr.write(" Cannot locate peak: %s in the peaklist\n" % (overlap_peaks[peakline[j][5]][i]))
                    sys.stderr.write(" PROGRAM ABORTED \n")
                    sys.exit(1)
                #
                temp={}
                temp['id']  =peak_name_to_no(overlap_peaks[peakline[j][5]][i])
                temp['lw1'] =def_linewidth[0]
                temp['lw2'] =def_linewidth[1]
                temp['lw3'] =def_linewidth[2]            
                temp['rad1']=def_radius[0]
                temp['rad2']=def_radius[1]
                temp['rad3']=def_radius[2]
                temp['shape']=shape
                temp['intensity']=1.
                temp['isotopeshift']=isotopeshift
                temp['restoreparameters']=0
                peak_fit_array.append(temp)
                temp={}

                """
                peak_fit_array.append([\
                    peak_name_to_no(overlap_peaks[peakline[j][5]][i]),\
                    def_linewidth[0],\
                    def_linewidth[1],\
                    def_radius[0],\
                    def_radius[1],\
                    isotopeshift,\
                    shape,\
                    1.,\
                    0])
                """
        #
        #Get ready to fit new set of peaks
        data_del_all()
        func_del_all()
        param_del_all()
        dtype_del_all()
        #
        # Explanatory parameters.
        param('f1',kind='EXPL')
        param('f2',kind='EXPL')
        param('f3',kind='EXPL')
        #
        # Auxiliary parameters
        param('Zero' ,value= 0.00,kind='CONST')
        param('Half' ,value= 0.50,kind='CONST')
        param('-Half',value=-0.50,kind='CONST')
        param('One'  ,value= 1.00,kind='CONST')
        param('-One' ,value=-1.00,kind='CONST')
        param('Two'  ,value= 2.00,kind='CONST')
        #
        #Check if the peak is in the array 'not_def_peak'
        for i in range (0,len(not_def_peak)):
            for k in range (0,len(peak_fit_array)):
                #if (not_def_peak[i][0] == peak_fit_array[k][0]):
                if (not_def_peak[i]['id'] == peak_fit_array[k]['id']):
                    #if not ( not_def_peak[i][8] == 0 ):
                    if not ( not_def_peak[i]['restoreparameters'] == 0 ):
                        if ( verboselevel > 2 ):
                            sys.stdout.write(' Uses non-default start values for: %s\n' %(peakline[not_def_peak[i]['id']][0]))
                            sys.stdout.write(' Start parameters will be read from parameter-file\n')
                            sys.stdout.write(' SHAPE=      %10s\n' %( not_def_peak[i]['shape']))
                    else:
                        if ( verboselevel > 2 ):
                            sys.stdout.write(' Uses non-default start values for: %s\n' %(peakline[not_def_peak[i]['id']][0]))
                            sys.stdout.write(' LINEWIDTH_F1=%10.4f\n' %( not_def_peak[i]['lw1'] ))
                            sys.stdout.write(' LINEWIDTH_F2=%10.4f\n' %( not_def_peak[i]['lw2'] ))
                            if ( zcoor == "3D" ):
                                sys.stdout.write(' LINEWIDTH_F3=%10.4f\n' %( not_def_peak[i]['lw3'] ))
                            sys.stdout.write(' RADIUS_F1=   %10.4f\n' %( not_def_peak[i]['rad1'] ))
                            sys.stdout.write(' RADIUS_F2=   %10.4f\n' %( not_def_peak[i]['rad2'] ))
                            if ( zcoor == "3D" ):
                                sys.stdout.write(' RADIUS_F3=   %10.4f\n' %( not_def_peak[i]['rad3'] ))
                            print ' ISOTOPESHIFT=   ',not_def_peak[i]['isotopeshift']
                            sys.stdout.write(' SHAPE=       %10s\n' %( not_def_peak[i]['shape']))
                            sys.stdout.write(' INTENSITY=   %10.4f\n' %( not_def_peak[i]['intensity']))
                    #
                    for sommer in ['id','lw1','lw2','lw3','rad1','rad2','rad3','isotopeshift','shape','intensity','restoreparameters']:
                        peak_fit_array[k][sommer]=not_def_peak[i][sommer]

        #
        for i in range (0,len(peak_fit_array)):
            #print peak_fit_array[i]['shape']
            if ( peak_fit_array[i]['shape'] == "VOLUME" ):
                #No parameters are defined.
                if ( zcoor=="2D" ):
                    sys.stderr.write(" Volume fitting in 2D spectra is not implemented yet!\n")
                    sys.exit(1)
                continue
            # Read parameters from parameterfile 
            rp=0 
            if not ( peak_fit_array[i]['restoreparameters'] == 0 ):
                rp=1
                rpfilename=''
                if ( peak_fit_array[i]['restoreparameters'] == 1 ):
                    rpfilename=sys.argv[2]+"/"+peakline[peak_fit_array[i]['id']][0]+".par" 
                else:
                    rpfilename=peak_fit_array[i]['restoreparameters']
                try:
                    rpfile=open(rpfilename,'r')
                except (IOError):
                    sys.stderr.write("\n The parameterfile \'%s\' cannot be opened.\n"  % (rpfilename ))
                    sys.stderr.write(" PROGRAM ABORTED \n")
                    sys.exit(1)
                    
                rpvals=[]
                rpnames=[]
                for plines in rpfile.readlines():
                    its=string.split(plines)
                    if ( len(its) > 1 and not its[0][0] == '#' ):
                        rpvals.append(float(its[1]))
                        rpnames.append(its[0])
                rpfile.close()
                #
                # Save the parameters in the intensity array, not to be used anyway
                peak_fit_array[i]['intensity']=[rpnames,rpvals]
            #
            # Freq 01
            if ( rp ):
                try:
                    myint=rpnames.index('f01')
                except (IndexError,ValueError):
                    sys.stderr.write("Parameter \'f01\' is not found in the parameter file\n")
                    sys.exit(10)
                #
                param('f01_'+str(i),value=rpvals[myint])
            else:
                # Take value from the peaklist
                if ( zcoor=="3D"):
                    param('f01_'+str(i),value=float(peakline[peak_fit_array[i]['id']][1]))
                else:
                    param('f01_'+str(i),value=float(peakline[peak_fit_array[i]['id']][2]))
            #
            # linewidth w1
            if ( rp ):
                try:
                    myint=rpnames.index('w1')
                except (IndexError,ValueError):
                    sys.stderr.write("Parameter \'w1\' is not found in the parameter file\n")
                    sys.exit(10)
                #
                param('w1_'+str(i),value=rpvals[myint])
            else:
                param('w1_'+str(i),value=float(peak_fit_array[i]['lw1']))
            #
            # isotopeshift
            if not ( peak_fit_array[i]['isotopeshift']=='n'):
                read_j1=1
                #if ( len(peak_fit_array[i]['isotopeshift']) > 2 ):
                #    if ( peak_fit_array[i]['isotopeshift'][2] == "fix" ):
                #        read_j1=0
                if ( "fix" in peak_fit_array[i]['isotopeshift'] ):
                    read_j1=0
                #
                if ( rp and read_j1 ):
                    try:
                        myint=rpnames.index('j1')
                    except (IndexError,ValueError):
                        sys.stderr.write("Parameter \'j1\' is not found in the parameter file\n")
                        sys.exit(1)
                    tempval=rpvals[myint]
                else:
                    tempval=peak_fit_array[i]['isotopeshift'][0]/sfrq[0]
                #
                #if ( len(peak_fit_array[i]['isotopeshift']) > 2 ):
                #    if ( peak_fit_array[i]['isotopeshift'][2] == "fix" ):
                #        param('j1_'+str(i),value=tempval,kind='CONST')
                #    elif ( peak_fit_array[i]['isotopeshift'][2] == "triplet" ):
                #        param('j1_'+str(i),value=tempval)
                #        if ( verboselevel > 3 ):
                #            sys.stdout.write("Defined triplet for peak: %s\n" %(peakline[peak_fit_array[i]['id']][0]))
                #    else:
                #        sys.stderr.write(" Isotopeshift parameter %s is unknown\n" % (peak_fit_array[i]['isotopeshift'][2],))
                #        sys.stderr.write("                                PROGRAM ABORTED\n")
                #        sys.exit(1)
                if ( "fix" in peak_fit_array[i]['isotopeshift'] ):
                    param('j1_'+str(i),value=tempval,kind='CONST')
                else:
                    param('j1_'+str(i),value=tempval)
            if ( peak_fit_array[i]['shape']  == "LORENTZIAN" ):
                param('g1_'+str(i),value=0, kind='CONST' )
                #param('p1_'+str(i),value=0 )
            elif ( peak_fit_array[i]['shape'] == "GAUSSIAN" ):
                param('g1_'+str(i),value=1, kind='CONST' )
            elif ( peak_fit_array[i]['shape'] == "GLORE" ):
                if ( rp ):
                    try:
                        myint=rpnames.index('g1')
                    except (IndexError,ValueError):
                        sys.stderr.write("Parameter \'g1\' is not found in the parameter file\n")
                        sys.exit(10)
                    param('g1_'+str(i),value=rpvals[myint])
                else:
                    param('g1_'+str(i),value=0.5)
            else:
                sys.stderr.write(" The shape: \'%s\' is unknown\n" % ( peak_fit_array[i]['shape'])) 
                sys.stderr.write("                      PROGRAM ABORTED\n")
                sys.exit(1)
                
            #
            # Frequencies in the second dimension
            if ( rp ):
                try:
                    myint=rpnames.index('f02')
                except (IndexError,ValueError):
                    sys.stderr.write("Parameter \'f02\' is not found in the parameter file\n")
                    sys.exit(10)
                #
                param('f02_'+str(i),value=rpvals[myint])
            else:                
                if ( zcoor == "3D" ):
                    param('f02_'+str(i),value=float(peakline[peak_fit_array[i]['id']][2]))
                else:
                    param('f02_'+str(i),value=float(peakline[peak_fit_array[i]['id']][3]))
            #
            # W2 linewidth
            if ( rp ):
                try:
                    myint=rpnames.index('w2')
                except (IndexError,ValueError):
                    sys.stderr.write("Parameter \'w2\' is not found in the parameter file\n")
                    sys.exit(10)
                #
                param('w2_'+str(i),value=rpvals[myint])
            else:
                param('w2_'+str(i),value=float(peak_fit_array[i]['lw2']))
            #
            # Check for isotopeshift
            if not ( peak_fit_array[i]['isotopeshift']=='n'):
                read_j2=1
                if ( "fix" in peak_fit_array[i]['isotopeshift']):
                    read_j2=0
                #if ( len(peak_fit_array[i]['isotopeshift']) > 2 ):
                #    if ( peak_fit_array[i]['isotopeshift'][2] == "fix" ):
                #        read_j2=0
                #
                if ( rp and read_j2 ):
                    try:
                        myint=rpnames.index('j2')
                    except (IndexError,ValueError):
                        sys.stderr.write("Parameter \'j2\' is not found in the parameter file\n")
                        sys.exit(1)              
                    tempval=rpvals[myint]
                else:
                    tempval=peak_fit_array[i]['isotopeshift'][1]/sfrq[1]
                #
                #if ( len(peak_fit_array[i]['isotopeshift']) > 2 ):
                #    if ( peak_fit_array[i]['isotopeshift'][2] == 'fix' ):
                #        param('j2_'+str(i),value=tempval, kind='CONST')
                #    elif ( peak_fit_array[i]['isotopeshift'][2] == 'triplet' ):
                #        param('j2_'+str(i),value=tempval)
                #    else:
                #        sys.stderr.write(" The isotopeshift parameter %s is not regonized!\n PROGARM ABORTED\n" % (  peak_fit_array[i]['isotopeshift'][2] ))
                #        sys.exit(2)
                if ( "fix" in peak_fit_array[i]['isotopeshift'] ):
                    param('j2_'+str(i),value=tempval, kind='CONST')
                else:
                    param('j2_'+str(i),value=tempval)            
            if ( peak_fit_array[i]['shape'] == "LORENTZIAN" ):
                param('g2_'+str(i),value=0, kind='CONST' )
            elif ( peak_fit_array[i]['shape'] == "GAUSSIAN" ):
                param('g2_'+str(i),value=1, kind='CONST' )                
            elif ( peak_fit_array[i]['shape'] == "GLORE" ):
               if ( rp):
                    try:
                        myint=rpnames.index('g2')
                    except (IndexError,ValueError):
                        sys.stderr.write("Parameter \'g2\' is not found in the parameter file\n")
                        sys.exit(1)              
                    param('g2_'+str(i),value=rpvals[myint] )
               else:
                   param('g2_'+str(i),value=0.5 )
            else:
                sys.stderr.write(" The shape \'%s\' is unknown\n" % (peak_fit_array[i]['shape']))
                sys.exit(1)
            #
            #Check for isotopeshift, and if so define the ratio
            if ( not peak_fit_array[i]['isotopeshift']=='n' ):
                if ( rp):
                    try:
                        myint=rpnames.index('Ratio')
                    except (IndexError,ValueError):
                        sys.stderr.write("Parameter \'Ratio\' is not found in the parameter file\n")
                        sys.exit(1)
                    param('Ratio_'+str(i),value=rpvals[myint])
                else:
                    param('Ratio_'+str(i),value=1.0)
            #
            #
            if ( zcoor == "3D" ):
                #
                # Start with defining frequency in the 3rd dimension:
                if ( rp ):
                    try:
                        myint=rpnames.index('f03')
                    except (IndexError,ValueError):
                        sys.stderr.write("Parameter \'f03\' is not found in the parameter file\n")
                        sys.exit(1)
                    #
                    param('f03_'+str(i),value=rpvals[myint])
                else:
                    # Take value from the peaklist
                    param('f03_'+str(i),value=float(peakline[peak_fit_array[i]['id']][3]))
                #
                # linewidth w3
                if ( rp ):
                    try:
                        myint=rpnames.index('w3')
                    except (IndexError,ValueError):
                        sys.stderr.write("Parameter \'w3\' is not found in the parameter file\n")
                        sys.exit(10)
                        #
                    param('w3_'+str(i),value=rpvals[myint])
                else:
                    param('w3_'+str(i),value=float(peak_fit_array[i]['lw3']))
                #
                # isotopeshift
                if not ( peak_fit_array[i]['isotopeshift']=='n'):
                    read_j3=1
                    if ( "fix" in peak_fit_array[i]['isotopeshift']):
                        read_j3=0
                    if ( rp and read_j3 ):
                        try:
                            myint=rpnames.index('j3')
                        except (IndexError,ValueError):
                            sys.stderr.write("Parameter \'j3\' is not found in the parameter file\n")
                            sys.exit(1)              
                        tempval=rpvals[myint]
                    else:
                        tempval=peak_fit_array[i]['isotopeshift'][2]/sfrq[2]
                    if ( "fix" in peak_fit_array[i]['isotopeshift'] ):
                        param('j3_'+str(i),value=tempval, kind='CONST')
                    else:
                        param('j3_'+str(i),value=tempval)            
                if ( peak_fit_array[i]['shape']  == "LORENTZIAN" ):
                    param('g3_'+str(i),value=0, kind='CONST' )
                elif ( peak_fit_array[i]['shape'] == "GAUSSIAN" ):
                    param('g3_'+str(i),value=1, kind='CONST' )
                elif ( peak_fit_array[i]['shape'] == "GLORE" ):
                    if ( rp ):
                        try:
                            myint=rpnames.index('g3')
                        except (IndexError,ValueError):
                            sys.stderr.write("Parameter \'g3\' is not found in the parameter file\n")
                            sys.exit(1)
                        param('g3_'+str(i),value=rpvals[myint])
                    else:
                        param('g3_'+str(i),value=0.5)
                else:
                    sys.stderr.write(" The shape: \'%s\' is unknown\n" % ( peak_fit_array[i]['shape'])) 
                    sys.stderr.write("                      PROGRAM ABORTED\n")
                    sys.exit(1)
        #
        #         
        # Do the mapping of the spectrum
        centrum=[]
        radius=[]
        for i in range (0,len(peak_fit_array)):
            if not zcoor=="3D":
                centrum.append([\
                   float(peakline[peak_fit_array[i]['id']][2]),\
                   float(peakline[peak_fit_array[i]['id']][3])]\
                )
            else:
                centrum.append([\
                   float(peakline[peak_fit_array[i]['id']][1]),\
                   float(peakline[peak_fit_array[i]['id']][2]),\
                   float(peakline[peak_fit_array[i]['id']][3])]\
                )
            #
            #Check for isotopeshift
            if not ( peak_fit_array[i]['isotopeshift']=='n'):
                if not zcoor=="3D":
                    centrum.append([\
                        float(peakline[peak_fit_array[i]['id']][2])+peak_fit_array[i]['isotopeshift'][0]/sfrq[0],\
                        float(peakline[peak_fit_array[i]['id']][3])+peak_fit_array[i]['isotopeshift'][1]/sfrq[1]]\
                    )
                else:
                    centrum.append([\
                        float(peakline[peak_fit_array[i]['id']][1])+peak_fit_array[i]['isotopeshift'][0]/sfrq[0],\
                        float(peakline[peak_fit_array[i]['id']][2])+peak_fit_array[i]['isotopeshift'][1]/sfrq[1],\
                        float(peakline[peak_fit_array[i]['id']][3])+peak_fit_array[i]['isotopeshift'][2]/sfrq[2]]\
                    )
                #
                if ( "triplet" in peak_fit_array[i]['isotopeshift']):
                    centrum.append([\
                          float(peakline[peak_fit_array[i]['id']][2])+2*peak_fit_array[i]['isotopeshift'][0]/sfrq[0],\
                          float(peakline[peak_fit_array[i]['id']][3])+2*peak_fit_array[i]['isotopeshift'][1]/sfrq[1]]\
                    )
            #
            if not zcoor=="3D":
                radius.append([\
                   int(float(peak_fit_array[i]['rad1'])*FILESIZE[0]/(first[0]-last[0])),\
                   int(float(peak_fit_array[i]['rad2'])*FILESIZE[1]/(first[1]-last[1]))]\
                   )
            else:
                radius.append([\
                   int(float(peak_fit_array[i]['rad1'])*FILESIZE[0]/(first[0]-last[0])),\
                   int(float(peak_fit_array[i]['rad2'])*FILESIZE[1]/(first[1]-last[1])),\
                   int(float(peak_fit_array[i]['rad3'])*FILESIZE[2]/(first[2]-last[2]))],\
                   )
            #
            #Check for isotopeshift
            if not ( peak_fit_array[i]['isotopeshift']=='n'):
                radius.append(radius[len(radius)-1])
                if ( "triplet" in peak_fit_array[i]['isotopeshift'] ):
                    radius.append(radius[len(radius)-1])
        #
        if ( verboselevel > 2 ):
            sys.stdout.write(" Radius  (points) : ")
            PrintArray(sys.stdout,radius,"%7d")
            sys.stdout.write(" Centrum    (ppm) : ")
            PrintArray(sys.stdout,centrum,"%7.3f")
            sys.stdout.write(" Freq.      (MHz) : ")
            PrintArray(sys.stdout,sfrq,"%.2f")
            sys.stdout.write(" IsotopeShift     : ")
            for i in range(len(peak_fit_array)):
                print "\t",peak_fit_array[i]['isotopeshift'],
            sys.stdout.write("\n")
        NoSlices = 0
        if ( zcoor == "2D" ):
            NoSlices = 1
        elif ( zcoor == "3D" ):
            # Warning.. NoSlices will be determined later.
            NoSlices=2*radius[0][0]+1
        else:
            NoSlices=FILESIZE[2]-len(discard_slices)
            SlicesToUse=[]
            for i in range(FILESIZE[2]):
                discard_slice=0
                for slice in discard_slices:
                    if ( int(slice) == i+1 ):
                        discard_slice=1
                if not ( discard_slice ):
                    SlicesToUse.append(i+1)
        #
        # 3D spectrum
        if ( zcoor == "3D" ):
            #
            Iname=[]
            for k in range (0,len(peak_fit_array)):
                if ( peak_fit_array[k]['shape']=="VOLUME" ):
                    continue
                #
                # Intensity
                Iname.append('I_%d' % (k) )
                #
                # check if we are reading in the intensities
                if not ( peak_fit_array[k]['restoreparameters'] == 0 ):
                    ThisIname="I"
                    rpnames=peak_fit_array[k]['intensity'][0]
                    rpvals=peak_fit_array[k]['intensity'][1]
                    try:
                        myint=rpnames.index(ThisIname)
                    except (IndexError, ValueError):
                        sys.stderr.write(" Parameter %s is not found in the parameter file\n" % (ThisIname))
                        sys.stderr.write("                                 PROGRAM ABORTED\n")
                        sys.exit(11)
                    param(Iname[k], value=rpvals[myint])
                else:
                    # The intensity slot has 'hopefully' not been overwritten
                    StartIntensity=3*noise*peak_fit_array[k]['intensity']
                    param(Iname[k], value=StartIntensity)
            #
            # Only if baseline
            if ( baseline == 'Y' or baseline == 'y' ):
                Bname = "B"
                param(Bname, value=noise)
            #
            # Dtypes
            dtname = 'dt'
            if not (dtype_exists(dtname)):
                dtype(dtname,('f3','f2','f1'))
            #
            # Function names
            Fname=[]
            for k in range (0,len(peak_fit_array)):
                Fname.append('peak_%d' % (k))
            Fname2 = 'base'
            #
            #The functions
            for k in range (0,len(peak_fit_array)):
                if ( peak_fit_array[k]['isotopeshift']=='n' ):
                    func(Fname[k],'gausslore3d',dtname,\
                         (\
                           Iname[k],
                           'f1','f01_'+str(k),'w1_'+str(k),'g1_'+str(k),\
                           'f2','f02_'+str(k),'w2_'+str(k),'g2_'+str(k),\
                           'f3','f03_'+str(k),'w3_'+str(k),'g3_'+str(k),\
                         )\
                     )
                else:
                    if not ( "triplet" in peak_fit_array[k]['isotopeshift'] ):
                        func(Fname[k],'glore3d_J',dtname,\
                             ( Iname[k],\
                                 'Ratio_'+str(k),'f1', 'One','j1_'+str(k),'f01_'+str(k),'w1_'+str(k),'g1_'+str(k), \
                                           'One','f2', 'One','j2_'+str(k),'f02_'+str(k),'w2_'+str(k),'g2_'+str(k), \
                                           'One','f3', 'One','j3_'+str(k),'f03_'+str(k),'w3_'+str(k),'g3_'+str(k), \
                                           'One','f1','Zero','j1_'+str(k),'f01_'+str(k),'w1_'+str(k),'g1_'+str(k), \
                                           'One','f2','Zero','j2_'+str(k),'f02_'+str(k),'w2_'+str(k),'g2_'+str(k), \
                                           'One','f3','Zero','j3_'+str(k),'f03_'+str(k),'w3_'+str(k),'g3_'+str(k)  \
                               ))
                    else:
                        sys.stderr.write(" Isotopeshift(triplet) not implemented for 3D spectral fitting yet\n")
                        sys.stderr.write("                                                    PROGRAM ABOTED\n")
                        sys.exit(-1)
            #
            #The base level
            if ( baseline == 'Y' or baseline == 'y' ):
                func(Fname2,'constant',dtname,(Bname,))
                if ( verboselevel > 4 and i == 0 ):
                    sys.stderr.write('Baseline added\n')
            #fudaIO_read_xregion(......shape='eer')
            # Read data region.
            for k in range(0,len(radius)):
                    fudaIO.data_read_xcenter(fd,\
                          (centrum[k][2],centrum[k][1],centrum[k][0]),\
                          ( radius[k][2], radius[k][1], radius[k][0]), \
                          noise,shape='eee')
            dtype_set_purge(dtname,1)
        #
        # Loop over the slices in the arraied 2D spectrum
        for i in range(1,NoSlices+1+len(discard_slices)):
            if ( zcoor=="3D" ):
                break
            #
            # Check if this slice should be excluded:
            discard_slice=0
            for slice in discard_slices:
                if ( int(slice) == i ):
                    discard_slice=1
                    if ( verboselevel > 3 ):
                        print ' Slice number ',i,' is discared from the fit.'
            if ( discard_slice ):
                continue
            #
            # Parameters (intensities+baseline)
            Iname=[]
            for k in range (0,len(peak_fit_array)):
                if( peak_fit_array[k]['shape'] == "VOLUME" ):
                    continue
                Iname.append('I%d_%d' % (i,k) )
                #
                # check if we are reading in the intensities
                if not ( peak_fit_array[k]['restoreparameters'] == 0 ):
                    ThisIname="I%d" % (i)
                    rpnames=peak_fit_array[k]['intensity'][0]
                    rpvals=peak_fit_array[k]['intensity'][1]
                    try:
                        myint=rpnames.index(ThisIname)
                    except (IndexError, ValueError):
                        sys.stderr.write(" Parameter %s is not found in the parameter file\n" % (ThisIname))
                        sys.stderr.write("                                 PROGRAM ABORTED\n")
                        sys.exit(10)
                    param(Iname[k], value=rpvals[myint])
                else:
                    # The intensity slot has 'hopefully' not been overwritten
                    StartIntensity=3*noise*peak_fit_array[k]['intensity']
                    param(Iname[k], value=StartIntensity)
            #
            # Only if baseline
            if ( baseline == 'Y' or baseline == 'y' ):
                Bname = 'B%d' % (i,)
                param(Bname, value=noise)
            #
            # Dtypes
            dtname = 'dt%d' % (i,)
            if not (dtype_exists(dtname)):
                if ( len(FILESIZE) == 2 ):
                    dtype(dtname,('f2','f1'))
                elif ( len(FILESIZE) == 3 ):
                    dtype(dtname,('f2','f1','f3'))
                else:
                    sys.stderr.write(" hmmm .. strange spectrum .. \n This error should have been caught earlier in the program\n INTERNAL ERROR #13 --- PROGRAM ABORTED")                    
                    sys.exit(13)
            #
            # Function names
            Fname=[]
            for k in range (0,len(peak_fit_array)):
                Fname.append('peak%d_%d' % (i,k))
            Fname2 = 'base%d' % (i,)
            #
            #The functions
            for k in range (0,len(peak_fit_array)):
                if ( peak_fit_array[k]['shape']=="VOLUME" ):
                    if k>0:
                        sys.stderr.write(" Cannot define overlapped peaks, while estimating Volume.\n")
                        sys.stderr.write("                                         PROGRAM ABORTED\n")
                        sys.exit(1)
                    FisseName= "Sommer%d" % (i,)
                    DummyName= "Frisk%d" % (i,)
                    param(DummyName,value=0.)                    
                    func(FisseName,'constant',dtname,(DummyName,))
                    continue
                #
                # Check for isotopeshift
                if not ( peak_fit_array[k]['isotopeshift']=='n'):
                    use_isotopeshift=1
                    #
                    Multiplicity=2
                    if (  len(peak_fit_array[k]['isotopeshift']) > 2 ):
                        if (  peak_fit_array[k]['isotopeshift'][2]=='triplet' ):
                            Multiplicity=3
                    if ( Multiplicity==3 ):
                        func(Fname[k],'glore2d_J-J',dtname,\
                             ( Iname[k],\
                                 'Ratio_'+str(k),'f1', 'Two','j1_'+str(k),'f01_'+str(k),'w1_'+str(k),'g1_'+str(k),\
                                 'Ratio_'+str(k),'f2', 'Two','j2_'+str(k),'f02_'+str(k),'w2_'+str(k),'g2_'+str(k),\
                                            \
                                 'Ratio_'+str(k),'f2', 'One','j2_'+str(k),'f02_'+str(k),'w2_'+str(k),'g2_'+str(k),\
                                           'Two','f1', 'One','j1_'+str(k),'f01_'+str(k),'w1_'+str(k),'g1_'+str(k),\
                                            \
                                           'One','f2','Zero','j2_'+str(k),'f02_'+str(k),'w2_'+str(k),'g2_'+str(k),\
                                           'One','f1','Zero','j1_'+str(k),'f01_'+str(k),'w1_'+str(k),'g1_'+str(k) \
                               ))
                    else:
                        func(Fname[k],'glore2d_J',dtname,\
                             ( Iname[k],\
                                 'Ratio_'+str(k),'f1', 'One','j1_'+str(k),'f01_'+str(k),'w1_'+str(k),'g1_'+str(k),\
                                           'One','f2', 'One','j2_'+str(k),'f02_'+str(k),'w2_'+str(k),'g2_'+str(k),\
                                           'One','f2','Zero','j2_'+str(k),'f02_'+str(k),'w2_'+str(k),'g2_'+str(k),\
                                           'One','f1','Zero','j1_'+str(k),'f01_'+str(k),'w1_'+str(k),'g1_'+str(k) \
                               ))
                else:
                    func(Fname[k],'gausslore2d',dtname,\
                         (\
                           Iname[k],
                           'f2','f02_'+str(k),'w2_'+str(k),'g2_'+str(k),\
                           'f1','f01_'+str(k),'w1_'+str(k),'g1_'+str(k)\
                         )\
                     )
            #
            #The base level
            if ( baseline == 'Y' or baseline == 'y' ):
                func(Fname2,'constant',dtname,(Bname,))
                if ( verboselevel > 4 and i == 0 ):
                    sys.stderr.write('Baseline added\n')
            #fudaIO_read_xregion(......shape='eer')
            # Read data region.
            for k in range(0,len(radius)):
                if ( len(FILESIZE)==3):
                    fudaIO.data_read_xcenter(fd,(centrum[k][1],centrum[k][0],float(i)),(radius[k][1],radius[k][0],0), noise,shape='eer')
                elif ( len(FILESIZE) == 2 ):
                    fudaIO.data_read_xcenter(fd,(centrum[k][1],centrum[k][0]),(radius[k][1],radius[k][0]),noise,shape='ee')
                else:
                    print 'boehh .. ' 
                    sys.exit(165)
            dtype_set_purge(dtname,1)
        
        # Minimize.
        eval_init()
        if ( verboselevel < 5 ):
            lm(nprint=0)
        lm(tol=lmtol,maxfev=lmmaxfev)
        #
        # Compensate for error in derivatives in fuda.        
        if ( ( len(FILESIZE) == 2 ) and use_isotopeshift ):
            if ( verboselevel > 2 ):
                sys.stderr.write("Warning: Using numerical derivatives!\n")
            lm(numderiv=1)
        else:
            lm(numderiv=0)
        # Print parameters.
        #for P in range(eval_get('nfree')):
        #    print eval_get_free(P),param_get(eval_get_free(P),'value')
        #sys.exit(10)
        #
        #
        # If there are more than one peak in the group,
        # then make a more sufisticated fitting procedure
        if ( len(peak_fit_array) > 1 ):
            for k in range(len(peak_fit_array)):
                if peak_fit_array[k]['shape']=="VOLUME":
                    sys.stderr.write(" Cannot determine VOLUME, when more than one peak is present\n")
                    sys.stderr.write("                                             PROGRAM ABORTED\n")
                    sys.exit(1)            
            if ( verboselevel > 4 ):
                sys.stdout.write("\n *** Step One ***\n  Fitting intensities and linewidths,\n  fixing other spectral parameters.\n\n")
            #
            FitIsOK=1
            #
            SParam=[]
            for P in range(eval_get('nfree')):
                Name=eval_get_free(P)
                if not ( ( Name[0] == "I" ) or (Name[0] == "w" ) or ( Name[0] == "R" ) ):
                    SParam.append(Name)
                    param(Name,free=0)
            #
            lm_minimize()
            if not ( lm_get('fit_converged') ):       
                FitIsOK=0
            #
            if ( FitIsOK ):
                if ( verboselevel > 4 ):
                    lm_report()
                lm_update_param()
                #
                # Invert the fixed values:
                #
                if ( verboselevel > 4 ):
                    sys.stdout.write("\n *** Step Two *** \n  Fixing intensities and linewidths\n  fitting other spectral parameters\n\n")
                for P in SParam:
                    param(P,free=1)
                #
                IParam=[]
                for P in range(eval_get('nfree')):
                    Name=eval_get_free(P)
                    if ( ( Name[0] == "I" ) or (Name[0]== "w" ) or ( Name[0] == "R" ) ) :
                        IParam.append(Name)
                        param(Name,free=0)      
                #
                lm_minimize()
                if not ( lm_get('fit_converged') ):       
                    FitIsOK=0
            if ( FitIsOK ):
                #
                if ( verboselevel > 4 ):
                    lm_report()
                lm_update_param()
                #
                # Now Free all parameters
                for P in IParam:
                    param(P,free=1)
                #
                if ( verboselevel > 4 ):
                    sys.stdout.write("\n *** Final Step *** \n  Free all parameters\n\n")
                #
                lm_minimize()
        else:
            if ( peak_fit_array[0]['shape'] == "VOLUME" ):
                if ( zcoor == "3D" ):
                    sys.stderr.write(" Sorry! Volume determination of 3D spectra is not implemented yet... \n")
                    sys.stderr.write(" Go an take a coffee and pray to Santa Claus!\n")
                    sys.stderr.write("                                                     PROGRAM ABORTED\n")
                Volumen=[]
                EsdVolumen=[]
                for s in range(NoSlices):
                    Volumen.append(0.)
                    EsdVolumen.append(0.)
                    for p in range(eval_get('ndata')/NoSlices):
                        if not ( zcoor == "2D" ):
                            Volumen[s]+=eval_get_data(s*eval_get('ndata')/NoSlices+p)[3]
                            EsdVolumen[s]+=eval_get_data(s*eval_get('ndata')/NoSlices+p)[4]
                        else:
                            Volumen[s]+=eval_get_data(s*eval_get('ndata')/NoSlices+p)[2]
                            EsdVolumen[s]+=eval_get_data(s*eval_get('ndata')/NoSlices+p)[3]
            else:
                lm_minimize()
        #
        #
        if ( verboselevel > 3 ):
            if not ( peak_fit_array[0]['shape'] == "VOLUME" ):
                lm_report()
        if ( verboselevel > 1 ):
            if ( lm_get('fit_converged') ):
                sys.stdout.write(" ***\n *** The fit converged successfully :) \n ***\n")
            else:
                sys.stdout.write(" ***\n *** ... NOT CONVERGED! :( \n ***\n")
            sys.stdout.write(" X2 of the fit: %12.5e\n" % ( math.pow(lm_get('enorm'),2.)))
            sys.stdout.write(" Nfree        : %12d\n" % ( eval_get('nfree')))
            sys.stdout.write(" Ndata        : %12d\n" % ( eval_get('ndata')))
        #
        # Print data vs fit, if asked for:
        if ( printdata=='y' or printdata=='Y'):
            datafile=sys.argv[2]+"/"+peakline[peak_fit_array[0]['id']][0]+".dat"
            dataout = open(datafile,'w')
            ndata=eval_get('ndata')
            #
            # print out the fitting statistics:
            if ( lm_get('fit_converged') ):
                dataout.write("# ***\n# *** The fit converged successfully \n# ***\n")
            else:
                dataout.write("# ***\n# *** The fit did not converge \n# ***\n")
            dataout.write("# X2 of the fit  : %12.5e\n" % ( math.pow(lm_get('enorm'),2.)))
            dataout.write("# Free Parameters: %12d\n" % ( eval_get('nfree')))
            dataout.write("# Data Points    : %12d\n" % ( eval_get('ndata')))
            # Debugging !
            #ExpData=[]
            #FitData=[]
            #for sommer in range(ndata):
            #    ExpData.append(eval_get_data(sommer)[2])
            #lm_update_param()
            #eval_data_recalc()
            #for sommer in range(ndata):
            #    FitData.append(eval_get_data(sommer)[2])
            #for sommer in range(ndata):
            #    print eval_get_data(sommer)[0],eval_get_data(sommer)[1],ExpData[sommer],FitData[sommer]
            #sys.exit(1)
            SpecData=[]
            width=0
            length=0
            maxwidth=0
            for datacount in range(ndata):
                if ( datacount+1 < ndata ):
                    if ( math.fabs(eval_get_data(datacount)[1] - eval_get_data(datacount+1)[1]) < 1e-10 ):
                        width=width+1
                    else:
                        if ( width > maxwidth ):
                            maxwidth=width
                        width=0.
                        length=length+1
                SpecData.append(eval_get_data(datacount)[NDIM])
            #
            #Calculate data from parameters
            if ( lm_get('fit_converged') ):
                lm_update_param()
                eval_data_recalc()
                if ( zcoor=="3D" ):
                    dataout.write('#%10s %11s %11s %11s %11s \n' % ('F3(ppm)','F2(ppm)','F1(ppm)','Data','Calc'))
                else:
                    dataout.write('#%10s %11s %11s %11s \n' % ('F2(ppm)','F1(ppm)','Data','Calc'))
            else:
                if ( zcoor=="3D"):
                    dataout.write('#%10s %11s %11s %11s %11s \n' % ('F3(ppm)','F2(ppm)','F1(ppm)','Data','Init'))
                else:
                    dataout.write('#%10s %11s %11s %11s \n' % ('F2(ppm)','F1(ppm)','Data','Init'))
                eval_data_recalc()
            #
            # Print out
            for slice in range(NoSlices):                
                if not ( zcoor=="3D" ):
                    for f1f2 in range(ndata/NoSlices):
                        #
                        # Only print out calculated spectrum, if fit is converged.
                        if ( lm_get('fit_converged') or ( 1==1) ):
                            dataout.write('%11.4e %11.4e %11.4e %11.4e \n' % ( \
                                eval_get_data(f1f2+slice*ndata/NoSlices)[0], \
                                eval_get_data(f1f2+slice*ndata/NoSlices)[1], \
                                SpecData[f1f2+slice*ndata/NoSlices], \
                                eval_get_data(f1f2+slice*ndata/NoSlices)[NDIM] \
                                ))
                        else:
                            dataout.write('%11.4e %11.4e %11.4e \n' % ( \
                                eval_get_data(f1f2+slice*ndata/NoSlices)[0], \
                                eval_get_data(f1f2+slice*ndata/NoSlices)[1], \
                                SpecData[f1f2+slice*ndata/NoSlices] \
                                ))
                        #
                        if ( f1f2+slice*ndata/NoSlices+1 < ndata ):
                            if ( math.fabs(\
                               eval_get_data(f1f2+slice*ndata/NoSlices)[1]- \
                               eval_get_data(f1f2+slice*ndata/NoSlices+1)[1] \
                               ) > 1e-10 ):
                                dataout.write('  \n')
                    dataout.write('  \n   \n')
            if not (zcoor=="3D" ):
                dataout.close()
            if ( zcoor=="3D"):
                NoSlices=0
                # First point is always written
                dataout.write('%11.4e %11.4e %11.4e %11.4e %11.4e \n' % ( \
                     eval_get_data(0)[0], \
                     eval_get_data(0)[1], \
                     eval_get_data(0)[2], \
                     SpecData[0], \
                     eval_get_data(0)[NDIM] \
                     ))
                for dc in range(1,ndata):
                    if not ( math.fabs(eval_get_data(dc)[1]-eval_get_data(dc-1)[1]) < 1e-10 ):
                        dataout.write("  \n")
                    if not ( math.fabs(eval_get_data(dc)[2]-eval_get_data(dc-1)[2]) < 1e-10 ):
                        dataout.write("  \n")
                        NoSlices+=1
                    dataout.write('%11.4e %11.4e %11.4e %11.4e %11.4e \n' % ( \
                         eval_get_data(dc)[0], \
                         eval_get_data(dc)[1], \
                         eval_get_data(dc)[2], \
                         SpecData[dc], \
                         eval_get_data(dc)[NDIM] \
                         ))              
                dataout.close()
                #sys.exit(10)
        #
        #If not converged, then set FIT_FLAG to -1.
        if ( not (lm_get('fit_converged'))) and ( not peak_fit_array[0]['shape']=="VOLUME") :
            if ( verboselevel > 0 ):
                sys.stdout.write('The peak(s): ')
            errorfile.write('The peak(s): ')
            for k in range (0,len(peak_fit_array)):
                if ( verboselevel >0 ):
                    sys.stdout.write('%s%s ' % (peakline[peak_fit_array[k]['id']][0],'  '))
                errorfile.write('%s%s ' % (peakline[peak_fit_array[k]['id']][0],'  '))
                peakline[peak_fit_array[k]['id']][4]=-1
            if ( verboselevel > 0 ):
                sys.stdout.write(' did not converge!!\n')
            errorfile.write(' did not converge!!\n')
            func_del_all()
            data_del_all()
            param_del_all()
            dtype_del_all()
        #
        if (lm_get('fit_converged') or peak_fit_array[0]['shape'] == "VOLUME" ):
            #Set the FIT_FLAG
            for k in range (0,len(peak_fit_array)):
                peakline[peak_fit_array[k]['id']][4]=1
            if not ( peak_fit_array[0]['shape']=="VOLUME" ):
                lm_update_param()
            #
            #Update spectral parameters and store in convenient array
            specparam=[]
            #
            NoParameters=[]
            for k in range(len(peak_fit_array)):
                NoParameters.append(0)
                if ( peak_fit_array[k]['shape'] == 'LORENTZIAN' or peak_fit_array[k]['shape']=='GAUSSIAN' ):
                    NoParameters[k]=4
                elif (peak_fit_array[k]['shape'] == "GLORE" ):
                    NoParameters[k]=6
                elif (peak_fit_array[k]['shape'] == "VOLUME" ):
                    NoParameters[k]=0
                else:
                    sys.stderr.write(" The shape %s is unknown \n" % peak_fit_array[k]['shape'] )
                    sys.stderr.write(" PROGRAM ABORTED \n")
                    sys.exit(1)
                if ( zcoor=="3D"):
                    NoParameters[k]=NoParameters[k]*3/2
                #
                #Check for isotopeshift
                if not ( peak_fit_array[k]['isotopeshift']=='n'):
                    if ( "fix" in peak_fit_array[k]['isotopeshift'] ):
                        NoParameters[k]+=1
                    else:
                        NoParameters[k]=NoParameters[k]+3
                        if ( zcoor =="3D"):
                            NoParameters[k]+=1
                    #
                    # If volume - no parameters!
                    if ( peak_fit_array[k]['shape']=="VOLUME" ):
                        NoParameters[k]=0
                    #
            NoSpecParam=0
            for k in range (0,len(peak_fit_array)):
                for parameter in range(NoParameters[k]):
                    specparam.append([eval_get_free(NoSpecParam),lm_get_value(NoSpecParam),lm_get_esd(NoSpecParam)])
                    NoSpecParam=NoSpecParam+1
            #
            #Store indentisies and errors in convenient array
            intensity=[]
            for k in range (0,len(peak_fit_array)):
                if ( zcoor == "3D" ):
                    intensity.append([\
                        lm_get_value(NoSpecParam+k),\
                        lm_get_esd(NoSpecParam+k)\
                        ])
                    continue
                #
                for i in range (0,NoSlices):
                    if not ( peak_fit_array[0]['shape'] == "VOLUME" ):
                        if ( baseline == 'y' or baseline == 'Y' ):
                            intensity.append([\
                               lm_get_value(NoSpecParam+k+(1+len(peak_fit_array))*i ), \
                               lm_get_esd(NoSpecParam+k+(1+len(peak_fit_array))*i )    \
                               ])    
                        else:
                            intensity.append([\
                                lm_get_value(NoSpecParam+k+(len(peak_fit_array))*i ), \
                                lm_get_esd(NoSpecParam+k+(len(peak_fit_array))*i )    \
                                ])
                    else:
                        intensity.append([Volumen[i],EsdVolumen[i]])
            #
            #
            for k in range (0,len(peak_fit_array)):
                #
                #Fit single exponential to data.    
                if (FITEXP == "Y" ):
                    data_del_all()
                    func_del_all()
                    dtype_del_all()
                    param_del_all()

                    # Clean up the data and func array
                    if not (param_exists('t')):
                        param('t',kind='EXPL')
                    if not (dtype_exists('dexp')):
                        dtype('dexp',('t',))
                    #
                    # Get appropriate start parameters as the slope:
                    R1_init=-(math.log(abs(intensity[taumin+NoSlices*k][0]))-math.log(abs(intensity[taumax+NoSlices*k][0])))/   \
                                (tau[taumin]-tau[taumax])
                    Amp_init=float(intensity[taumin+NoSlices*k][0])
                    #
                    # Define the parameters for the single exp. decay                
                    param('Amp',value=Amp_init)
                    param('R1',value=R1_init)
                    #
                    # Define the function
                    func('expdecay','exp_decay','dexp', ('Amp','R1','t'))
                    for i in range (0,NoSlices):
                        data(tau[i],\
                             float(intensity[i+NoSlices*k][0]),\
                             float(intensity[i+NoSlices*k][1])\
                             )
                    lm_minimize()
                    fitexp_amp=0
                    fitexp_r1=0
                    if (lm_get('fit_converged')):
                        lm_update_param()
                        #
                        # Write results to total logfile.
                        if ( lm_get('sd') < 1e-5 ):
                            errorfile.write(" WARNING! .. STD=%10.5e for single exp fit of %s \n" % ( lm_get('sd'),peakline[peak_fit_array[k]['id']][0]))
                            fitlog.write('%-15s %10.2f %10.2f %10.5f %10.5f %10.5f %s' % ( peakline[peak_fit_array[k]['id']][0],lm_get_value(0),lm_get_esd(0),lm_get_value(1),lm_get_esd(1),lm_get('sd'),'\n'))
                        else:
                            fitlog.write('%-15s %10.2f %10.2f %10.5f %10.5f %10.5f %s' % ( peakline[peak_fit_array[k]['id']][0],lm_get_value(0),lm_get_esd(0),lm_get_value(1),lm_get_esd(1),lm_get('sd'),'\n'))
                        #
                        # Store parameters
                        fitexp_amp=lm_get_value(0)
                        fitexp_r1=lm_get_value(1)
                        #
                        #Get the parameters to calculate the F-statistic .. if needed
                        single_nf=eval_get('ndata')-eval_get('nfree')
                        single_X2=math.pow(lm_get('enorm'),2)
                        #
                    else:
                        fitlog.write('#%s %s' % ( peakline[peak_fit_array[k]['id']][0],'  Single exponentiel fit did not converge .. !\n'))
                #
                #Fit single exponential 180 to data.    
                if (FIT180EXP == "Y" ):
                    data_del_all()
                    func_del_all()
                    dtype_del_all()
                    #
                    # Clean up the data and func array
                    if not (param_exists('t')):
                        param('t',kind='EXPL')
                    if not (dtype_exists('dexp')):
                        dtype('dexp',('t',))
                    #
                    # Get appropriate start parameters as the slope:
                    R1_init=(math.log(abs(intensity[taumin+NoSlices*k][0]))-2*math.log(abs(intensity[taumax+NoSlices*k][0])))/ \
                             (tau[taumin]-tau[taumax])
                    Amp1_init=float(intensity[taumax+FILESIZE[2]*k][0])
                    Amp0_init=float(intensity[taumin+FILESIZE[2]*k][0])-Amp1_init
                    #
                    # Define the parameters for the single exp. decay
                    param('Amp0',value=Amp0_init)
                    param('Amp1',value=Amp1_init)
                    param('R1',value=R1_init)
                    #
                    func('expdecay','exp_decay','dexp', ('Amp0','R1','t'))
                    func('cons','constant','dexp', ('Amp1',))
                    for i in range (0,NoSlices):
                        data(tau[i],float(intensity[i+NoSlices*k][0]),float(intensity[i+NoSlices*k][1]))
                    lm_minimize()
                    fitexp_amp0=0
                    fitexp_amp1=0
                    fitexp_r1=0
                    #
                    if (lm_get('fit_converged')):
                        lm_update_param()
                        # Write results to total logfile.
                        fit180log.write('%-10s %10.2f  %10.2f  %10.2f  %10.2f  %10.5f %10.5f %10.5f %s' % ( peakline[peak_fit_array[k]['id']][0],lm_get_value(0),lm_get_esd(0),lm_get_value(1),lm_get_esd(1),lm_get_value(2),lm_get_esd(2),lm_get('sd'),'\n'))
                        fitexp_amp0=lm_get_value(0)
                        fitexp_amp1=lm_get_value(1)
                        fitexp_r1=lm_get_value(2)
                        #
                        #Get the parameters to calculate the F-statistic
                        single_nf=eval_get('ndata')-eval_get('nfree')
                        single_X2=math.pow(lm_get('enorm'),2)
                    else:
                        fit180log.write('%s %s' % ( peakline[peak_fit_array[k]['id']][0],'  Single exponentiel fit did not converge .. !\n'))
                #
                #Fit the bi-exponential, if flag is set:
                if (FITBIEXP =="Y"):
                    #Start bi fixing the known parameters C1,L1
                    data_del_all()
                    func_del_all()
                    #param_del_all()
                    if not (param_exists('t')):
                        param('t',kind='EXPL')
                    if not (dtype_exists('dexp')):
                        dtype('dexp',('t',))
                    #
                    # Get appropriate start parameters as the slope:
                    R1_init=-(math.log(abs(intensity[tausort[0]+NoSlices*k][0]))-math.log(abs(intensity[tausort[1]+NoSlices*k][0])))/ \
                             (tau[tausort[0]]-tau[tausort[1]])
                    Amp_init=float(intensity[taumin+NoSlices*k][0])
                    #
                    param('C1_init',value=Amp_init)
                    param('L1_init',value=R1_init)
                    #
                    # First function in the bi-exp decay
                    func('expdecay1','exp_decay','dexp', ('C1_init','L1_init','t'))
                    for i in range (0,4):
#                    for i in range (len(tau)-6,len(tau)-1):
                        data(tau[tausort[i]],intensity[tausort[i]+NoSlices*k][0],intensity[tausort[i]+NoSlices*k][1])
                    #
                    lm(tol=1.0e-12,maxfev=100)
                    lm_minimize()
                    #
                    if (lm_get('fit_converged')):
                        lm_update_param()
                        fitshortlog.write('%-10s %10.2f  %10.2f  %10.5f  %10.5f %10.5f %s' % ( peakline[peak_fit_array[k]['id']][0],lm_get_value(0),lm_get_esd(0),lm_get_value(1),lm_get_esd(1),lm_get('sd'),'\n'))
                        #Store the values to put in gnuplot file
                        shortamp=[]
                        shortexp=[]
                        shortamp.append(lm_get_value(0))
                        shortamp.append(lm_get_esd(0))
                        shortexp.append(lm_get_value(1))
                        shortexp.append(lm_get_esd(1))
                    #
                    #Clean up to next fit.
                    data_del_all()
                    func_del_all()
                    #
                    if (param_exists('C1_2')):
                        param_del('C1_2')
                    if (param_exists('L1_2')):
                        param_del('L1_2')
                    #
                    param('C1_2',value=lm_get_value(0),kind='CONST')
                    param('L1_2',value=lm_get_value(1),kind='CONST')
                    param('C2_init',value=-lm_get_value(0)/10)
                    param('L2_init',value=lm_get_value(1)*5)
                    #
                    func_del_all()
                    data_del_all()
                    #
                    func('expdecay1','exp_decay','dexp', ('C1_2','L1_2','t'))
                    func('expdecay2','exp_decay','dexp', ('C2_init','L2_init','t'))
                    #
                    for i in range (0,NoSlices):
                        data(tau[i],float(intensity[i+NoSlices*k][0]),float(intensity[i+NoSlices*k][1]))
                    #
                    lm_minimize()
                    #
                    if (lm_get('fit_converged')):
                        lm_update_param()
                    #
                    func_del_all()
                    data_del_all()
                    #
                    param('C1',value=param_get('C1_2','value'))
                    param('L1',value=param_get('L1_2','value'))
                    param('C2',value=lm_get_value(0))
                    param('L2',value=lm_get_value(1))
                    #
                    func('expdecay1','exp_decay','dexp', ('C1','L1','t'))
                    func('expdecay2','exp_decay','dexp', ('C2','L2','t'))
                    #
                    for i in range (0,NoSlices):
                        data(tau[i],float(intensity[i+NoSlices*k][0]),float(intensity[i+NoSlices*k][1]))
                    #
                    lm(tol=1.0e-10,maxfev=50)
                    lm_minimize()
                    #
                    fitbiexp_c1=0
                    fitbiexp_l1=0
                    fitbiexp_c2=0
                    fitbiexp_l2=0
                    #
                    if (lm_get('fit_converged')):
                        lm_update_param()
                        #
                        #Get data to calculate F-Statistic. 
                        bi_nf=eval_get('ndata')-eval_get('nfree')
                        bi_X2=math.pow(lm_get('enorm'),2)
                        #
                        #Calculate F:
                        F=bi_nf*(single_X2-bi_X2)/((single_nf-bi_nf)*bi_X2)
                        #
                        #
                        if (F>0.02):
                            fitbilog.write('%-10s ' % ( peakline[peak_fit_array[k]['id']][0]))
                            fitbilog.write(' %10.2f %10.2f '  % ( lm_get_value(0),lm_get_esd(0)))
                            fitbilog.write(' %10.5f %10.5f '  % ( lm_get_value(1),lm_get_esd(1)))
                            fitbilog.write(' %10.2f %10.2f '  % ( lm_get_value(2),lm_get_esd(2)))
                            fitbilog.write(' %10.5f %10.5f '  % ( lm_get_value(3),lm_get_esd(3)))
                            fitbilog.write(' %10.5f %s ' % (lm_get('sd'),'  '))
                            fitbilog.write(' %10.5f %s' % (F,'\n'))
                        else:
                            print '%s  Bi exp. fit has F= %10.2f \n' % ( peakline[peak_fit_array[k]['id']][0],F)
                        #
                        fitbiexp_c1=lm_get_value(0)
                        fitbiexp_l1=lm_get_value(1)
                        fitbiexp_c2=lm_get_value(2)
                        fitbiexp_l2=lm_get_value(3)
                    #
                    if not (lm_get('fit_converged')):
                        print '%s  Bi exp. fit did not converge\n' % ( peakline[peak_fit_array[k]['id']][0])
                #
                # Chek if we are to dump parameters
                if ( dumpparameters ):
                   outputfile=sys.argv[2]+"/"+peakline[peak_fit_array[k]['id']][0]+".par" 
                   parfile = open(outputfile,'w')
                   SpecStart=0
                   for i in range(0,k):
                       SpecStart=SpecStart+NoParameters[i]
                   for i in range(NoParameters[k]):
                        parfile.write("%s %20.13e\n" % (\
                            string.split(specparam[i+SpecStart][0],'_')[0],\
                            specparam[i+SpecStart][1],\
                            ))
                   #
                   if ( not ( zcoor=="2D" or zcoor=="3D" )):
                       ic=0
                       for i in range(1,NoSlices+1+len(discard_slices)):
                           #
                           # Check if this slice should be excluded:
                           discard_slice=0
                           for slice in discard_slices:
                               if ( int(slice) == i ):
                                   discard_slice=1
                           if ( discard_slice ):
                               continue
                           parfile.write("I%d %20.13e \n" % (i,intensity[ic+NoSlices*k][0])) 
                           ic=ic+1
                           #
                   else:
                       if ( zcoor=="3D" ):
                           parfile.write("I1 %20.13e \n" % ( intensity[k][0]))
                       else:
                           parfile.write("I1 %20.13e \n" % (\
                                  intensity[NoSlices*k][0],\
                                  ))

                #
                # Check if we are writing individual outputfiles for the peaks
                if ( write_individual ):
                    outputfile=sys.argv[2]+"/"+peakline[peak_fit_array[k]['id']][0]+".out"
                    outfile = open(outputfile,'w')
                    outfile.write('#\n')
                    outfile.write("# %15s%15s\n" %( "Peak Name","Overlap_group"))
                    outfile.write("# %15s%15d\n" % ( \
                        peakline[peak_fit_array[k]['id']][0],
                        int(peakline[peak_fit_array[k]['id']][5]),
                        ))
                    outfile.write('#\n')
                    outfile.write("# Input frequencies:\n")
                    outfile.write("# %10s%10s%10s\n" % ("Omega1","Omega2","Omega3"))
                    outfile.write("# %10.3f%10.3f%10.3f\n" % (\
                        float(peakline[peak_fit_array[k]['id']][1]),
                        float(peakline[peak_fit_array[k]['id']][2]),
                        float(peakline[peak_fit_array[k]['id']][3]),
                        ))
                    outfile.write('#\n# --------- Results of the fit -------------\n#\n')
                    outfile.write('# Parameter        Value           Esd      \n')   
                    SpecStart=0
                    for i in range(0,k):
                        SpecStart=SpecStart+NoParameters[i]
                    #print 'We have #Parameters=',NoParameters[k]
                    #for P in range(eval_get('nfree')):
                    #    print eval_get_free(P),param_get(eval_get_free(P),'value')
                    #sys.exit(10)
                    for i in range(NoParameters[k]):
                        prename=specparam[i+SpecStart][0][0]
                        if ( prename == "j" or prename == "w" ):
                            prename=string.split(specparam[i+SpecStart][0],'_')[0]+"(Hz)"
                            dimension=int(specparam[i+SpecStart][0][1])
                            if ( dimension==1 ):
                                val=specparam[i+SpecStart][1]*sfrq[0]
                                esd=specparam[i+SpecStart][2]*sfrq[0]
                            elif ( dimension == 2 ):                                
                                val=specparam[i+SpecStart][1]*sfrq[1]
                                esd=specparam[i+SpecStart][2]*sfrq[1]
                            elif ( dimension == 3):
                                val=specparam[i+SpecStart][1]*sfrq[2]
                                esd=specparam[i+SpecStart][2]*sfrq[2]
                            else:
                                sys.stderr.write("\n The parameter %s cannot be decomposed into dimension and linewidth/coupling\n PROGRAM ABORTED\n" % ( specparam[i+SpecStart][0] ))
                                sys.exit(1)
                        elif ( prename == "f" ):
                            prename=string.split(specparam[i+SpecStart][0],'_')[0]+"(ppm)"
                            val=specparam[i+SpecStart][1]
                            esd=specparam[i+SpecStart][2]
                        else:
                            prename=string.split(specparam[i+SpecStart][0],'_')[0]
                            val=specparam[i+SpecStart][1]
                            esd=specparam[i+SpecStart][2]
                            
                        if ( not ( zcoor=="2D" or zcoor == "3D" )):
                            outfile.write( '%s %-12s %14.7e %14.7e %s' % ('#',\
                                                                          prename,
                                                                          val,\
                                                                          esd,'\n') )
                        else:
                            outfile.write( '%s %-12s %14.7e %14.7e %s' % (' ',\
                                                                          prename,
                                                                          val,\
                                                                          esd,'\n') )                            
                        #
                    if ( not ( zcoor=="2D" or zcoor=="3D" )):
                        for i in range (0,44):
                            outfile.write('#')
                        outfile.write('\n')
                        #
                        #outfile.write('# Z-coordinate      Intensity      Esd(Int.)\n')
                        outfile.write('# %10s        Intensity      Esd(Int.)\n' % (zcoor))
                        for i in range(0,NoSlices):
                            outfile.write('%12.3e   %14.7e %14.7e' % (float(tau[i]),intensity[i+NoSlices*k][0],intensity[i+NoSlices*k][1])+'\n') 
                        outfile.close()
                    else:    
                        if ( zcoor=="3D" ):
                            outfile.write('%s %-12s %14.7e %14.7e %s' % (' ', \
                                 'Intensity',
                                  intensity[k][0],\
                                  intensity[k][1],\
                                 '\n'\
                                 ) )
                        else:
                            outfile.write('%s %-12s %14.7e %14.7e %s' % (' ', \
                                 'Intensity',
                                  intensity[NoSlices*k][0],\
                                  intensity[NoSlices*k][1],\
                                 '\n'\
                                 ) )
                    outfile.close()
                if ( LOGFILENR=="Y" ):
                    SpecStart=0
                    for i in range(0,k):
                        SpecStart=SpecStart+NoParameters[i]
                    #
                    # Parameter order {f1,f2,w1,w2,intensity,isotope(y/n)}
                    ParamOrder=[]
                    for i in range(NoParameters[k]):
                        ParamOrder.append(-1)
                    # 
                    for i in range(NoParameters[k]):
                        if ( zcoor=="3D" ):
                            if ( specparam[SpecStart+i][0][0:3] == 'f01' ):
                                ParamOrder[0]=i+SpecStart
                            if ( specparam[SpecStart+i][0][0:3] == 'f02' ):
                                ParamOrder[1]=i+SpecStart
                            if ( specparam[SpecStart+i][0][0:3] == 'f03' ):
                                ParamOrder[2]=i+SpecStart
                            if ( specparam[SpecStart+i][0][0:2] == 'w1' ):
                                ParamOrder[3]=i+SpecStart
                            if ( specparam[SpecStart+i][0][0:2] == 'w2' ):
                                ParamOrder[4]=i+SpecStart                            
                            if ( specparam[SpecStart+i][0][0:2] == 'w3' ):
                                ParamOrder[5]=i+SpecStart                            
                        elif ( zcoor=="2D" ):
                            if ( specparam[SpecStart+i][0][0:3] == 'f01' ):
                                ParamOrder[0]=i+SpecStart
                            if ( specparam[SpecStart+i][0][0:3] == 'f02' ):
                                ParamOrder[1]=i+SpecStart
                            if ( specparam[SpecStart+i][0][0:2] == 'w1' ):
                                ParamOrder[2]=i+SpecStart
                            if ( specparam[SpecStart+i][0][0:2] == 'w2' ):
                                ParamOrder[3]=i+SpecStart
                        else:
                            sys.stderr.write("\n Cannot write non-array logfile\n")
                            sys.stderr.write("                      PROGARM ABORTED\n")
                    #
                    fitlog.write("%14s " % ( peakline[peak_fit_array[k]['id']][0]))
                    if ( zcoor=="3D" ):
                        fitlog.write("%10.5f %10.5f " % ( specparam[ParamOrder[0]][1],specparam[ParamOrder[0]][2]))
                        fitlog.write("%10.5f %10.5f " % ( specparam[ParamOrder[1]][1],specparam[ParamOrder[1]][2]))
                        fitlog.write("%10.5f %10.5f " % ( specparam[ParamOrder[2]][1],specparam[ParamOrder[2]][2]))
                        fitlog.write("%10.3e %10.3e " % ( specparam[ParamOrder[3]][1]*sfrq[0],specparam[ParamOrder[3]][2]*sfrq[0]))
                        fitlog.write("%10.3e %10.3e " % ( specparam[ParamOrder[4]][1]*sfrq[1],specparam[ParamOrder[4]][2]*sfrq[1]))
                        fitlog.write("%10.3e %10.3e " % ( specparam[ParamOrder[5]][1]*sfrq[2],specparam[ParamOrder[5]][2]*sfrq[2]))
                        fitlog.write("%10.3e %10.3e " % ( intensity[k][0],intensity[k][1]))
                    else:
                        fitlog.write("%10.5f %10.5f " % ( specparam[ParamOrder[0]][1],specparam[ParamOrder[0]][2]))
                        fitlog.write("%10.5f %10.5f " % ( specparam[ParamOrder[1]][1],specparam[ParamOrder[1]][2]))
                        fitlog.write("%10.3e %10.3e " % ( specparam[ParamOrder[2]][1]*sfrq[0],specparam[ParamOrder[2]][2]*sfrq[0]))
                        fitlog.write("%10.3e %10.3e " % ( specparam[ParamOrder[3]][1]*sfrq[1],specparam[ParamOrder[3]][2]*sfrq[1]))
                        fitlog.write("%10.3e %10.3e " % ( intensity[NoSlices*k][0],intensity[NoSlices*k][1]))
                    #
                    if ( peak_fit_array[k]['isotopeshift'] == 'n' ):
                        fitlog.write("%5s" % ( 'n' ))
                    else:
                        fitlog.write("%5s" % ( 'y' ))
                    #    
                    fitlog.write("\n")
                #
                #
                # Write GnuplotScript file; single exponentiel
                if (GNUPLOT == "Y" and ( FITEXP == "Y" or FIT180EXP == "Y") ):
                    gnuoutfile=sys.argv[2]+"/"+peakline[peak_fit_array[k]['id']][0]+".gnu"
                    gnuout = open(gnuoutfile,'w')
                    gnuout.write('%s %s '% ("# Gnuplot script to fit relaxation rate:","\n"))
                    if ( FIT180EXP == "Y" ):
                        gnuout.write('%s %s '% ("f(x)=m2-(m2-m1)*exp(-r1*x)","\n"))
                        gnuout.write('%s %9.4f %s '% ( "m1=",intensity[taumin+NoSlices*k][0],"\n" ))
                        gnuout.write('%s %9.4f %s '% ( "m2=",intensity[taumax+NoSlices*k][0],"\n" ))
                    else:
                        gnuout.write('%s %s '% ("f(x)=m1*exp(-r1*x)","\n"))
                        gnuout.write('%s %9.4f %s '% ( "m1=",intensity[taumin+NoSlices*k][0],"\n" ))
                    #
                    gnuout.write('%s %s '% ("set zero 1e-40","\n"))
                    gnuout.write('%s %s '% ("FIT_LIMIT=1e-20","\n"))
                    gnuout.write('%s %s '% ("r1=1.5","\n"))
                    #gnuout.write('%s %s ' %("set logscale x","\n"))
                    ## Single exponential fit
                    #
                    if ( (FITEXP =="Y" ) and (fitexp_amp+fitexp_r1 != 0)):
                        gnuout.write('%s%s%s%s '% ("fit f(x) '",peakline[peak_fit_array[k]['id']][0],".out' u 1:2:3 via m1,r1","\n")) 
                        gnuout.write('%s%s%s%s'% ("plot '",peakline[peak_fit_array[k]['id']][0],".out' u 1:2:3 w e, f(x) t ''","\n")) 
                        gnuout.write('%s%s'% ("pause -1 'Press [ENTER] key to exit...'","\n"))
                    #
                    # Single exponential fit
                    #####
                    if ( (FIT180EXP =="Y" ) and (fitexp_amp0+fitexp_amp1+fitexp_r1 != 0)):
                        gnuout.write('%s%s%s%s '% ("fit f(x) '",peakline[peak_fit_array[k]['id']][0],".out' u 1:2:3 via m2,m1,r1","\n")) 
                        gnuout.write('%s%s%s%s'% ("plot '",peakline[peak_fit_array[k]['id']][0],".out' u 1:2:3 w e, f(x) t ''","\n")) 
                        gnuout.write('%s%s'% ("pause -1 'Press [ENTER] key to exit...'","\n"))
                    #
                    if (FITBIEXP == "Y" and (fitbiexp_c1+fitbiexp_c2+fitbiexp_l1+fitbiexp_l2 != 0)):
                        gnuout.write('%s%s%s%10.5f%s%10.5f%s%10.5f%s%10.5f%s %s' % ("plot '",peakline[peak_fit_array[k]['id']][0],".out' u 1:2:3 w e,",fitbiexp_c1,"*exp(-",fitbiexp_l1,"*x)+ ",fitbiexp_c2,"*exp(-",fitbiexp_l2,"*x) t 'Bi-exp fit'","\n"))
                        gnuout.write('%s%s'% ("pause -1 'Press [ENTER] key to exit...'","\n"))
                        gnuout.write('%s%s%s%10.5f%s%10.5f%s %s' % ("plot '",peakline[peak_fit_array[k]['id']][0],".out' u 1:2:3 w e,",shortamp[0],"*exp(-",shortexp[0],"*x) t 'Short-exp fit'","\n"))
                        gnuout.write('%s%s'% ("pause -1 'Press [ENTER] key to exit...'","\n"))
                    if ( verboselevel > 4 ):
                        sys.stdout.write(" Gnuplot ScriptFile has been written to the file %s\n" % (gnuoutfile))
                    gnuout.close()
                #
                # Write GnuplotScript file; single exponentiel
                if (GNUPLOT == "Y" and FITLIN == "Y" ):
                    gnuoutfile=sys.argv[2]+"/"+peakline[peak_fit_array[k]['id']][0]+".gnu"
                    gnuout = open(gnuoutfile,'w')
                    gnuout.write('%s %s '% ("# Gnuplot script to fit linear:","\n"))
                    gnuout.write('%s %s '% ("f(x)=a+b*x","\n"))
                    gnuout.write('%s %s '% ("set zero 1e-40","\n"))
                    gnuout.write('%s %s '% ("FIT_LIMIT=1e-20","\n"))
                    gnuout.write('%s%s%s%s '% ("fit f(x) '",peakline[peak_fit_array[k]['id']][0],".out' u 1:2:3 via a,b","\n")) 
                    gnuout.write('%s%s%s%s'% ("plot '",peakline[peak_fit_array[k]['id']][0],".out' u 1:2:3 w e, f(x) t ''","\n"))
                    gnuout.write('%s%s'% ("pause -1 'Press [ENTER] key to exit...'","\n"))
                    gnuout.close()
            #
            if (( printdata=="Y" or printdata == "y" ) and GNUPLOT == "Y" ):
                gnuoutfile=sys.argv[2]+"/"+peakline[peak_fit_array[0]['id']][0]+".dat.gnu"
                gnuout=open(gnuoutfile,'w')
                gnuout.write('%s %s '% ("# Gnuplot script to plot data:","\n"))
                #gnuout.write('unset ztics\n')
                datafile=peakline[peak_fit_array[0]['id']][0]+".dat"
                groupname=peakline[peak_fit_array[0]['id']][0]
                for slices in range(NoSlices):
                    if ( 0 == 0 ):
                        if ( not ( zcoor=="2D" or zcoor == "3D" )):
                            gnuout.write('splot \'%s\' i %d u 1:2:3 t \'Experimental\' w linesp, \'%s\' i %d u 1:2:4 w lines t \'Group: %s; Tau %13.6es\' \n' %(datafile,slices,datafile,slices,groupname,tau[slices]))
                        else:
                            if ( zcoor=="2D" ):
                                gnuout.write('splot \'%s\' i %d u 1:2:3 t \'Experimental\' w linesp, \'%s\' i %d u 1:2:4 w lines t \'Group: %s \' \n' %(datafile,slices,datafile,slices,groupname))
                            else:
                                gnuout.write('splot \'%s\' i %d u 1:2:4 t \'Experimental\' w linesp, \'%s\' i %d u 1:2:5 w lines t \'Group: %s; plane %13d(ppm)\' \n' %(datafile,slices,datafile,slices,groupname,slices))                        
                        gnuout.write('pause -1 \'Press [ENTER] key to continue\'\n')
                gnuout.close()

#
# Check gnuplot version
os.system('gnuplot --version >& templine ')
gpversion=string.split(open('templine','r').readline())
if ( gpversion[0] == "gnuplot" ):
    if ( float(gpversion[1][0:3]) > 3.7 ):
        if gpversion[1][3:4] < 'i' and float(gpversion[1][0:3]) < 3.9 :
            sys.stderr.write("\n\n *******************************************************************\n")
            sys.stderr.write(" You have installed GnuPlot version %s. \n" % (gpversion[1]))
            sys.stderr.write(" To properly run the scripts generated by this program,\n")
            sys.stderr.write(" you will need GnuPlot version > 3.8j, which can be downloaded from \n")
            sys.stderr.write(" http://www.gnuplot.info/ \n")
            sys.stderr.write(" *******************************************************************\n")
else:
    os.system('/usr/bin/gnuplot --version >& templine ')
    gpversion=string.split(open('templine','r').readline())
    if ( gpversion[0] == "gnuplot" ):
        if ( float(gpversion[1][0:3]) > 3.7 ):
            if gpversion[1][3:4] > 'b' :
                os.system('which gnuplot >& templine ')
                gppointer=string.split(open('templine','r').readline())
                sys.stderr.write("\n\n *******************************************************************\n")
                sys.stderr.write(" It look like the reference to your gnuplot executeable is wrong!\n")
                sys.stderr.write(" The reference \'gnuplot' points to %s,\n" % (gppointer[0]))
                sys.stderr.write(" which is an extremely old version of the GnuPlot program.\n")
                sys.stderr.write(" However, you have a never version of GnuPlot located in /usr/bin/gnuplot\n")
                sys.stderr.write(" Please update your alias/pointers to fix this problem,\n")
                sys.stderr.write(" in order to probably run the scripts generated by this program\n")                
                sys.stderr.write(" *******************************************************************\n")
    else:
        if ( verboselevel > 2 ):
            sys.stderr.write("\n\n *******************************************************************\n")
            sys.stderr.write(" It looks like you do not have GnuPlot installed, or an extremely\n")
            sys.stderr.write(" old version is installed. This means that you will probably not be\n")
            sys.stderr.write(" able to run the scripts generated by this program\n" )
            sys.stderr.write(" You will need GnuPlot version > 3.8j, which can be downloaded from \n")
            sys.stderr.write(" http://www.gnuplot.info/ \n")
            sys.stderr.write(" *******************************************************************\n")
os.system('rm -fR templine')
