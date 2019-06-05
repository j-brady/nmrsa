#!/Users/jacobbrady/virtual_envs/py37/bin/python
"""
	Usage: spec.py <yaml_file>
           spec.py make <new_yaml_file>

	Plot NMRPipe spectra overlays using nmrglue and matplotlib. This is my attempt to make a general script for
    plotting NMR data.
    
    
    Below is an example yaml file for input
 
        # This first block is global parameters which can be overridden by adding the desired argument
        # to your list of spectra. One exception is "colors" which if set in global params overrides the 
        # color option set for individual spectra as the colors will now cycle through the chosen matplotlib 
        # colormap
		cs: 10e5                        # contour start
		contour_num: 10                 # number of contours
		contour_factor: 1.2             # contour factor
		colors: Set1                    # must be matplotlib.cm colormap

		outname: ["overlay_wt_s95a.pdf","overlay_wt_s95a.png"] # either single value or list of output names 

        # Here is where your list of spectra to plot goes
		spectra: 

				- fname: test.ft2
				  label: shHTL5(S95A) pH 6.5 - 30 degrees 
				  contour_num: 1
				  linewidths: 1

	Options:
	    -h --help
        -v --version


    Dependencies:

        -- python3
        -- matplotlib, pyyaml, numpy, nmrglue and docopt
        -- pandas (if plotting clusters)


    TODO:

        -- enable N dim spectra to be plotted
        -- make smart contour start setting


"""
import os

import yaml
import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm
from docopt import docopt

yaml_file = """
##########################################################################################################
#  This first block is global parameters which can be overridden by adding the desired argument          #
#  to your list of spectra. One exception is "colors" which if set in global params overrides the        #
#  color option set for individual spectra as the colors will now cycle through the chosen matplotlib    #
#  colormap                                                                                              #         
##########################################################################################################

cs: 10e5                        # contour start
contour_num: 10                 # number of contours
contour_factor: 1.2             # contour factor
colors: Set1                    # must be matplotlib.cm colormap

outname: ["overlay_wt_s95a.pdf","overlay_wt_s95a.png"] # either single value or list of output names 
ncol: 1 #  tells matplotlib how many columns to give the figure legend - if not set defaults to 2

# Here is where your list of spectra to plot goes
spectra: 

        - fname: test.ft2
          label: shHTL5(S95A) pH 6.5 - 30 degrees 
          contour_num: 1
          linewidths: 1
"""


def make_yaml_file(name, yaml_file=yaml_file):

    if os.path.exists(name):
        print(f"Copying {name} to {name}.bak")
        os.system(f"cp {name} {name}.bak")

    print(f"Making yaml file ... {name}")
    with open(name, "w") as new_yaml_file:
        new_yaml_file.write(yaml_file)


if __name__ == "__main__":
    arguments = docopt(__doc__, version="Spec 0.1")
    if arguments["make"]:
        make_yaml_file(name=arguments["<new_yaml_file>"])
        exit()

    params = yaml.load(open(arguments["<yaml_file>"], "r"))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    cs_g = float(params["cs"])
    spectra = params["spectra"]
    contour_num_g = params.get("contour_num", 10)
    contour_factor_g = params.get("contour_factor", 1.2)
    nspec = len(spectra)
    notes = []

    for num, spec in enumerate(spectra):

        # unpack spec specific parameters
        fname = spec["fname"]

        if params.get("colors"):
            # currently overrides color option
            color = np.linspace(0, 1, nspec)[num]
            colors = cm.get_cmap(params.get("colors"))(color)
            # print("Colors set to cycle though %s from Matplotlib"%params.get("colors"))
            # print(colors)
            colors = colors[:-1]

        else:
            colors = spec["colors"]

        neg_colors = spec.get("neg_colors")
        label = spec["label"]
        cs = float(spec.get("cs", cs_g))
        contour_num = spec.get("contour_num", contour_num_g)
        contour_factor = spec.get("contour_factor", contour_factor_g)
        #  append cs and colors to notes
        notes.append((cs, colors))

        # read spectra
        dic, data = ng.pipe.read(fname)
        udic = ng.pipe.guess_udic(dic, data)

        ndim = udic["ndim"]

        if ndim == 1:
            uc_f1 = ng.pipe.make_uc(dic, data, dim=0)

        elif ndim == 2:
            f1, f2 = params.get("dims", [0, 1])
            uc_f1 = ng.pipe.make_uc(dic, data, dim=f1)
            uc_f2 = ng.pipe.make_uc(dic, data, dim=f2)

            ppm_f1 = uc_f1.ppm_scale()
            ppm_f2 = uc_f2.ppm_scale()

            ppm_f1_0, ppm_f1_1 = uc_f1.ppm_limits()  # max,min
            ppm_f2_0, ppm_f2_1 = uc_f2.ppm_limits()  # max,min

        elif ndim == 3:
            f1, f2, f3 = params.get("dims", [0, 1, 2])
            uc_f1 = ng.pipe.make_uc(dic, data, dim=f1)
            uc_f2 = ng.pipe.make_uc(dic, data, dim=f2)
            uc_f3 = ng.pipe.make_uc(dic, data, dim=f3)
            # need to make more robust
            ppm_f1 = uc_f2.ppm_scale()
            ppm_f2 = uc_f3.ppm_scale()

            ppm_f1_0, ppm_f1_1 = uc_f2.ppm_limits()  # max,min
            ppm_f2_0, ppm_f2_1 = uc_f3.ppm_limits()  # max,min

            if f1 == 0:
                data = data[f1]

            # x and y are set to f2 and f1
            f1, f2 = f2, f3
            #elif f1 == 1:
            #    data = data[:,0,:]
            #else:
            #    data = data[:,:,0]

        # plot parameters
        contour_start = cs  # contour level start value
        contour_num = contour_num  # number of contour levels
        contour_factor = contour_factor  # scaling factor between contour levels

        # calculate contour levels
        cl = contour_start * contour_factor ** np.arange(contour_num)
        ax.contour(
            data,
            cl,
            colors=[colors for _ in cl],
            linewidths=spec.get("linewidths", 0.5),
            extent=(ppm_f2_0, ppm_f2_1, ppm_f1_0, ppm_f1_1),
        )

        if neg_colors:
            ax.contour(
                data * -1,
                cl,
                colors=[neg_colors for _ in cl],
                linewidths=spec.get("linewidths", 0.5),
                extent=(ppm_f2_0, ppm_f2_1, ppm_f1_0, ppm_f1_1),
            )

        else:  # if no neg color given then plot with 0.5 alpha
            ax.contour(
                data * -1,
                cl,
                colors=[colors for _ in cl],
                linewidths=spec.get("linewidths", 0.5),
                extent=(ppm_f2_0, ppm_f2_1, ppm_f1_0, ppm_f1_1),
                alpha=0.5,
            )

        # hack for legend
        ax.plot([], [], c=colors, label=label)

    # plt.xlim(ppm_f2_0, ppm_f2_1)
    ax.invert_xaxis()
    ax.set_xlabel(udic[f2]["label"] + " ppm")
    if params.get("xlim"):
        ax.set_xlim(*params.get("xlim"))

    # plt.ylim(ppm_f1_0, ppm_f1_1)
    ax.invert_yaxis()
    ax.set_ylabel(udic[f1]["label"] + " ppm")

    plt.legend(
        loc="upper center", bbox_to_anchor=(0.5, 1.20), ncol=params.get("ncol", 2)
    )
    # plt.tight_layout()
    #  add a list of outfiles
    y = 0.025
    # only write cs levels if show_cs: True in yaml file
    if params.get("show_cs"):
        for num, j in enumerate(notes):
            col = j[1]
            con_strt = j[0]
            ax.text(0.025, y, "cs=%.2e" % con_strt, color=col, transform=ax.transAxes)
            y += 0.05

    if params.get("clusters"):
        import pandas as pd
        clusters = pd.read_pickle(params.get("clusters"))
        groups = clusters.groupby("CLUSTID")
        for ind, group in groups:
            if len(group)==1:
                plt.plot(group.X_PPM, group.Y_PPM, "ko", markersize=2)
            else:
                plt.plot(group.X_PPM, group.Y_PPM, "o",markersize=2)

    if params.get("outname") and (type(params.get("outname")) == list):
        for i in params.get("outname"):
            plt.savefig(i, bbox_inches="tight",dpi=300)
    else:
        plt.savefig(params.get("outname", "test.pdf"), bbox_inches="tight")
    plt.show()
