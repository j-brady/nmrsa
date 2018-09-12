import os
from math import log10
import nmrglue as ng
import argparse as ap
from jinja2 import Environment, PackageLoader


if __name__ == "__main__":

    parser = ap.ArgumentParser(description="Script to read bruker experiment directory and extract paramers for human readability")
    parser.add_argument("--dirname","-d",
            type=str,
            help="Name of directory containing experimental data")
    parser.add_argument("--outname","-o",
            type=str,
            default="params.txt",
            help="Name of output text file. Default = params.txt")

    parser.add_argument("--notes","-n",
            type=str,
            help="Write some notes to remind yourself what you were doing")

    args = parser.parse_args()

    env = Environment(loader=PackageLoader('bruk2par', 'templates'))
    template = env.get_template('params2col.txt')

    dirname = args.dirname
    notes = args.notes
    out_name = args.outname

    dic, data = ng.bruker.read(dirname)
    #print(dic["acqus"])
    pp_dic = ng.bruker.read_pprog(os.path.join(dirname,"pulseprogram"))
    pp_vars = pp_dic["var"]
    print(pp_vars)
    shape_pulses = [int(i.strip("spw")) for i in pp_vars.values() if i.startswith("spw")]
    pws = [i for i in pp_dic["var"].keys() if i.startswith("pw")]
    pw_numbers = [int(pp_vars[i].strip()[1:]) for i in pws]
    pws = [[i,
            num,dic["acqus"]["P"][num],
            dic["acqus"]["PLW"][num],
            -10.*log10(float(dic["acqus"]["PLW"][num]))] for num,i in zip(pw_numbers,pws) if dic["acqus"]["PLW"][num]!=0]

    delays = [i for i in pp_dic["var"].keys() if i.startswith("d")]
    plws = [(num,["%.2f"%(-10*log10(float(i))),i]) for num,i in enumerate(dic["acqus"]["PLW"]) if i !=0 ]
    grads = [[num,x[0],x[1],x[2]] for num,x in enumerate(zip(dic["acqus"]["GPX"],dic["acqus"]["GPY"],dic["acqus"]["GPZ"])) \
            if (sum([x[0],x[1],x[2]])) != 0]
    cpdprg = [[num,i] for num,i in enumerate(dic["acqus"]["CPDPRG"]) if i!=""]
    render_dic = {}
    for i in dic.keys():
        if i.startswith("acqu"):
            render_dic[i] = dic[i]
        #print(dic.keys())
    #print(dic["acqu2s"])#.keys())

    out = template.render(dic=render_dic,
                          pp_vars=pp_vars,
                          cpdprg=cpdprg,
                          pws=pws,
                          delays=delays,
                          plws=plws,
                          shape_pulses=shape_pulses,
                          grads=grads)
    print(out)
    try:
        outfile = open(os.path.join(dirname,out_name),"w")
        print("Parameters written into %s"%out_name)
        outfile.write(out)
    except PermissionError:
        sp.call("gedit",shell=True)
    outfile.close()
