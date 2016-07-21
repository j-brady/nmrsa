#!/usr/bin/python

import nmrglue as ng

dic,data = ng.fileio.bruker.read()

with open("parameters.txt",'w') as o:
    comments = dic["acqus"]["_comments"] 
    #max_len = max([len(i) for i in comments])
    for i in comments:
        o.write("%s\n"%i)

    o.write("\n\n")
    o.write("%s\n\n"%dic["acqus"]["PULPROG"])

    o.write("%8s    %8s    %8s\n"%("Param","bruk_name","value"))
    o.write("%8s    %8s    %8f us\n"%("Proton90","p1",dic['acqus']['P'][1]))
    o.write("%8s    %8s    %8f us\n"%("Carbon90","p2",dic['acqus']['P'][2]))
    o.write("%8s    %8s    %8f s\n"%("Recycle","d1",dic['acqus']['D'][1]))
    o.write("%8s    %8s    %8f s\n"%("Scans","NS",dic['acqus']['NS']))
    o.write("%8s    %8s    %8f us\n"%("delta","p53",dic['acqus']['P'][53]))
    o.write("%8s    %8s    %8f s\n"%("T_diff","d9",dic['acqus']['D'][9]))
    o.write("%8s    %8s    %s\n"%("options","ZGOPTNS",dic["acqus"]["ZGOPTNS"]))
    

print("Params output in parameters.txt")
    

