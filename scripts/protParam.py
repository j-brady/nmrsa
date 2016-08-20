#!/usr/bin/python
import argparse
import pprint
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

# taken from Gill and von Hipple Anal Biochem 182,319-326 (1989) 
extinction_coeffs = {"W":5690,"Y":1280,"C":120}

def protParam(seq):
    params = ProteinAnalysis(seq)
    mw = params.molecular_weight()
    c_aa = params.count_amino_acids()
    p_aa = params.get_amino_acids_percent()
    gravy = params.gravy()
    aromaticity = params.aromaticity()
    isoelectric_point = params.isoelectric_point()
    ext_coeff = sum([c_aa["W"]*5690,c_aa["Y"]*1280,c_aa["C"]*120])
    mgml = ext_coeff * (1./mw)
    
    print("Amino acid count")
    pprint.pprint(c_aa)
    print("Amino acid percent")
    pprint.pprint(p_aa)
    print("Molecular weight")
    print("%f Da"%mw)
    print("Gravy")
    print(gravy)
    print("Isoelectric point")
    print(isoelectric_point)
    print("Aromaticity")
    print(aromaticity)
    print("Extinction coefficient: %d M-1cm-1 (Assuming reduced)"%ext_coeff)
    print("")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Jazzy J's Protparams")

    parser.add_argument("-s","--sequence",type=str)
    parser.add_argument("-f","--fasta", action='store_true')
   
    args = parser.parse_args()
    if args.fasta:
        seq = SeqIO.read(args.sequence,"fasta")
        print seq.name
        seq = str(seq.seq)
        print seq
        protParam(seq)



