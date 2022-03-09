import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse as ap
import pathlib as pl

def create_parser():
    parser = ap.ArgumentParser(description='Recursive partition function \
            solver. Computes desired quantities from the PF. Designed with \
            using GNU paralell in mind and storing output "manually".',
            formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-b", "--bosons", type=pl.Path, help="File containing \
            boson data.")
    parser.add_argument("-f", "--fermions", type=pl.Path, help="File \
            containing fermion data.")
    parser.add_argument("-log", "--logy", action='store_true', help="Use \
            numpy's semilogy to plot the results.")

    return parser

def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = create_parser()
    args = parser.parse_args(argv[1:])

    bosons = args.bosons
    fermions = args.fermions

    with open(bosons, "r") as bos:
        boson_data = np.loadtxt(bos)
    
    with open(fermions, "r") as fer:
        fermion_data = np.loadtxt(fer)
    
    bosons = boson_data[:,1]
    fermions = fermion_data[:,1]

    for i in range(len(bosons)):
        print(fermions[i] - bosons[i])
    
    temps = boson_data[:,0]#np.linspace(1, len(boson_data), len(boson_data))
    if not args.logy:
        plt.plot(temps, bosons, 'b.', label="Boson <E>, N = 2")
        plt.plot(temps, fermions, 'r.', label="Fermion <E>, N = 2")

    else:
        plt.semilogy(temps, bosons, 'b.')
        plt.semilogy(temps, fermions, 'r.')

    plt.plot(temps, np.abs(fermions - bosons), 'k.', label="Difference")
    plt.xlabel("T (K)")
    plt.ylabel("Energy (K)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
