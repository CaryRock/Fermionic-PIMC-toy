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
    parser.add_argument("-T", "--Tmax", type=float, help="Maximum T to go run")
    parser.add_argument("-log", "--logy", action='store_true', help="Use \
            numpy's semilogy to plot the results.")
    parser.add_argument("-c", "--cutoff", type=int, help="Cutoff temperature")
    parser.add_argument("-N", "--N", type=int, help="Number of particles")
    parser.add_argument("-p", "--production", type=pl.Path, help="File \
            containing the production code's values.")
    parser.add_argument("-hf", "--high-temp-fermion", action='store_true', 
            help="Plot high-temperature expansion of exact fermion energy.")
    parser.add_argument("-lf", "--low-temp-fermion", action='store_true',
            help="Plot low-temperature expansion of exact boson energy.")
    parser.add_argument("-hb", "--high-temp-boson", action='store_true',
            help="Plot high-temperature expansion of exact boson energy.")
    parser.add_argument("-lb", "--low-temp-boson", action='store_true',
            help="Plot low-temperature expansion of exact boson energy.")
    return parser

def main(argv=None):
    if argv is None:
        argv = sys.argv
    
    parser = create_parser()
    args = parser.parse_args(argv[1:])

    tmax = int(args.Tmax)
    production = args.production

    L = args.cutoff
    N = args.N

    E_f = N*N/2.0 * np.ones(100*tmax)
    E_b = N/2.0 * np.ones(100*tmax)
    E_fhT = np.copy(E_f)
    E_flT = np.copy(E_f)
    E_bhT = np.copy(E_b)
    E_blT = np.copy(E_b)

    temp = np.linspace(0.0001, tmax, 100*tmax)
    for t in range(len(temp)):
        for k in range(1,N+1):
            E_f[t] += k / (np.exp(k / temp[t]) - 1.0)
            E_b[t] += k / (np.exp(k / temp[t]) - 1.0)
            E_fhT[t] += 1.0 / (1.0/temp[t] * (1.0 + k/(2.0 * temp[t])))
            E_bhT[t] += 1.0 / (1.0/temp[t] * (1.0 + k/(2.0 * temp[t])))
        E_flT[t] += 1.0 / (np.exp(1.0 / temp[t]) - 1)
        E_blT[t] += 1.0 / (np.exp(1.0 / temp[t]) - 1)

    if args.production is not None:
        with open(production, "r") as pro:
            prod_data = np.loadtxt(pro)

        prod_temp = prod_data[:-1, 0]
        production = prod_data[:-1,1]

    if not args.logy:
        plt.plot(temp, E_b, 'b', label="Boson <E>, N = 2")
        plt.plot(0, N/2, 'b.')
        plt.plot(temp, E_f, 'r', label="Fermion <E>, N = 2")
        plt.plot(0, N*N/2, 'r.')

        if args.high_temp_fermion:
            plt.fill_between(temp, E_fhT, E_f, alpha=0.33)
        if args.low_temp_fermion:
            plt.fill_between(temp, E_flT, E_f, alpha=0.5)
        if args.high_temp_boson:
            plt.fill_between(temp, E_bhT, E_b, alpha=0.33)
        if args.low_temp_boson:
            plt.fill_between(temp, E_blT, E_b, alpha=0.5)

        if args.production is not None:
            plt.plot(prod_temp, production, 'gx', label="Prod. Boson <E>, N = 2")

    else:
        plt.semilogy(temps, bosons, 'b.')
        plt.semilogy(temps, fermions, 'r.')
    
    plt.axvline(12.57, color='blue', marker="|", linestyle="dashed")
    plt.axvline(2.0, color='red', marker="|", linestyle="dashed")
    plt.xlabel("T (K)")
    plt.ylabel("Energy (K)")
    #plt.axvline(1.5, color="red")
    #plt.axvline(2.5, color="blue")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
