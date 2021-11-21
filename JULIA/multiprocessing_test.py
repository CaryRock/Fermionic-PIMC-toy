import sys
import multiprocessing as mp
import concurrent.futures
from time import perf_counter
import numpy as np
import matplotlib.pyplot as plt
import mpmath as mpm

def compute(y,x):
    # Analyze some portion of the data given to it
    if (len(y) != len(x)): exit()
    for i in range(len(x)):
        y[i] = mpm.cos(x[i])
    return y

def setY(yy,y):
    for i in range(len(yy)):
        for j in range(yy[0]):
            y[len(i*len(yy[0]) + j)] = yy[i][j]
    return y

def main():
    if sys.argv[1] == '':
        numThreads = 11
    else:
        numThreads = int(sys.argv[1])
    numThreads = 11
    div = 1e-6
    xa = -5.0
    xb = 5.0
    
    comp_start = perf_counter()
    x = np.arange(xa, xb, div)
    y = np.zeros_like(x)
    N = len(x)
    n = int(N/numThreads)
    rem = N % numThreads    # Because the last little bit might not be captured

    xx = np.zeros([numThreads,n])
    for i in range(numThreads):
        for j in range(n):
            xx[i][j] = x[i*n + j]
    yy = np.zeros_like(xx)
    
    fig = plt.figure()
    mult_start = perf_counter()
    '''
    # Step 1: Init multiprocessing.Pool()
    pool = mp.Pool(numThreads)  #mp.Pool(int(mp.cpu_count()/2))
    # Step 2: 'pool.apply' the function
    for k in range(numThreads):
        yy[k] = pool.apply(compute, args=(yy[k], xx[k]))
#        pool.starmap(compute, [(yy[k:k+1], xx[k:k+1])])
#        pool.starmap(setY, [(yy[k:k+1], y[k:k+1])])
#        pool.apply(plt.plot, args=(yy[k], xx[k], 'b.'))
#    y = pool.apply(setY, args=(yy, y))

    # Step 3: Pool's closed
    pool.close()
    '''
#    zz = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        zz = executor.map(compute, yy, xx)
    
    zz = list(zz)
    print("All threads threaded.")
    mult_end = perf_counter()
    print(f"It took {mult_end - mult_start :0.2f} seconds to do the threading.")
    for l in range(numThreads):
        plt.plot(xx[l], zz[l],'b.')

    comp_end = perf_counter()
    print(f"It took {comp_end - comp_start :0.2f} seconds to do the computing \
and plotting")
    plt.show()

if __name__ == "__main__":
    main()
