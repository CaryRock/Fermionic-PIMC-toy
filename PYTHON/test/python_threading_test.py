import threading
from time import perf_counter
import numpy as np
import matplotlib.pyplot as plt
import mpmath as mpm

class myThread(threading.Thread):
    def __init__(self, threadID, arraySliceY, arraySliceX):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.arraySliceY = arraySliceY
        self.arraySliceX = arraySliceX

def run(y,x):
    # Analyze some portion of the data given to it
    for i in range(len(x)):
        y[i] = mpm.cos(x[i])


def main():
    div = 1e-6
    mul = 1/div
    xa = -5.0
    xb = 5.0
    x = np.arange(xa, xb, div)
    y = np.zeros_like(x)
    N = len(x)
    numThreads = 11
    n = int(N/numThreads)
    rem = N % numThreads


    threads = [threading.Thread(target=run,args=(y[int(i*(n+1)):int((i+1)*(n+1))],x[int(i*(n+1)):int((i+1)*(n+1))])) for i in range(numThreads)]

    for i in range(numThreads):
        print(f"Ranges: {int(i*(n+1))}:{int((i+1)*(n+1))}")

    for thread in threads:
        thread.start()

    for thread in threads: 
        thread.join()

    print("All threads threaded.")

    fig = plt.figure()
    plt.plot(x, y,'b.')
    plt.show()

if __name__ == "__main__":
    start_time = perf_counter()
    main()
    end_time = perf_counter()
    print(f"It took {end_time - start_time :0.2f} seconds to complete.")
