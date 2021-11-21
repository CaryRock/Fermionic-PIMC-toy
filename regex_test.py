import numpy as np
import array as arr
import re
import sys

def parse_file(name):
    with open(name, "r") as f:
        raw_read = f.read().splitlines()[3:]
    f.close()
    
    N = len(raw_read)
    M = len(raw_read[0])
    print(f"N, M = {N}, {M}")

    for i in range(N):
        raw_read[i] = re.sub(r'\s+',',',raw_read[i])

    M = len(raw_read[0])
    print(f"N, M == {N}, {M}")
    #data = np.zeros(int(M/15))

    out = ['']*N
    for i in range(N):
        out[i] = raw_read[i].split(',')
        out[i] = out[i][1:]

    L = len(out[0])
    O = len(out)
    data = np.zeros(L * O)
    for i in range(O):
        for j in range(L):
            data[O*i + j] = float(out[i][j])

    return data

def main():
    name = sys.argv[1]
    
    production = ((re.search(r'gce-[a-zA-Z0-9-\.+]*',name) != None))
    if production:
        data = parse_file(name)
    else:
        data = 

if __name__ == "__main__":
    main()
