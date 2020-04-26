import sys


def list_to_fa(data, col, outfile):
    n = []
    with open(outfile, "w") as f:
        for i, d in enumerate(data):
            id = ">M_" + str(i + 1)
            f.write(id + "\n")
            n.append(id)
            f.write(d[col] + "\n")
    return n


if __name__ == "__main__":
    c = -1
    with open(sys.argv[1]) as f:
        for line in f:
            c += 1
            if c == 0:
                continue
            w = line.strip("\t").split(",")
            print(">M_" + str(c))
            print(w[1])
