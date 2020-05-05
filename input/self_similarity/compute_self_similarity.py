import math


def read_K1(mfile="BLOSUM62-2.matrix.txt"):
    bd = {}
    colid = []
    rowid = []

    c = 0
    with open(mfile) as f:
        for line in f:
            c += 1
            if c == 1:
                colid = line.strip("\n").split(" ")
                continue
            w = line.strip("\n").split(" ")
            id = w[0]
            v = [float(x) for x in w[1:]]
            rowid.append(id)
            x = bd.get(id, {})
            for i, vi in enumerate(v):
                x[colid[i]] = vi
            bd[id] = x

    # print bd["A"]["Y"]
    # print bd["Y"]["A"]
    # print bd["C"]["A"], bd["A"]["C"]
    # print bd["Y"]["C"], bd["C"]["Y"]

    BLOSUM62_2 = bd
    K1 = {}
    beta = 0.11387
    for i in list(BLOSUM62_2.keys()):
        x = K1.get(i, {})
        for j in list(BLOSUM62_2[i].keys()):
            x[j] = math.pow(BLOSUM62_2[i][j], beta)
        K1[i] = x
    return K1


def compute_K2k(u, v, K1):
    if len(u) != len(v):
        return None
    k = len(u)
    p = K1[u[0]][v[0]]
    for i in range(1, k):
        p = p * K1[u[i]][v[i]]
    return p


def compute_K3(f, g, K1):
    max_k = min(len(f), len(g))
    s = 0
    for k in range(1, max_k + 1):
        for i in range(len(f) - (k - 1)):
            u = f[i:i + k]
            for j in range(len(g) - (k - 1)):
                v = g[j:j + k]
                s += compute_K2k(u, v, K1)
    return s


def compute_K_hat_3(x, y, K1):
    return compute_K3(x, y, K1) / math.sqrt(compute_K3(x, x, K1) * compute_K3(y, y, K1))


class selfsim():
    def __init__(self, mfile):
        self.K1 = read_K1(mfile)

    def compute_K_hat_3(self, x, y):  # K^3
        return compute_K_hat_3(x, y, self.K1)


if __name__ == "__main__":
    K1 = read_K1()
    print(compute_K_hat_3("AAAAA", "AAAAA", K1))
    for i in range(5):
        print(i, compute_K_hat_3("AAAAA", "WWWWW" * (i + 1), K1))
    for i in list(K1.keys()):
        print(i, compute_K_hat_3("AAAAA", "AA" + i + "AA", K1))
    print()
    s = selfsim("BLOSUM62-2.matrix.txt")
    print(s.compute_K_hat_3("AAAAA", "AAAAA"))
