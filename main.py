import random

from utils import int_to_hex, get_uncompressed_public_key, homogeneous_to_affine


def add_points(x1, y1, z1, x2, y2, z2, a, b, p):
    """
    Add two points on an elliptic curve in homogeneous coordinates.
    Source: https://eprint.iacr.org/2015/1060 (alg 1)
    """

    b3 = 3 * b
    t0 = (x1 * x2) % p
    t1 = (y1 * y2) % p
    t2 = (z1 * z2) % p
    t3 = (x1 + y1) % p
    t4 = (x2 + y2) % p
    t3 = (t3 * t4) % p
    t4 = (t0 + t1) % p
    t3 = (t3 - t4) % p
    t4 = (x1 + z1) % p
    t5 = (x2 + z2) % p
    t4 = (t4 * t5) % p
    t5 = (t0 + t2) % p
    t4 = (t4 - t5) % p
    t5 = (y1 + z1) % p
    x3 = (y2 + z2) % p
    t5 = (t5 * x3) % p
    x3 = (t1 + t2) % p
    t5 = (t5 - x3) % p
    z3 = (a * t4) % p
    x3 = (b3 * t2) % p
    z3 = (x3 + z3) % p
    x3 = (t1 - z3) % p
    z3 = (t1 + z3) % p
    y3 = (x3 * z3) % p
    t1 = (t0 + t0) % p
    t1 = (t1 + t0) % p
    t2 = (a * t2) % p
    t4 = (b3 * t4) % p
    t1 = (t1 + t2) % p
    t2 = (t0 - t2) % p
    t2 = (a * t2) % p
    t4 = (t4 + t2) % p
    t0 = (t1 * t4) % p
    y3 = (y3 + t0) % p
    t0 = (t5 * t4) % p
    x3 = (t3 * x3) % p
    x3 = (x3 - t0) % p
    t0 = (t3 * t1) % p
    z3 = (t5 * z3) % p
    z3 = (z3 + t0) % p
    return x3, y3, z3


def double_point(x, y, z, a, b, p):
    """
    Double a point on an elliptic curve in homogeneous projective coordinates.
    Source: https://eprint.iacr.org/2015/1060 (alg 3)
    """
    b3 = 3 * b
    t0 = (x * x) % p
    t1 = (y * y) % p
    t2 = (z * z) % p
    t3 = (x * y) % p
    t3 = (t3 + t3) % p
    z3 = (x * z) % p
    z3 = (z3 + z3) % p
    x3 = (a * z3) % p
    y3 = (b3 * t2) % p
    y3 = (x3 + y3) % p
    x3 = (t1 - y3) % p
    y3 = (t1 + y3) % p
    y3 = (x3 * y3) % p
    x3 = (t3 * x3) % p
    z3 = (b3 * z3) % p
    t2 = (a * t2) % p
    t3 = (t0 - t2) % p
    t3 = (a * t3) % p
    t3 = (t3 + z3) % p
    z3 = (t0 + t0) % p
    t0 = (z3 + t0) % p
    t0 = (t0 + t2) % p
    t0 = (t0 * t3) % p
    y3 = (y3 + t0) % p
    t2 = (y * z) % p
    t2 = (t2 + t2) % p
    t0 = (t2 * t3) % p
    x3 = (x3 - t0) % p
    z3 = (t2 * t1) % p
    z3 = (z3 + z3) % p
    z3 = (z3 + z3) % p

    return x3, y3, z3


def double_and_add_method(x1, y1, z1, k, a, b, p):
    """
    Multiply a point on an elliptic curve in homogeneous projective coordinates using double-and-add method.
    """

    xt = x1
    yt = y1
    zt = z1

    x = 0
    y = 1
    z = 0

    for i in range(k.bit_length()):
        if k & (1 << i):
            x, y, z = add_points(x, y, z, xt, yt, zt, a, b, p)
        xt, yt, zt = double_point(xt, yt, zt, a, b, p)

    return x, y, z


def sliding_window_method(x1, y1, z1, k, a, b, p):
    """
    Compute multiplication using sliding window method
    """

    x = 0
    y = 1
    z = 0

    w = 4
    m = k.bit_length()

    tP = []
    for i in range(2**w):
        tP.append(double_and_add_method(x1, y1, z1, i, a, b, p))

    i = m

    while i > 0:
        di = k & (1 << (i - 1))
        if di == 0:
            x, y, z = double_point(x, y, z, a, b, p)
            i -= 1
        else:
            j = min(w, i)
            mask = (2**j - 1) << (i - j)
            t = (k & mask) >> (i - j)
            i -= j
            for _ in range(j):
                x, y, z = double_point(x, y, z, a, b, p)
            tPx, tPy, tPz = tP[t]
            x, y, z = add_points(x, y, z, tPx, tPy, tPz, a, b, p)

    return x, y, z


if __name__ == "__main__":
    # Secp256k1 params
    p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
    a = 0x0000000000000000000000000000000000000000000000000000000000000000
    b = 0x0000000000000000000000000000000000000000000000000000000000000007
    Gx = int(0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798)
    Gy = int(0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)

    for _ in range(100):
        private_key = random.randrange(1, 2**256)
        x1, y1, z1 = double_and_add_method(Gx, Gy, 1, private_key, a, b, p)
        ax1, ay1 = homogeneous_to_affine(x1, y1, z1, p)
        x2, y2, z2 = sliding_window_method(Gx, Gy, 1, private_key, a, b, p)
        ax2, ay2 = homogeneous_to_affine(x2, y2, z2, p)
        assert (ax1 == ax2) and (ay1 == ay2), "The x and y values obtained using different methods are not equal."
        print("Private key (HEX):", int_to_hex(private_key))
        print("Uncompressed public key (double and add method):", get_uncompressed_public_key(ax1, ay1))
        print("Uncompressed public key (sliding window method):", get_uncompressed_public_key(ax2, ay2))
        print()
