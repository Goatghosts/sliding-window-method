def int_to_hex(a: int):
    return int.to_bytes(a, 32, "big").hex().upper()


def get_uncompressed_public_key(x: int, y: int):
    return (b"\x04" + x.to_bytes(32, "big") + y.to_bytes(32, "big")).hex().upper()


def mod_inv(a, m):
    return pow(a, m - 2, m)


def homogeneous_to_affine(x, y, z, p):
    """
    Convert a point on an elliptic curve in homogeneous coordinates to affine coordinates.
    """

    try:
        z_inv = mod_inv(z, p)
        x = (x * z_inv) % p
        y = (y * z_inv) % p
    except:
        return 0, 0

    return x, y
