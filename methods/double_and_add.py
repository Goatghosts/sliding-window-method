from methods.secp256k1 import Secp256k1


class DoubleAndAdd(Secp256k1):
    """
    Multiply a point on an elliptic curve in homogeneous projective coordinates using double-and-add method.
    """

    def __init__(self) -> None:
        super().__init__()

    def double_and_add_method(self, x1, y1, z1, k):
        xt = x1
        yt = y1
        zt = z1
        x = 0
        y = 1
        z = 0
        for i in range(k.bit_length()):
            if k & (1 << i):
                x, y, z = self.add_points(x, y, z, xt, yt, zt)
            xt, yt, zt = self.double_point(xt, yt, zt)
        return x, y, z

    def multiplication(self, private_key):
        x, y, z = self.double_and_add_method(self.Gx, self.Gy, 1, private_key)
        return self.homogeneous_to_affine(x, y, z)
