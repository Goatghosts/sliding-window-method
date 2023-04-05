from methods.secp256k1 import Secp256k1
from methods.double_and_add import DoubleAndAdd


class SlidingWindow(Secp256k1):
    """
    Compute multiplication using sliding window method
    """

    def __init__(self, w) -> None:
        super().__init__()
        self.w = w
        self.tP = []
        double_and_add = DoubleAndAdd()
        for i in range(2**w):
            self.tP.append(double_and_add.double_and_add_method(self.Gx, self.Gy, 1, i))

    def sliding_window_method(self, x1, y1, z1, k):
        x = 0
        y = 1
        z = 0
        i = k.bit_length()
        while i > 0:
            di = k & (1 << (i - 1))
            if di == 0:
                x, y, z = self.double_point(x, y, z)
                i -= 1
            else:
                j = min(self.w, i)
                mask = (2**j - 1) << (i - j)
                t = (k & mask) >> (i - j)
                i -= j
                for _ in range(j):
                    x, y, z = self.double_point(x, y, z)
                tPx, tPy, tPz = self.tP[t]
                x, y, z = self.add_points(x, y, z, tPx, tPy, tPz)

        return x, y, z

    def scalar_mult(self, private_key):
        x, y, z = self.sliding_window_method(self.Gx, self.Gy, 1, private_key)
        return self.homogeneous_to_affine(x, y, z)
