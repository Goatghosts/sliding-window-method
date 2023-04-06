from utils import mod_inv


class Secp256k1:
    """
    Main Secp256k1 params and methods
    """

    def __init__(self) -> None:
        self.p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
        self.a = 0x0000000000000000000000000000000000000000000000000000000000000000
        self.b = 0x0000000000000000000000000000000000000000000000000000000000000007
        self.Gx = int(0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798)
        self.Gy = int(0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)

    def add_points(self, x1, y1, z1, x2, y2, z2):
        """
        Add two points on an elliptic curve in homogeneous coordinates.
        Source: https://eprint.iacr.org/2015/1060 (alg 1)
        """
        b3 = 3 * self.b
        t0 = (x1 * x2) % self.p
        t1 = (y1 * y2) % self.p
        t2 = (z1 * z2) % self.p
        t3 = (x1 + y1) % self.p
        t4 = (x2 + y2) % self.p
        t3 = (t3 * t4) % self.p
        t4 = (t0 + t1) % self.p
        t3 = (t3 - t4) % self.p
        t4 = (x1 + z1) % self.p
        t5 = (x2 + z2) % self.p
        t4 = (t4 * t5) % self.p
        t5 = (t0 + t2) % self.p
        t4 = (t4 - t5) % self.p
        t5 = (y1 + z1) % self.p
        x3 = (y2 + z2) % self.p
        t5 = (t5 * x3) % self.p
        x3 = (t1 + t2) % self.p
        t5 = (t5 - x3) % self.p
        z3 = (self.a * t4) % self.p
        x3 = (b3 * t2) % self.p
        z3 = (x3 + z3) % self.p
        x3 = (t1 - z3) % self.p
        z3 = (t1 + z3) % self.p
        y3 = (x3 * z3) % self.p
        t1 = (t0 + t0) % self.p
        t1 = (t1 + t0) % self.p
        t2 = (self.a * t2) % self.p
        t4 = (b3 * t4) % self.p
        t1 = (t1 + t2) % self.p
        t2 = (t0 - t2) % self.p
        t2 = (self.a * t2) % self.p
        t4 = (t4 + t2) % self.p
        t0 = (t1 * t4) % self.p
        y3 = (y3 + t0) % self.p
        t0 = (t5 * t4) % self.p
        x3 = (t3 * x3) % self.p
        x3 = (x3 - t0) % self.p
        t0 = (t3 * t1) % self.p
        z3 = (t5 * z3) % self.p
        z3 = (z3 + t0) % self.p
        return x3, y3, z3

    def double_point(self, x, y, z):
        """
        Double a point on an elliptic curve in homogeneous projective coordinates.
        Source: https://eprint.iacr.org/2015/1060 (alg 3)
        """
        b3 = 3 * self.b
        t0 = (x * x) % self.p
        t1 = (y * y) % self.p
        t2 = (z * z) % self.p
        t3 = (x * y) % self.p
        t3 = (t3 + t3) % self.p
        z3 = (x * z) % self.p
        z3 = (z3 + z3) % self.p
        x3 = (self.a * z3) % self.p
        y3 = (b3 * t2) % self.p
        y3 = (x3 + y3) % self.p
        x3 = (t1 - y3) % self.p
        y3 = (t1 + y3) % self.p
        y3 = (x3 * y3) % self.p
        x3 = (t3 * x3) % self.p
        z3 = (b3 * z3) % self.p
        t2 = (self.a * t2) % self.p
        t3 = (t0 - t2) % self.p
        t3 = (self.a * t3) % self.p
        t3 = (t3 + z3) % self.p
        z3 = (t0 + t0) % self.p
        t0 = (z3 + t0) % self.p
        t0 = (t0 + t2) % self.p
        t0 = (t0 * t3) % self.p
        y3 = (y3 + t0) % self.p
        t2 = (y * z) % self.p
        t2 = (t2 + t2) % self.p
        t0 = (t2 * t3) % self.p
        x3 = (x3 - t0) % self.p
        z3 = (t2 * t1) % self.p
        z3 = (z3 + z3) % self.p
        z3 = (z3 + z3) % self.p
        return x3, y3, z3

    # Функция для сложения двух точек в проективных координатах
    def add_points_projective(self, x1, y1, z1, x2, y2, z2):
        t0 = (y1 * z2) % self.p
        t1 = (y2 * z1) % self.p
        u0 = (x1 * z2) % self.p
        u1 = (x2 * z1) % self.p
        t = (t0 - t1) % self.p
        u = (u0 - u1) % self.p
        u2 = (u * u) % self.p
        v = (z1 * z2) % self.p
        w = (t * t * v - u2 * (u0 + u1)) % self.p
        u3 = (u * u2) % self.p
        rx = (u * w) % self.p
        ry = (t * (u0 * u2 - w) - t0 * u3) % self.p
        rz = (u3 * v) % self.p
        return rx, ry, rz

    # Функция для удвоения точки в проективных координатах
    def double_point_projective(self, x, y, z):
        t = (x * x * 3 + self.a * z * z) % self.p
        u = (y * z * 2) % self.p
        v = (u * x * y * 2) % self.p
        w = (t * t - v * 2) % self.p
        rx = (u * w) % self.p
        ry = (t * (v - w) - u * u * y * y * 2) % self.p
        rz = (u * u * u) % self.p
        return rx, ry, rz

    def double_point_affine(self, x, y):
        s = ((3 * x * x + self.a) * mod_inv(2 * y, self.p)) % self.p
        x3 = (s * s - 2 * x) % self.p
        y3 = (s * (x - x3) - y) % self.p
        return x3, y3

    def add_points_affine(self, x1, y1, x2, y2):
        s = ((y2 - y1) * mod_inv(x2 - x1, self.p)) % self.p
        x3 = (s * s - x1 - x2) % self.p
        y3 = (s * (x1 - x3) - y1) % self.p
        return x3, y3

    def point_subtract_affine(self, x1, y1, x2, y2):
        return self.add_points_affine(x1, y1, x2, -y2 % self.p)

    def is_on_secp256k1_affine(self, x, y):
        return (y**2 - x**3 - self.a * x - self.b) % self.p == 0

    def find_y_pairs_affine(self, x):
        y_squared = (x**3 + self.a * x + self.b) % self.p
        y1 = pow(y_squared, (self.p + 1) // 4, self.p)
        if (y1**2) % self.p != y_squared:
            raise ValueError("No y value for this x on this curve")
        y2 = self.p - y1
        return y1, y2

    def homogeneous_to_affine(self, x, y, z):
        """
        Convert a point on an elliptic curve in homogeneous coordinates to affine coordinates.
        """
        try:
            z_inv = mod_inv(z, self.p)
            x = (x * z_inv) % self.p
            y = (y * z_inv) % self.p
        except:
            return 0, 0

        return x, y
