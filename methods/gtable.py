import pickle

from utils import get_8_bits_parts, get_16_bits_parts
from methods.secp256k1 import Secp256k1


class GtableMethod(Secp256k1):
    """
    Compute multiplication using gtables
    """

    def __init__(self, init_table=True) -> None:
        super().__init__()
        self.NUM_GTABLE_CHUNK = 32
        self.NUM_GTABLE_VALUE = 256
        self.CHUNK_FIRST_ELEMENT = [self.NUM_GTABLE_VALUE * i for i in range(self.NUM_GTABLE_CHUNK)]
        self.g_table_x = [0] * self.NUM_GTABLE_CHUNK * self.NUM_GTABLE_VALUE
        self.g_table_y = [0] * self.NUM_GTABLE_CHUNK * self.NUM_GTABLE_VALUE
        if init_table:
            self.get_g_table()
            with open("g_table_data.pkl", "wb") as f:
                pickle.dump((self.g_table_x, self.g_table_y), f)
        else:
            with open("g_table_data.pkl", "rb") as f:
                self.g_table_x, self.g_table_y = pickle.load(f)

    def get_g_table(self):
        N = (self.Gx, self.Gy)
        for i in range(self.NUM_GTABLE_CHUNK):
            element = i * self.NUM_GTABLE_VALUE
            self.g_table_x[element] = N[0]
            self.g_table_y[element] = N[1]
            point = N
            N = self.double_point_affine(N[0], N[1])
            for j in range(1, self.NUM_GTABLE_VALUE - 1):
                element = (i * self.NUM_GTABLE_VALUE) + j
                self.g_table_x[element] = N[0]
                self.g_table_y[element] = N[1]
                N = self.add_points_affine(N[0], N[1], point[0], point[1])

    def g_table_method(self, privKey):
        converted_private = get_8_bits_parts(privKey)
        for chunk in range(self.NUM_GTABLE_CHUNK):
            if converted_private[chunk] > 0:
                index = self.CHUNK_FIRST_ELEMENT[chunk] + (converted_private[chunk] - 1)
                qx = self.g_table_x[index]
                qy = self.g_table_y[index]
                chunk += 1
                break

        qz = 1
        for chunk in range(chunk, self.NUM_GTABLE_CHUNK):
            if converted_private[chunk] > 0:
                index = self.CHUNK_FIRST_ELEMENT[chunk] + (converted_private[chunk] - 1)
                gx = self.g_table_x[index]
                gy = self.g_table_y[index]
                qx, qy, qz = self.add_points_projective(qx, qy, qz, gx, gy)

        return qx, qy, qz

    def scalar_mult(self, private_key):
        x, y, z = self.g_table_method(private_key)
        return self.homogeneous_to_affine(x, y, z)
