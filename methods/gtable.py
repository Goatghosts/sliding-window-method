import pickle

from methods.secp256k1 import Secp256k1


class GtableMethod(Secp256k1):
    """
    Compute multiplication using gtables
    """

    def __init__(self, init_table=True) -> None:
        super().__init__()
        self.NUM_GTABLE_CHUNK = 16
        self.NUM_GTABLE_VALUE = 65536
        self.CHUNK_FIRST_ELEMENT = [65536 * i for i in range(16)]
        self.g_table_x = [0] * self.NUM_GTABLE_CHUNK * self.NUM_GTABLE_VALUE
        self.g_table_y = [0] * self.NUM_GTABLE_CHUNK * self.NUM_GTABLE_VALUE
        if init_table:
            self.load_g_table()
            with open("g_table_data.pkl", "wb") as f:
                pickle.dump((self.g_table_x, self.g_table_y), f)
        else:
            with open("g_table_data.pkl", "rb") as f:
                self.g_table_x, self.g_table_y = pickle.load(f)

    def get_g_table(self):
        GTable = [(0, 0)] * self.NUM_GTABLE_CHUNK * self.NUM_GTABLE_VALUE
        N = (self.Gx, self.Gy)
        for i in range(self.NUM_GTABLE_CHUNK):
            GTable[i * self.NUM_GTABLE_VALUE] = N
            N = self.double_point_affine(N[0], N[1])
            for j in range(1, self.NUM_GTABLE_VALUE - 1):
                GTable[(i * self.NUM_GTABLE_VALUE) + j] = N
                point = GTable[i * self.NUM_GTABLE_VALUE]
                N = self.add_points_affine(N[0], N[1], point[0], point[1])
        return GTable

    def load_g_table(self):
        g_table = self.get_g_table()
        for i in range(self.NUM_GTABLE_CHUNK):
            for j in range(self.NUM_GTABLE_VALUE - 1):
                element = (i * self.NUM_GTABLE_VALUE) + j
                p = g_table[element]
                self.g_table_x[element] = p[0]
                self.g_table_y[element] = p[1]
        g_table.clear()

    def get_16_bits_parts(self, number):
        parts = [0] * 16
        for i in range(16):
            parts[i] = number & 0xFFFF
            number >>= 16
        return parts

    def g_table_method(self, privKey):
        converted_private = self.get_16_bits_parts(privKey)
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
                qx, qy, qz = self.add_points(qx, qy, qz, gx, gy, 1)

        return qx, qy, qz

    def scalar_mult(self, private_key):
        x, y, z = self.g_table_method(private_key)
        return self.homogeneous_to_affine(x, y, z)
