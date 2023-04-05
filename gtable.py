import pickle

# Параметры кривой secp256k1
p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
a = 0x0000000000000000000000000000000000000000000000000000000000000000
b = 0x0000000000000000000000000000000000000000000000000000000000000007
Gx = int(0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798)
Gy = int(0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141


NUM_GTABLE_CHUNK = 16
NUM_GTABLE_VALUE = 65536
SIZE_GTABLE_POINT = 32
CHUNK_FIRST_ELEMENT = [65536 * i for i in range(16)]
mod_inv = lambda a, m: pow(a, m - 2, m)


def DoubleDirect(x, y):
    s = ((3 * x * x + a) * mod_inv(2 * y, p)) % p
    x3 = (s * s - 2 * x) % p
    y3 = (s * (x - x3) - y) % p
    return (x3, y3)


def AddDirect(x1, y1, x2, y2):
    s = ((y2 - y1) * mod_inv(x2 - x1, p)) % p
    x3 = (s * s - x1 - x2) % p
    y3 = (s * (x1 - x3) - y1) % p
    return (x3, y3)


def get_g_table():
    GTable = [(0, 0)] * NUM_GTABLE_CHUNK * NUM_GTABLE_VALUE
    N = (Gx, Gy)
    for i in range(NUM_GTABLE_CHUNK):
        GTable[i * NUM_GTABLE_VALUE] = N
        N = DoubleDirect(N[0], N[1])
        for j in range(1, NUM_GTABLE_VALUE - 1):
            GTable[(i * NUM_GTABLE_VALUE) + j] = N
            point = GTable[i * NUM_GTABLE_VALUE]
            N = AddDirect(N[0], N[1], point[0], point[1])
    return GTable


def load_g_table(g_table_x, g_table_y):
    print("loadGTable started")

    g_table = get_g_table()

    for i in range(NUM_GTABLE_CHUNK):
        for j in range(NUM_GTABLE_VALUE - 1):
            element = (i * NUM_GTABLE_VALUE) + j
            p = g_table[element]
            g_table_x[element] = p[0]
            g_table_y[element] = p[1]

    print("loadGTable finished!")


def get_16_bits_parts(number):
    parts = [0] * 16
    for i in range(16):
        parts[i] = number & 0xFFFF
        number >>= 16
    return parts


def _PointMultiSecp256k1(privKey, g_table_x, g_table_y):
    converted_private = get_16_bits_parts(privKey)
    for chunk in range(NUM_GTABLE_CHUNK):
        if converted_private[chunk] > 0:
            index = CHUNK_FIRST_ELEMENT[chunk] + (converted_private[chunk] - 1)
            qx = g_table_x[index]
            qy = g_table_y[index]
            chunk += 1
            break

    for chunk in range(chunk, NUM_GTABLE_CHUNK):
        if converted_private[chunk] > 0:
            index = CHUNK_FIRST_ELEMENT[chunk] + (converted_private[chunk] - 1)
            gx = g_table_x[index]
            gy = g_table_y[index]
            qx, qy = AddDirect(qx, qy, gx, gy)

    print((b"\x04" + qx.to_bytes(32, "big") + qy.to_bytes(32, "big")).hex().upper())


# g_table_x = [0] * NUM_GTABLE_CHUNK * NUM_GTABLE_VALUE
# g_table_y = [0] * NUM_GTABLE_CHUNK * NUM_GTABLE_VALUE
# load_g_table(g_table_x, g_table_y)
# with open("g_table_data.pkl", "wb") as f:
#     pickle.dump((g_table_x, g_table_y), f)

# Загрузка данных из файла
with open("g_table_data.pkl", "rb") as f:
    g_table_x, g_table_y = pickle.load(f)


def int_to_hex(a: int):
    return int.to_bytes(a, 32, "big").hex().upper()


privKey = 1231231241343234123
print(int_to_hex(privKey))
_PointMultiSecp256k1(privKey, g_table_x, g_table_y)
