def int_to_hex(a: int):
    return int.to_bytes(a, 32, "big").hex().upper()


def get_uncompressed_public_key(x: int, y: int):
    return (b"\x04" + x.to_bytes(32, "big") + y.to_bytes(32, "big")).hex().upper()


def mod_inv(a, m):
    return pow(a, m - 2, m)


def get_8_bits_parts(number):
    parts = [0] * 32
    for i in range(8):
        parts[i] = number & 0xFF
        number >>= 8
    return parts


def get_16_bits_parts(number):
    parts = [0] * 16
    for i in range(16):
        parts[i] = number & 0xFFFF
        number >>= 16
    return parts
