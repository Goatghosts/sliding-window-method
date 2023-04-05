import random

from methods.double_and_add import DoubleAndAdd
from methods.sliding_window import SlidingWindow
from utils import int_to_hex, get_uncompressed_public_key


if __name__ == "__main__":
    double_and_add = DoubleAndAdd()
    sliding_window = SlidingWindow(4)

    for _ in range(100):
        private_key = random.randrange(1, 2**256)
        x1, y1 = double_and_add.multiplication(private_key)
        x2, y2 = sliding_window.multiplication(private_key)
        assert (x1 == x2) and (y1 == y2), "The x and y values obtained using different methods are not equal."
        print("Private key (HEX):", int_to_hex(private_key))
        print("Uncompressed public key (double and add method):", get_uncompressed_public_key(x1, y1))
        print("Uncompressed public key (sliding window method):", get_uncompressed_public_key(x2, y2))
        print()
