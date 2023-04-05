# Sliding Window Method

## [Original](https://github.com/stephancill/msm-hardware-acceleration/blob/master/sim/ecc.py)

---

Thanks to [stephancill](https://github.com/stephancill), I finally found a decent implementation of the sliding window method as described on the [Wikipedia page](https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication#Sliding-window_method).

This code demonstrates the process of multiplying the base point of the secp256k1 curve by a scalar (private key) and converts the resulting point into the public key format.

At present, this code only compares the results of multiplying the same scalar using the double and add method, the sliding window method and gtable method.