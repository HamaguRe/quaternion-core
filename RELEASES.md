# Version 0.5.2 (2024-08-03)

* The polynomial coefficients of the `sinc` function (in `src/pfs`) were re-calculated. The number of terms remained the same and the maximum error went from 1.55e-15 to 2.22e-16.
* I created a `norm2` function (in `src/pfs`) and organized the implementation.

# Version 0.5.1 (2024-07-15)

* The implementation of `src/pfs`(private module) was reviewed.
* The calculation results are not affected by this change.

# Version 0.5.0 (2023-11-30)

* Added `serde-serialize` feature (See [Pull request #2](https://github.com/HamaguRe/quaternion-core/pull/2)).

# Version 0.4.3 (2023-11-13)

* Added `identity` function. I thought I could just write `(1.0, [0.0; 3])` myself, but it is still more convenient to call the `identity` function.
* Active use of the FMA instruction when computing dot products (when `fma` feature is enabled).

# Version 0.4.2 (2023-06-02)

* Fixed a bug in the `to_rotation_vector` function. With that, version `0.4.0` and `0.4.1` were yanked.

# Version 0.4.1 (2023-05-28)

* Replaced the original pythag() function with the standard .hypot() method (because it was slow, even though the calculation results were the same).

# Version 0.4.0 (2023-02-05)

* The norm of the rotation vector  calculated by the `to_rotation_vector` function is now `[0, π]`.
* The singularity at `theta=0` in the `to_rotation_vector` and `from_rotation_vector` functions has been removed.
* The behavior of the zero vector input to the `norm` function has been changed. Previous: zero-vector --> new: NaN.
* `rotate_a_to_b_param` function has been renamed to `rotate_a_to_b_shortest`.
* `rotate_a_to_b` and `rotate_a_to_b_shortest` functions now return `Option<Quaternion3<T>>`.
* `norm-sqrt` feature is added.