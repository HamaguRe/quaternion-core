# Version 0.6.0 (2025-10-04)

* Within the `from_euler_angles` function, I have removed the use of the `FloatConst` trait. Whilst the `FloatConst` trait was intended to catch out-of-range values via `debug_assert!`, I did not wish to employ this trait solely for debugging purposes and have therefore eliminated it.
* The `from_rotation_vector` function now performs normalisation internally. This enables correct conversion to a Versor even when the norm of the rotation_vector argument exceeds 2π.
* The argument for adjusting the vector's movement amount has been removed from the `rotate_a_to_b` function. Should you wish to adjust the vector's movement amount, please use the method described in the function documentation.
* By employing the sinc function, singularities arising when inputting a quaternion with a zero vector component to the `exp` and `pow` function have been eliminated.

# Version 0.5.4 (2025-05-11)

* Added documentation comment for the `to_axis_angle` function.
* The internal implementation of the `to_rotation_vector` function has been improved. Before the change, the accuracy of `x.asin()` was degraded near `x=1`, so `atan` is used to prevent the accuracy degradation.
* The implementation of the `rotate_a_to_b_shortest` function has been fundamentally changed so that it can be calculated with good accuracy no matter what the positional relationship between the vectors `a` and `b` is (I think this is pretty amazing...).

# Version 0.5.3 (2024-10-27)

* The implementation related to Euler angles was reviewed. The calculation results are the same as before.
* Improved function documentation.

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