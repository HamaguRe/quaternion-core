//! Quaternion library written in Rust.
//! 
//! This provides Quaternion operations and interconversion with several attitude 
//! representations as generic functions (Supports `f32` & `f64`).
//! 
//! ## Generics
//! 
//! Functions implementing the `QuaternionOps` trait can take both `Quaternion<T>` 
//! and `Vector3<T>` as arguments. In this case, `Vector3<T>` is treated as a Pure Quaternion.
//! 
//! For example:
//! ```
//! use quaternion_core::{Vector3, Quaternion, add};
//! 
//! // --- Vector3 --- //
//! let v1: Vector3<f32> = [1.0, 2.0, 3.0];
//! let v2: Vector3<f32> = [0.1, 0.2, 0.3];
//! println!("{:?}", add(v1, v2));  // <--- It's [1.1, 2.2, 3.3]
//! 
//! // --- Quaternion --- //
//! let q1: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
//! let q2: Quaternion<f64> = (0.1, [0.2, 0.3, 0.4]);
//! println!("{:?}", add(q1, q2));  // <--- It's (1.1, [2.2, 3.3, 4.4])
//! ```
//! 
//! ## Versor
//! 
//! Versor refers to a Quaternion representing a rotation, the norm of which is 1.
//! 
//! The documentation for this crate basically writes Versor instead of Unit Quaternion, 
//! but the difference in usage is not clear.
//! Please think Versor = Unit Quaternion.

#![no_std]
#[cfg(feature = "std")]
extern crate std;

use num_traits::{Float, FloatConst};

mod euler;
mod generics;
pub use generics::QuaternionOps;

/// Vector3 (Pure Quaternion)
/// 
/// `[i, j, k]`
pub type Vector3<T> = [T; 3];

/// Quaternion
/// 
/// `(1, [i, j, k])`
pub type Quaternion<T> = (T, Vector3<T>);

/// Direction Cosine Matrix
/// 
/// `mij`: row `i`, column `j`
/// 
/// `
/// [
///     [m11, m12, m13],
///     [m21, m22, m23],
///     [m31, m32, m33]
/// ]
/// `
pub type DCM<T> = [Vector3<T>; 3];

/// Specify the rotation type of Euler angles.
/// 
/// Considering a fixed `Reference frame` and a rotating `Body frame`, 
/// `Intrinsic rotation` and `Extrinsic rotation` represent the following rotations:
/// 
/// * `Intrinsic`: Rotate around the axes of the `Body frame`
/// * `Extrinsic`: Rotate around the axes of the `Reference frame`
#[derive(Debug, Clone, Copy)]
pub enum RotationType {
    Intrinsic,
    Extrinsic,
}

/// Represents 12 different rotations.
/// 
/// Each variant reads from left to right.
/// For example, `RotationSequence::XYZ` represents rotation around the X axis first, 
/// then the Y axis, and finally the Z axis in that order (X ---> Y ---> Z).
#[derive(Debug, Clone, Copy)]
pub enum RotationSequence {
    // Proper (z-x-z, x-y-x, y-z-y, z-y-z, x-z-x, y-x-y)
    ZXZ,
    XYX,
    YZY,
    ZYZ,
    XZX,
    YXY,
    // Tait–Bryan (x-y-z, y-z-x, z-x-y, x-z-y, z-y-x, y-x-z)
    XYZ,
    YZX,
    ZXY,
    XZY,
    ZYX,
    YXZ,
}

/// Generate Versor by specifying rotation `angle`\[rad\] and `axis` vector.
/// 
/// The `axis` vector does not have to be a unit vector.
/// 
/// If you enter a zero vector, it returns an identity quaternion.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{from_axis_angle, point_rotation, sub};
/// # let PI = std::f64::consts::PI;
/// // Generates a quaternion representing the
/// // rotation of π/2[rad] around the y-axis.
/// let q = from_axis_angle([0.0, 1.0, 0.0], PI/2.0);
/// 
/// // Rotate the point.
/// let r = point_rotation(q, [2.0, 2.0, 0.0]);
/// 
/// // Check if the calculation is correct.
/// let diff = sub([0.0, 2.0, -2.0], r);
/// for val in diff {
///     assert!( val.abs() < 1e-12 );
/// }
/// ```
#[inline]
pub fn from_axis_angle<T>(axis: Vector3<T>, angle: T) -> Quaternion<T>
where T: Float + FloatConst {
    let theta = angle % ( T::PI() + T::PI() );  // limit to (-2π, 2π)
    let f = ( theta * cast(0.5) ).sin_cos();
    let coef = f.0 / norm(axis);
    if coef.is_infinite() {
        IDENTITY()
    } else {
        ( f.1, scale(coef, axis) )
    }
}

/// Calculate the rotation `axis` (unit vector) and the rotation `angle`\[rad\] 
/// around the `axis` from the Versor.
/// 
/// If identity quaternion is entered, `angle` returns zero and 
/// the `axis` returns a zero vector.
/// 
/// Range of `angle`: `(-π, π]`
#[inline]
pub fn to_axis_angle<T>(q: Quaternion<T>) -> (Vector3<T>, T)
where T: Float {
    let norm_q_v = norm(q.1);
    let coef = norm_q_v.recip();
    if coef.is_infinite() {
        ( ZERO_VECTOR(), T::zero() )
    } else {
        // 少しの誤差は見逃す．
        let tmp = norm_q_v.min( T::one() ).asin();
        ( scale(coef, q.1), (tmp + tmp).copysign(q.0) ) // theta = 2*tmp
    }
}

/// Convert a DCM to a Versor representing 
/// the `q v q*` rotation (Point Rotation - Frame Fixed).
/// 
/// When convert to a DCM representing `q* v q` rotation
/// (Frame Rotation - Point Fixed) to a Versor, do the following:
/// 
/// ```
/// # use quaternion_core::{from_dcm, to_dcm, conj};
/// # let dcm = to_dcm((1.0, [0.0; 3]));
/// let q = conj( from_dcm(dcm) );
/// ```
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{
/// #     from_axis_angle, dot, conj, negate, to_dcm, from_dcm,
/// #     matrix_product, point_rotation, frame_rotation
/// # };
/// # let PI = std::f64::consts::PI;
/// // Make these as you like.
/// let v = [1.0, 0.5, -8.0];
/// let q = from_axis_angle([0.2, 1.0, -2.0], PI/4.0);
/// 
/// // --- Point rotation --- //
/// {
///     let m = to_dcm(q);
///     let q_check = from_dcm(m);
///     
///     assert!( (q.0    - q_check.0).abs() < 1e-12 );
///     assert!( (q.1[0] - q_check.1[0]).abs() < 1e-12 );
///     assert!( (q.1[1] - q_check.1[1]).abs() < 1e-12 );
///     assert!( (q.1[2] - q_check.1[2]).abs() < 1e-12 );
/// }
/// 
/// // --- Frame rotation --- //
/// {
///     let m = to_dcm( conj(q) );
///     let q_check = conj( from_dcm(m) );
///     
///     assert!( (q.0    - q_check.0).abs() < 1e-12 );
///     assert!( (q.1[0] - q_check.1[0]).abs() < 1e-12 );
///     assert!( (q.1[1] - q_check.1[1]).abs() < 1e-12 );
///     assert!( (q.1[2] - q_check.1[2]).abs() < 1e-12 );
/// }
/// ```
#[inline]
pub fn from_dcm<T>(m: DCM<T>) -> Quaternion<T>
where T: Float {
    // ゼロ除算を避けるために，4通りの式で求めたうちの最大値を係数として使う．
    let (index, max_num) = max4([
        m[0][0] + m[1][1] + m[2][2],
        m[0][0] - m[1][1] - m[2][2],
       -m[0][0] + m[1][1] - m[2][2],
       -m[0][0] - m[1][1] + m[2][2],
    ]);

    let half: T = cast(0.5);
    let tmp = ( max_num + T::one() ).sqrt();
    let coef = half / tmp;

    let (q0, [q1, q2, q3]): Quaternion<T>;
    match index {
        0 => {
            q0 = half * tmp;
            q1 = (m[2][1] - m[1][2]) * coef;
            q2 = (m[0][2] - m[2][0]) * coef;
            q3 = (m[1][0] - m[0][1]) * coef;
        },
        1 => {
            q0 = (m[2][1] - m[1][2]) * coef;
            q1 = half * tmp;
            q2 = (m[0][1] + m[1][0]) * coef;
            q3 = (m[0][2] + m[2][0]) * coef;
        },
        2 => {
            q0 = (m[0][2] - m[2][0]) * coef;
            q1 = (m[0][1] + m[1][0]) * coef;
            q2 = half * tmp;
            q3 = (m[1][2] + m[2][1]) * coef;
        },
        3 => {
            q0 = (m[1][0] - m[0][1]) * coef;
            q1 = (m[0][2] + m[2][0]) * coef;
            q2 = (m[1][2] + m[2][1]) * coef;
            q3 = half * tmp;
        },
        _ => unreachable!(),
    };

    (q0, [q1, q2, q3])
}

/// Convert a Versor to a DCM representing 
/// the `q v q*` rotation (Point Rotation - Frame Fixed).
/// 
/// When convert to a DCM representing the 
/// `q* v q` rotation (Frame Rotation - Point Fixed), do the following:
/// 
/// ```
/// # use quaternion_core::{to_dcm, conj};
/// # let q = (1.0, [0.0; 3]);
/// let dcm = to_dcm( conj(q) );
/// ```
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{
/// #     from_axis_angle, to_dcm, conj, 
/// #     matrix_product, point_rotation, frame_rotation
/// # };
/// # let PI = std::f64::consts::PI;
/// // Make these as you like.
/// let v = [1.0, 0.5, -8.0];
/// let q = from_axis_angle([0.2, 1.0, -2.0], PI/4.0);
/// 
/// // --- Point rotation --- //
/// {
///     let m = to_dcm(q);
/// 
///     let rm = matrix_product(m, v);
///     let rq = point_rotation(q, v);
///     assert!( (rm[0] - rq[0]).abs() < 1e-12 );
///     assert!( (rm[1] - rq[1]).abs() < 1e-12 );
///     assert!( (rm[2] - rq[2]).abs() < 1e-12 );
/// }
/// 
/// // --- Frame rotation --- //
/// {
///     let m = to_dcm( conj(q) );
/// 
///     let rm = matrix_product(m, v);
///     let rq = frame_rotation(q, v);
///     assert!( (rm[0] - rq[0]).abs() < 1e-12 );
///     assert!( (rm[1] - rq[1]).abs() < 1e-12 );
///     assert!( (rm[2] - rq[2]).abs() < 1e-12 );
/// }
/// ```
#[inline]
pub fn to_dcm<T>(q: Quaternion<T>) -> DCM<T>
where T: Float {
    let neg_one = -T::one();
    let two = cast(2.0);

    // Compute these value only once.
    let (q0_q0, [q0_q1, q0_q2, q0_q3]) = scale(q.0, q);
    let q1_q2 = q.1[0] * q.1[1];
    let q1_q3 = q.1[0] * q.1[2];
    let q2_q3 = q.1[1] * q.1[2];

    let m11 = mul_add(mul_add(q.1[0], q.1[0], q0_q0), two, neg_one);
    let m12 = (q1_q2 - q0_q3) * two;
    let m13 = (q1_q3 + q0_q2) * two;
    let m21 = (q1_q2 + q0_q3) * two;
    let m22 = mul_add(mul_add(q.1[1], q.1[1], q0_q0), two, neg_one);
    let m23 = (q2_q3 - q0_q1) * two;
    let m31 = (q1_q3 - q0_q2) * two;
    let m32 = (q2_q3 + q0_q1) * two;
    let m33 = mul_add(mul_add(q.1[2], q.1[2], q0_q0), two, neg_one);

    [
        [m11, m12, m13],
        [m21, m22, m23],
        [m31, m32, m33],
    ]
}

/// Convert Euler angles to Versor.
/// 
/// The type of rotation (Intrinsic or Extrinsic) is specified by `RotationType` enum, 
/// and the rotation sequence (XZX, XYZ, ...) is specified by `RotationSequence` enum.
/// 
/// Each element of `angles` should be specified in the range: `[-2π, 2π]`.
/// 
/// Sequences: `angles[0]` ---> `angles[1]` ---> `angles[2]`
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{from_axis_angle, mul, from_euler_angles, point_rotation};
/// # let PI = std::f64::consts::PI;
/// use quaternion_core::{RotationType::*, RotationSequence::XYZ};
/// 
/// let angles = [PI/6.0, 1.6*PI, -PI/4.0];
/// let v = [1.0, 0.5, -0.4];
/// 
/// // Quaternions representing rotation around each axis.
/// let x = from_axis_angle([1.0, 0.0, 0.0], angles[0]);
/// let y = from_axis_angle([0.0, 1.0, 0.0], angles[1]);
/// let z = from_axis_angle([0.0, 0.0, 1.0], angles[2]);
/// 
/// // ---- Intrinsic (X-Y-Z) ---- //
/// // These represent the same rotation.
/// let q_in = mul( mul(x, y), z );
/// let e2q_in = from_euler_angles(Intrinsic, XYZ, angles);
/// // Confirmation
/// let a_in = point_rotation(q_in, v);
/// let b_in = point_rotation(e2q_in, v);
/// assert!( (a_in[0] - b_in[0]).abs() < 1e-12 );
/// assert!( (a_in[1] - b_in[1]).abs() < 1e-12 );
/// assert!( (a_in[2] - b_in[2]).abs() < 1e-12 );
/// 
/// // ---- Extrinsic (X-Y-Z) ---- //
/// // These represent the same rotation.
/// let q_ex = mul( mul(z, y), x );
/// let e2q_ex = from_euler_angles(Extrinsic, XYZ, angles);
/// // Confirmation
/// let a_ex = point_rotation(q_ex, v);
/// let b_ex = point_rotation(e2q_ex, v);
/// assert!( (a_ex[0] - b_ex[0]).abs() < 1e-12 );
/// assert!( (a_ex[1] - b_ex[1]).abs() < 1e-12 );
/// assert!( (a_ex[2] - b_ex[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn from_euler_angles<T>(rt: RotationType, rs: RotationSequence, angles: Vector3<T>) -> Quaternion<T>
where T: Float + FloatConst {
    debug_assert!( angles[0].abs() <= T::PI() + T::PI(), "angles[0] is out of range!");
    debug_assert!( angles[1].abs() <= T::PI() + T::PI(), "angles[1] is out of range!");
    debug_assert!( angles[2].abs() <= T::PI() + T::PI(), "angles[2] is out of range!");

    match rt {
        RotationType::Intrinsic => euler::from_intrinsic_euler_angles(rs, angles),
        RotationType::Extrinsic => euler::from_extrinsic_euler_angles(rs, angles),
    }
}

/// Convert Versor to Euler angles.
/// 
/// The type of rotation (Intrinsic or Extrinsic) is specified by `RotationType` enum, 
/// and the rotation sequence (XZX, XYZ, ...) is specified by `RotationSequence` enum.
/// 
/// ```
/// # use quaternion_core::{RotationType::Intrinsic, RotationSequence::XYZ, to_euler_angles};
/// # let q = (1.0, [0.0; 3]);
/// let angles = to_euler_angles(Intrinsic, XYZ, q);
/// ```
/// 
/// Sequences: `angles[0]` ---> `angles[1]` ---> `angles[2]`
/// 
/// # Singularity
/// 
/// ## RotationType::Intrinsic
/// 
/// For Proper Euler angles (ZXZ, XYX, YZY, ZYZ, XZX, YXY), the singularity is reached 
/// when the sine of the second rotation angle is 0 (angle = 0, ±π, ...), and for 
/// Tait-Bryan angles (XYZ, YZX, ZXY, XZY, ZYX, YXZ), the singularity is reached when 
/// the cosine of the second rotation angle is 0 (angle = ±π/2).
/// 
/// At the singularity, the third rotation angle is set to 0\[rad\].
/// 
/// ## RotationType::Extrinsic
/// 
/// As in the case of Intrinsic rotation, for Proper Euler angles, the singularity occurs 
/// when the sine of the second rotation angle is 0 (angle = 0, ±π, ...), and for 
/// Tait-Bryan angles, the singularity occurs when the cosine of the second rotation angle 
/// is 0 (angle = ±π/2).
/// 
/// At the singularity, the first rotation angle is set to 0\[rad\].
/// 
/// # Examples
/// 
/// Depending on the rotation angle of each axis, it may not be possible to recover the 
/// same rotation angle as the original. However, they represent the same rotation in 3D space.
/// 
/// ```
/// # use quaternion_core::{from_euler_angles, to_euler_angles, point_rotation};
/// # let PI = std::f64::consts::PI;
/// use quaternion_core::{RotationType::*, RotationSequence::XYZ};
/// 
/// let angles = [PI/6.0, PI/4.0, PI/3.0];
/// 
/// // ---- Intrinsic (X-Y-Z) ---- //
/// let q_in = from_euler_angles(Intrinsic, XYZ, angles);
/// let e_in = to_euler_angles(Intrinsic, XYZ, q_in);
/// assert!( (angles[0] - e_in[0]).abs() < 1e-12 );
/// assert!( (angles[1] - e_in[1]).abs() < 1e-12 );
/// assert!( (angles[2] - e_in[2]).abs() < 1e-12 );
/// 
/// // ---- Extrinsic (X-Y-Z) ---- //
/// let q_ex = from_euler_angles(Extrinsic, XYZ, angles);
/// let e_ex = to_euler_angles(Extrinsic, XYZ, q_ex);
/// assert!( (angles[0] - e_ex[0]).abs() < 1e-12 );
/// assert!( (angles[1] - e_ex[1]).abs() < 1e-12 );
/// assert!( (angles[2] - e_ex[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn to_euler_angles<T>(rt: RotationType, rs: RotationSequence, q: Quaternion<T>) -> Vector3<T>
where T: Float + FloatConst {
    match rt {
        RotationType::Intrinsic => euler::to_intrinsic_euler_angles(rs, q),
        RotationType::Extrinsic => euler::to_extrinsic_euler_angles(rs, q),
    }
}

/// Convert Rotation vector to Versor.
/// 
/// The Rotation vector itself represents the axis of rotation, 
/// and the norm represents the angle of rotation around the axis.
/// 
/// Range of the norm of the rotation vector: `[0, 2π]`
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{from_rotation_vector, scale, point_rotation};
/// # let PI = std::f64::consts::PI;
/// let angle = PI / 2.0;
/// let axis = [1.0, 0.0, 0.0];
/// 
/// // This represents a rotation of π/2 around the x-axis.
/// let rot_vec = scale(angle, axis);  // Rotation vector
/// 
/// // Rotation vector ---> Quaternion
/// let q = from_rotation_vector(rot_vec);
/// 
/// let r = point_rotation(q, [1.0, 1.0, 0.0]);
/// 
/// assert!( (r[0] - 1.0).abs() < 1e-12 );
/// assert!( (r[1] - 0.0).abs() < 1e-12 );
/// assert!( (r[2] - 1.0).abs() < 1e-12 );
/// ```
#[inline]
pub fn from_rotation_vector<T>(v: Vector3<T>) -> Quaternion<T>
where T: Float {
    let theta = norm(v);
    let f = ( theta * cast(0.5) ).sin_cos();
    let coef = f.0 / theta;
    if coef.is_infinite() {
        IDENTITY()
    } else {
        ( f.1, scale(coef, v) )
    }
}

/// Convert Versor to Rotation vector.
/// 
/// The Rotation vector itself represents the axis of rotation, 
/// and the norm represents the angle of rotation around the axis.
/// 
/// Range of the norm of the rotation vector: `[0, 2π]`
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{from_axis_angle, to_rotation_vector, scale};
/// # let PI = std::f64::consts::PI;
/// let angle = PI / 2.0;
/// let axis = [1.0, 0.0, 0.0];
/// 
/// // These represent the same rotation.
/// let rv = scale(angle, axis);  // Rotation vector
/// let q = from_axis_angle(axis, angle);  // Quaternion
/// 
/// // Quaternion ---> Rotation vector
/// let q2rv = to_rotation_vector(q);
/// 
/// assert!( (rv[0] - q2rv[0]).abs() < 1e-12 );
/// assert!( (rv[1] - q2rv[1]).abs() < 1e-12 );
/// assert!( (rv[2] - q2rv[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn to_rotation_vector<T>(q: Quaternion<T>) -> Vector3<T>
where T: Float {
    let tmp = acos_safe(q.0);
    let coef = (tmp + tmp) / norm(q.1);  // 2*tmp
    if coef.is_infinite() {
        ZERO_VECTOR()
    } else {
        scale(coef, q.1)
    }
}

/// Product of DCM and Vector3
#[inline]
pub fn matrix_product<T>(m: DCM<T>, v: Vector3<T>) -> Vector3<T>
where T: Float {
    [
        dot(m[0], v),
        dot(m[1], v),
        dot(m[2], v),
    ]
}

/// Calculate the sum of each element of Quaternion or Vector3.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, sum};
/// // --- Vector3 --- //
/// let v: Vector3<f64> = [1.0, 2.0, 3.0];
/// 
/// assert!( (6.0 - sum(v)).abs() < 1e-12 );
/// 
/// // --- Quaternion --- //
/// let q: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// 
/// assert!( (10.0 - sum(q)).abs() < 1e-12 );
/// ```
#[inline]
pub fn sum<T, U>(a: U) -> T
where T: Float, U: QuaternionOps<T> {
    a.sum()
}

/// `a + b`
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, add};
/// // --- Vector3 --- //
/// let v1: Vector3<f64> = [1.0, 2.0, 3.0];
/// let v2: Vector3<f64> = [0.1, 0.2, 0.3];
/// let v_result = add(v1, v2);
/// 
/// assert!( (1.1 - v_result[0]).abs() < 1e-12 );
/// assert!( (2.2 - v_result[1]).abs() < 1e-12 );
/// assert!( (3.3 - v_result[2]).abs() < 1e-12 );
/// 
/// // --- Quaternion --- //
/// let q1: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// let q2: Quaternion<f64> = (0.1, [0.2, 0.3, 0.4]);
/// let q_result = add(q1, q2);
/// 
/// assert!( (1.1 - q_result.0).abs() < 1e-12 );
/// assert!( (2.2 - q_result.1[0]).abs() < 1e-12 );
/// assert!( (3.3 - q_result.1[1]).abs() < 1e-12 );
/// assert!( (4.4 - q_result.1[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn add<T, U>(a: U, b: U) -> U
where T: Float, U: QuaternionOps<T> {
    a.add(b)
}

/// `a - b`
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, sub};
/// // --- Vector3 --- //
/// let v1: Vector3<f64> = [1.0, 2.0, 3.0];
/// let v2: Vector3<f64> = [0.1, 0.2, 0.3];
/// let v_result = sub(v1, v2);
/// 
/// assert!( (0.9 - v_result[0]).abs() < 1e-12 );
/// assert!( (1.8 - v_result[1]).abs() < 1e-12 );
/// assert!( (2.7 - v_result[2]).abs() < 1e-12 );
/// 
/// // --- Quaternion --- //
/// let q1: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// let q2: Quaternion<f64> = (0.1, [0.2, 0.3, 0.4]);
/// let q_result = sub(q1, q2);
/// 
/// assert!( (0.9 - q_result.0).abs() < 1e-12 );
/// assert!( (1.8 - q_result.1[0]).abs() < 1e-12 );
/// assert!( (2.7 - q_result.1[1]).abs() < 1e-12 );
/// assert!( (3.6 - q_result.1[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn sub<T, U>(a: U, b: U) -> U
where T: Float, U: QuaternionOps<T> {
    a.sub(b)
}

/// `s * a`
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, scale};
/// // --- Vector3 --- //
/// let v: Vector3<f64> = [1.0, 2.0, 3.0];
/// let v_result = scale(2.0, v);
/// 
/// assert!( (2.0 - v_result[0]).abs() < 1e-12 );
/// assert!( (4.0 - v_result[1]).abs() < 1e-12 );
/// assert!( (6.0 - v_result[2]).abs() < 1e-12 );
/// 
/// // --- Quaternion --- //
/// let q: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// let q_result = scale(2.0, q);
/// 
/// assert!( (2.0 - q_result.0).abs() < 1e-12 );
/// assert!( (4.0 - q_result.1[0]).abs() < 1e-12 );
/// assert!( (6.0 - q_result.1[1]).abs() < 1e-12 );
/// assert!( (8.0 - q_result.1[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn scale<T, U>(s: T, a: U) -> U
where T: Float, U: QuaternionOps<T> {
    a.scale(s)
}

/// `s * a + b`
/// 
/// If the `fma` feature is enabled, the FMA calculation is performed using 
/// the `mul_add` method (`s.mul_add(a, b)`). 
/// If not enabled, it's computed by unfused multiply-add (`s * a + b`).
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, scale_add};
/// // --- Vector3 --- //
/// let v1: Vector3<f64> = [1.0, 2.0, 3.0];
/// let v2: Vector3<f64> = [0.1, 0.2, 0.3];
/// let v_result = scale_add(2.0, v1, v2);
/// 
/// assert!( (2.1 - v_result[0]).abs() < 1e-12 );
/// assert!( (4.2 - v_result[1]).abs() < 1e-12 );
/// assert!( (6.3 - v_result[2]).abs() < 1e-12 );
/// 
/// // --- Quaternion --- //
/// let q1: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// let q2: Quaternion<f64> = (0.1, [0.2, 0.3, 0.4]);
/// let q_result = scale_add(2.0, q1, q2);
/// 
/// assert!( (2.1 - q_result.0).abs() < 1e-12 );
/// assert!( (4.2 - q_result.1[0]).abs() < 1e-12 );
/// assert!( (6.3 - q_result.1[1]).abs() < 1e-12 );
/// assert!( (8.4 - q_result.1[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn scale_add<T, U>(s: T, a: U, b: U) -> U
where T: Float, U: QuaternionOps<T> {
    a.scale_add(s, b)
}

/// `a ∘ b`
/// 
/// Hadamard product of Vector3 or Quaternion.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, hadamard};
/// // --- Vector3 --- //
/// let v1: Vector3<f64> = [1.0, 2.0, 3.0];
/// let v2: Vector3<f64> = [0.1, 0.2, 0.3];
/// let v_result = hadamard(v1, v2);
/// 
/// assert!( (0.1 - v_result[0]).abs() < 1e-12 );
/// assert!( (0.4 - v_result[1]).abs() < 1e-12 );
/// assert!( (0.9 - v_result[2]).abs() < 1e-12 );
/// 
/// // --- Quaternion --- //
/// let q1: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// let q2: Quaternion<f64> = (0.1, [0.2, 0.3, 0.4]);
/// let q_result = hadamard(q1, q2);
/// 
/// assert!( (0.1 - q_result.0).abs() < 1e-12 );
/// assert!( (0.4 - q_result.1[0]).abs() < 1e-12 );
/// assert!( (0.9 - q_result.1[1]).abs() < 1e-12 );
/// assert!( (1.6 - q_result.1[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn hadamard<T, U>(a: U, b: U) -> U
where T: Float, U: QuaternionOps<T> {
    a.hadamard(b)
}

/// `a ∘ b + c`
/// 
/// Hadamard product and addiction of Quaternion or Vector3.
/// 
/// If the `fma` feature is enabled, the FMA calculation is performed using 
/// the `mul_add` method (`s.mul_add(a, b)`). 
/// If not enabled, it's computed by unfused multiply-add (`s * a + b`).
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, hadamard_add};
/// // --- Vector3 --- //
/// let v1: Vector3<f64> = [1.0, 2.0, 3.0];
/// let v2: Vector3<f64> = [0.1, 0.2, 0.3];
/// let v3: Vector3<f64> = [0.4, 0.5, 0.6];
/// let v_result = hadamard_add(v1, v2, v3);
/// 
/// assert!( (0.5 - v_result[0]).abs() < 1e-12 );
/// assert!( (0.9 - v_result[1]).abs() < 1e-12 );
/// assert!( (1.5 - v_result[2]).abs() < 1e-12 );
/// 
/// // --- Quaternion --- //
/// let q1: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// let q2: Quaternion<f64> = (0.1, [0.2, 0.3, 0.4]);
/// let q3: Quaternion<f64> = (0.5, [0.6, 0.7, 0.8]);
/// let q_result = hadamard_add(q1, q2, q3);
/// 
/// assert!( (0.6 - q_result.0).abs() < 1e-12 );
/// assert!( (1.0 - q_result.1[0]).abs() < 1e-12 );
/// assert!( (1.6 - q_result.1[1]).abs() < 1e-12 );
/// assert!( (2.4 - q_result.1[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn hadamard_add<T, U>(a: U, b: U, c: U) -> U
where T: Float, U: QuaternionOps<T> {
    a.hadamard_add(b, c)
}

/// `a · b`
/// 
/// Dot product of Vector3 or Quaternion.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, dot};
/// // --- Vector3 --- //
/// let v1: Vector3<f64> = [1.0, 2.0, 3.0];
/// let v2: Vector3<f64> = [0.1, 0.2, 0.3];
/// 
/// assert!( (1.4 - dot(v1, v2)).abs() < 1e-12 );
/// 
/// // --- Quaternion --- //
/// let q1: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// let q2: Quaternion<f64> = (0.1, [0.2, 0.3, 0.4]);
/// 
/// assert!( (3.0 - dot(q1, q2)).abs() < 1e-12 );
/// ```
#[inline]
pub fn dot<T, U>(a: U, b: U) -> T 
where T: Float, U: QuaternionOps<T> {
    sum( hadamard(a, b) )
}

/// Cross product (vector product): `a × b`
/// 
/// The product order is `a × b (!= b × a)`
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, scale, cross};
/// let v1: Vector3<f64> = [0.5, -1.0, 0.8];
/// let v2: Vector3<f64> = scale(2.0, v1);
/// let v_result = cross(v1, v2);
/// 
/// // The cross product of parallel vectors is a zero vector.
/// assert!( v_result[0].abs() < 1e-12 );
/// assert!( v_result[1].abs() < 1e-12 );
/// assert!( v_result[2].abs() < 1e-12 );
/// ```
#[inline]
pub fn cross<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T>
where T: Float {
    [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    ]
}

/// Calculate L2 norm of Vector3 or Quaternion.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, norm};
/// // --- Vector3 --- //
/// let v: Vector3<f64> = [1.0, 2.0, 3.0];
/// 
/// assert!( (14.0_f64.sqrt() - norm(v)).abs() < 1e-12 );
/// 
/// // --- Quaternion --- //
/// let q: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// 
/// assert!( (30.0_f64.sqrt() - norm(q)).abs() < 1e-12 );
/// ```
#[inline]
pub fn norm<T, U>(a: U) -> T 
where T: Float, U: QuaternionOps<T> {
    dot(a, a).sqrt()
}

/// Normalization of Vector3 or Quaternion.
/// 
/// If you enter a zero vector, it returns a zero vector.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, norm, normalize};
/// // --- Vector3 --- //
/// // This norm is not 1.
/// let v: Vector3<f64> = [1.0, 2.0, 3.0];
/// assert!( (1.0 - norm(v)).abs() > 1e-12 );
/// 
/// // Now that normalized, this norm is 1!
/// let v_n = normalize(v);
/// assert!( (1.0 - norm(v_n)).abs() < 1e-12 );
/// 
/// // --- Quaternion --- //
/// // This norm is not 1.
/// let q: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// assert!( (1.0 - norm(q)).abs() > 1e-12 );
/// 
/// // Now that normalized, this norm is 1!
/// let q_n = normalize(q);
/// assert!( (1.0 - norm(q_n)).abs() < 1e-12 );
/// ```
#[inline]
pub fn normalize<T, U>(a: U) -> U
where T: Float, U: QuaternionOps<T> {
    a.normalize()
}

/// Invert the sign of a Vector3 or Quaternion.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, negate};
/// // --- Vector3 --- //
/// let v: Vector3<f64> = [1.0, 2.0, 3.0];
/// let v_n = negate(v);
/// 
/// assert_eq!(-v[0], v_n[0]);
/// assert_eq!(-v[1], v_n[1]);
/// assert_eq!(-v[2], v_n[2]);
/// 
/// // --- Quaternion --- //
/// let q: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// let q_n = negate(q);
/// 
/// assert_eq!(-q.0,    q_n.0);
/// assert_eq!(-q.1[0], q_n.1[0]);
/// assert_eq!(-q.1[1], q_n.1[1]);
/// assert_eq!(-q.1[2], q_n.1[2]);
/// ```
#[inline]
pub fn negate<T, U>(a: U) -> U
where T: Float, U: QuaternionOps<T> {
    a.negate()
}

/// Hamilton product (Product of Quaternion or Pure Quaternion)
/// 
/// The product order is `ab (!= ba)`
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, inv, mul};
/// // ---- Pure Quaternion (Vector3) ---- //
/// let v: Vector3<f64> = [1.0, 2.0, 3.0];
/// 
/// // Identity quaternion
/// let id = mul( v, inv(v) );  // = mul( inv(v), v );
/// 
/// assert!( (1.0 - id.0).abs() < 1e-12 );
/// assert!( id.1[0].abs() < 1e-12 );
/// assert!( id.1[1].abs() < 1e-12 );
/// assert!( id.1[2].abs() < 1e-12 );
/// 
/// // ---- Quaternion ---- //
/// let q: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// 
/// // Identity quaternion
/// let id = mul( q, inv(q) );  // = mul( inv(q), q );
/// 
/// assert!( (1.0 - id.0).abs() < 1e-12 );
/// assert!( id.1[0].abs() < 1e-12 );
/// assert!( id.1[1].abs() < 1e-12 );
/// assert!( id.1[2].abs() < 1e-12 );
/// ```
#[inline]
pub fn mul<T, U>(a: U, b: U) -> Quaternion<T>
where T: Float, U: QuaternionOps<T> {
    a.mul(b)
}

/// Calculate the conjugate of Quaternion.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Quaternion, conj};
/// let q: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// let q_conj = conj(q);
/// 
/// assert_eq!( q.0,    q_conj.0);
/// assert_eq!(-q.1[0], q_conj.1[0]);
/// assert_eq!(-q.1[1], q_conj.1[1]);
/// assert_eq!(-q.1[2], q_conj.1[2]);
/// ```
#[inline]
pub fn conj<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( q.0, negate(q.1) )
}

/// Calculate the inverse of Quaternion or Pure Quaternion (Vector3).
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, inv, mul};
/// // ---- Pure Quaternion (Vector3) ---- //
/// let v: Vector3<f64> = [1.0, 2.0, 3.0];
/// 
/// // Identity quaternion
/// let id = mul( v, inv(v) );  // = mul( inv(v), v );
/// 
/// assert!( (id.0 - 1.0).abs() < 1e-12 );
/// assert!( id.1[0].abs() < 1e-12 );
/// assert!( id.1[1].abs() < 1e-12 );
/// assert!( id.1[2].abs() < 1e-12 );
/// 
/// // ---- Quaternion ---- //
/// let q: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// 
/// // Identity quaternion
/// let id = mul( q, inv(q) );  // = mul( inv(q), q );
/// 
/// assert!( (id.0 - 1.0).abs() < 1e-12 );
/// assert!( id.1[0].abs() < 1e-12 );
/// assert!( id.1[1].abs() < 1e-12 );
/// assert!( id.1[2].abs() < 1e-12 );
/// ```
#[inline]
pub fn inv<T, U>(a: U) -> U
where T: Float, U: QuaternionOps<T> {
    a.inv()
}

// acosは[-π/2, π/2]の範囲でしか値を返さないので、qのとり方によってはlnで完全に復元できない。
// q == ln(exp(q)) が成り立つのはcos(norm(q.1))が[-π/2, π/2]の範囲内にある場合のみ。
/// Exponential function of Quaternion or Pure Quaternion (Vector3).
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, scale, exp, ln};
/// // ---- Pure Quaternion (Vector3) ---- //
/// let v: Vector3<f64> = [0.1, 0.2, 0.3];
/// let v_r = ln( exp(v) );
/// 
/// assert!( v_r.0.abs() < 1e-12 );
/// assert!( (v[0] - v_r.1[0]).abs() < 1e-12 );
/// assert!( (v[1] - v_r.1[1]).abs() < 1e-12 );
/// assert!( (v[2] - v_r.1[2]).abs() < 1e-12 );
/// 
/// // ---- Quaternion ---- //
/// let q: Quaternion<f64> = (0.1, [0.2, 0.3, 0.4]);
/// let q_r = ln( exp(q) );
/// 
/// assert!( (q.0    - q_r.0).abs() < 1e-12 );
/// assert!( (q.1[0] - q_r.1[0]).abs() < 1e-12 );
/// assert!( (q.1[1] - q_r.1[1]).abs() < 1e-12 );
/// assert!( (q.1[2] - q_r.1[2]).abs() < 1e-12 );
/// 
/// // The relationship between exp(q) and exp(q.1)
/// let exp_q = exp(q);
/// let exp_check = scale( q.0.exp(), exp(q.1) );
/// assert!( (exp_q.0    - exp_check.0).abs() < 1e-12 );
/// assert!( (exp_q.1[0] - exp_check.1[0]).abs() < 1e-12 );
/// assert!( (exp_q.1[1] - exp_check.1[1]).abs() < 1e-12 );
/// assert!( (exp_q.1[2] - exp_check.1[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn exp<T, U>(a: U) -> Quaternion<T>
where T: Float, U: QuaternionOps<T> {
    a.exp()
}

// acosは[-π/2, π/2]の範囲でしか値を返さないので、qのとり方によってはlnで完全に復元できない。
// q == ln(exp(q)) が成り立つのはcos(norm(q.1))が[-π/2, π/2]の範囲内にある場合のみ。
/// Natural logarithm of Quaternion.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, exp, ln};
/// let q: Quaternion<f64> = (0.1, [0.2, 0.3, 0.4]);
/// let q_r = ln( exp(q) );
/// 
/// assert!( (q.0    - q_r.0).abs() < 1e-12 );
/// assert!( (q.1[0] - q_r.1[0]).abs() < 1e-12 );
/// assert!( (q.1[1] - q_r.1[1]).abs() < 1e-12 );
/// assert!( (q.1[2] - q_r.1[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn ln<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    let tmp = dot(q.1, q.1);
    let norm_q = mul_add(q.0, q.0, tmp).sqrt();
    let coef = (q.0 / norm_q).acos() / tmp.sqrt();
    ( norm_q.ln(), scale(coef, q.1) )
}

/// Natural logarithm of Versor.
/// 
/// If the argument `q` is guaranteed to be a Versor,
/// it is less computationally expensive than the `ln(...)` function.
/// 
/// Only the vector part is returned since the real part is always zero.
#[inline]
pub fn ln_versor<T>(q: Quaternion<T>) -> Vector3<T>
where T: Float {
    scale( acos_safe(q.0) / norm(q.1), q.1)
}

/// Power function of Quaternion.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Quaternion, mul, inv, pow, sqrt};
/// let q: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// 
/// let q_q = mul(q, q);
/// let q_pow_2 = pow(q, 2.0);
/// assert!( (q_q.0    - q_pow_2.0).abs() < 1e-12 );
/// assert!( (q_q.1[0] - q_pow_2.1[0]).abs() < 1e-12 );
/// assert!( (q_q.1[1] - q_pow_2.1[1]).abs() < 1e-12 );
/// assert!( (q_q.1[2] - q_pow_2.1[2]).abs() < 1e-12 );
/// 
/// let q_sqrt = sqrt(q);
/// let q_pow_0p5 = pow(q, 0.5);
/// assert!( (q_sqrt.0    - q_pow_0p5.0).abs() < 1e-12 );
/// assert!( (q_sqrt.1[0] - q_pow_0p5.1[0]).abs() < 1e-12 );
/// assert!( (q_sqrt.1[1] - q_pow_0p5.1[1]).abs() < 1e-12 );
/// assert!( (q_sqrt.1[2] - q_pow_0p5.1[2]).abs() < 1e-12 );
/// 
/// let q_inv = inv(q);
/// let q_pow_m1 = pow(q, -1.0);
/// assert!( (q_inv.0    - q_pow_m1.0).abs() < 1e-12 );
/// assert!( (q_inv.1[0] - q_pow_m1.1[0]).abs() < 1e-12 );
/// assert!( (q_inv.1[1] - q_pow_m1.1[1]).abs() < 1e-12 );
/// assert!( (q_inv.1[2] - q_pow_m1.1[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn pow<T>(q: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float {
    let tmp = dot(q.1, q.1);
    let norm_q = mul_add(q.0, q.0, tmp).sqrt();
    let omega = (q.0 / norm_q).acos();
    let (sin, cos) = (t * omega).sin_cos();
    let coef = norm_q.powf(t);
    ( coef * cos, scale((coef * sin) / tmp.sqrt(), q.1) )
}

/// Power function of Versor.
/// 
/// If the argument `q` is guaranteed to be a Versor, it is less 
/// computationally expensive than the `pow(...)` function. 
#[inline]
pub fn pow_versor<T>(q: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float {
    let (sin, cos) = (t * acos_safe(q.0)).sin_cos();
    ( cos, scale(sin / norm(q.1), q.1) )
}

/// Square root of Quaternion.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Quaternion, mul, pow, sqrt};
/// let q: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// let q_sqrt = sqrt(q);
/// 
/// let result = mul(q_sqrt, q_sqrt);
/// assert!( (q.0    - result.0).abs() < 1e-12 );
/// assert!( (q.1[0] - result.1[0]).abs() < 1e-12 );
/// assert!( (q.1[1] - result.1[1]).abs() < 1e-12 );
/// assert!( (q.1[2] - result.1[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn sqrt<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    let half = cast(0.5);
    let dot_v = dot(q.1, q.1);
    let norm_q = mul_add(q.0, q.0, dot_v).sqrt();
    let coef = (((norm_q - q.0) * half) / dot_v).sqrt();
    ( ((norm_q + q.0) * half).sqrt(), scale(coef, q.1) )
}

/// Rotation of point (Point Rotation - Frame Fixed)
/// 
/// `q v q*  (||q|| = 1)`
/// 
/// Since it is implemented with an optimized formula, 
/// it can be calculated with the amount of operations shown in the table below:
/// 
/// | Operation    | Num |
/// |:------------:|:---:|
/// | Multiply     | 18  |
/// | Add/Subtract | 12  |
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{from_axis_angle, point_rotation, mul, conj};
/// # let PI = std::f64::consts::PI;
/// // Make these as you like.
/// let v = [1.0, 0.5, -8.0];
/// let q = from_axis_angle([0.2, 1.0, -2.0], PI);
/// 
/// let r = point_rotation(q, v);
/// 
/// // This makes a lot of wasted calculations.
/// let r_check = mul( mul(q, (0.0, v)), conj(q) ).1;
/// 
/// assert!( (r[0] - r_check[0]).abs() < 1e-12 );
/// assert!( (r[1] - r_check[1]).abs() < 1e-12 );
/// assert!( (r[2] - r_check[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn point_rotation<T>(q: Quaternion<T>, v: Vector3<T>) -> Vector3<T>
where T: Float {
    let tmp = scale_add(q.0, v, cross(q.1, v));
    scale_add(cast(2.0), cross(q.1, tmp), v)
}

/// Rotation of frame (Frame Rotation - Point Fixed)
/// 
/// `q* v q  (||q|| = 1)`
/// 
/// Since it is implemented with an optimized formula, 
/// it can be calculated with the amount of operations shown in the table below:
/// 
/// | Operation    | Num |
/// |:------------:|:---:|
/// | Multiply     | 18  |
/// | Add/Subtract | 12  |
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{from_axis_angle, point_rotation, mul, conj};
/// # let PI = std::f64::consts::PI;
/// // Make these as you like.
/// let v = [1.0, 0.5, -8.0];
/// let q = from_axis_angle([0.2, 1.0, -2.0], PI);
/// 
/// let r = point_rotation(q, v);
/// 
/// // This makes a lot of wasted calculations.
/// let r_check = mul( mul(conj(q), (0.0, v)), q ).1;
/// 
/// assert!( (r[0] - r_check[0]).abs() < 1e-12 );
/// assert!( (r[1] - r_check[1]).abs() < 1e-12 );
/// assert!( (r[2] - r_check[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn frame_rotation<T>(q: Quaternion<T>, v: Vector3<T>) -> Vector3<T>
where T: Float {
    let tmp = scale_add(q.0, v, cross(v, q.1));
    scale_add(cast(2.0), cross(tmp, q.1), v)
}

/// Calculate a Versor to rotate from vector `a` to `b`.
/// 
/// If you enter a zero vector, it returns an identity quaternion.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, cross, rotate_a_to_b, point_rotation};
/// let a: Vector3<f64> = [1.5, -0.5, 0.2];
/// let b: Vector3<f64> = [0.1, 0.6, 1.0];
/// 
/// let q = rotate_a_to_b(a, b);
/// let b_check = point_rotation(q, a);
/// 
/// let cross = cross(b, b_check);
/// assert!( cross[0].abs() < 1e-12 );
/// assert!( cross[1].abs() < 1e-12 );
/// assert!( cross[2].abs() < 1e-12 );
/// ```
#[inline]
pub fn rotate_a_to_b<T>(a: Vector3<T>, b: Vector3<T>) -> Quaternion<T>
where T: Float {
    let half: T = cast(0.5);

    let t = dot(a, b);
    let s_square = dot(a, a) * dot(b, b);
    let e_half = half * (t / s_square.sqrt());
    let v = ((half - e_half) / (s_square - t * t)).sqrt();

    // vがfiniteならeもfiniteである．
    if v.is_finite() {
        ( (half + e_half).sqrt(), scale(v, cross(a, b)) )
    } else {
        IDENTITY()
    }
}

/// Calculate a Versor to rotate from vector `a` to `b`.
/// 
/// The parameter `t` adjusts the amount of movement from `a` to `b`, 
/// so that When `t = 1`, it moves to position `b` completely.
/// 
/// If you enter a zero vector, it returns an identity quaternion.
/// 
/// If `t = 1` at all times, it is less computationally expensive to 
/// use `rotate_a_to_b(...)` function.
#[inline]
pub fn rotate_a_to_b_param<T>(a: Vector3<T>, b: Vector3<T>, t: T) -> Quaternion<T>
where T: Float {
    let dot_ab = dot(a, b);
    let norm_ab_square = dot(a, a) * dot(b, b);
    let tmp_acos = dot_ab / norm_ab_square.sqrt();
    if tmp_acos.is_infinite() {
        IDENTITY()
    } else {
        let theta = acos_safe(tmp_acos);
        let (sin, cos) = ( t * theta * cast(0.5) ).sin_cos();
        let coef_v = sin / (norm_ab_square - dot_ab * dot_ab).sqrt();
        if coef_v.is_finite() {
            ( cos, scale(coef_v, cross(a, b)) )
        } else {
            IDENTITY()
        }
    }
}

/// Lerp (Linear interpolation)
/// 
/// Generate a Versor that interpolate the shortest path from `a` to `b`.
/// The argument `t (0 <= t <= 1)` is the interpolation parameter.
/// 
/// The arguments `a` and `b` must be Versor.
/// 
/// Normalization is not performed internally because 
/// it increases the computational complexity.
#[inline]
pub fn lerp<T>(a: Quaternion<T>, b: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float {
    debug_assert!(
        t >= T::zero() && t <= T::one(), 
        "Parameter `t` is out of range!"
    );

    // 最短経路で補間する
    if dot(a, b).is_sign_negative() {
        // bの符号を反転
        if cfg!(feature = "fma") {
            scale_add(-t, add(a, b), a)
        } else {
            sub( a, scale(t, add(a, b)) )
        }
    } else {
        scale_add( t, sub(b, a), a)
    }
}

/// Slerp (Spherical linear interpolation)
/// 
/// Generate a Versor that interpolate the shortest path from `a` to `b`.
/// The argument `t (0 <= t <= 1)` is the interpolation parameter.
/// 
/// The arguments `a` and `b` must be Versor.
#[inline]
pub fn slerp<T>(a: Quaternion<T>, mut b: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float {
    debug_assert!(
        t >= T::zero() && t <= T::one(), 
        "Parameter `t` is out of range!"
    );
    
    // 最短経路で補間する
    let mut dot = dot(a, b);
    if dot.is_sign_negative() {
        b = negate(b);
        dot = -dot;
    }
    // If the distance between quaternions is close enough, use lerp.
    if dot > cast(0.9995) {  // Approximation error < 0.017%
        normalize( scale_add(t, sub(b, a), a) )  // lerp
    } else {
        let omega = dot.acos();  // Angle between the two quaternions.
        let tmp = t * omega;
        let s1 = (omega - tmp).sin();
        let s2 = tmp.sin();
        let coef = (T::one() - dot*dot).sqrt().recip();
        let term1 = scale(s1 * coef, a);
        let term2 = scale(s2 * coef, b);
        add(term1, term2)
    }
}

// ============================================================================= //
// Private functions
// ============================================================================= //
/// Identity quaternion
#[inline(always)]
#[allow(non_snake_case)]
fn IDENTITY<T: Float>() -> Quaternion<T> {
    (T::one(), [T::zero(); 3])
}

#[inline(always)]
#[allow(non_snake_case)]
fn ZERO_VECTOR<T: Float>() -> Vector3<T> {
    [T::zero(); 3]
}

/// 定数呼び出し以外に使わないのでエラー処理を省略．
#[inline(always)]
fn cast<T: Float>(x: f64) -> T {
    num_traits::cast::<f64, T>(x).unwrap()
}

/// `fma` featureを有効にした場合は`s.mul_add(a, b)`として展開され，
/// 有効にしなかった場合は単純な積和`s*a + b`に展開してコンパイルされる．
#[inline(always)]
fn mul_add<T: Float>(s: T, a: T, b: T) -> T {
    if cfg!(feature = "fma") {
        s.mul_add(a, b)
    } else {
        s * a + b
    }
}

/// 配列内の最大値とそのインデックスを返す．
#[inline(always)]
fn max4<T: Float>(nums: [T; 4]) -> (usize, T) {
    let mut index = 0;
    let mut max_num = nums[0];
    for (i, num) in nums.iter().enumerate().skip(1) {
        if *num > max_num {
            max_num = *num;
            index = i;
        }
    }
    (index, max_num)
}

/// 定義域外の値をカットして未定義動作を防ぐ．
#[inline(always)]
fn acos_safe<T: Float>(x: T) -> T {
    x.abs().min( T::one() ).copysign(x).acos()
}