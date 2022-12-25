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

use core::{mem, mem::MaybeUninit};
#[cfg(any(feature = "std", feature = "libm"))]
use num_traits::float::{Float, FloatConst};

mod euler;
mod generics;
pub use generics::QuaternionOps;

/// Vector3 (Pure Quaternion)
/// 
/// The type `[q1, q2, q3]` is equivalent to the expression `q1i + q2j + q3k`,
/// where `i`, `j`, `k` are basis of quaternions and satisfy the following equality:
/// 
/// `i^2 = j^2 = k^2 = ijk = -1`
pub type Vector3<T> = [T; 3];

/// Quaternion
/// 
/// The type `(q0, [q1, q2, q3])` is equivalent to the expression `q0 + q1i + q2j + q3k`, 
/// where `1`, `i`, `j`, `k` are basis of quaternions and satisfy the following equality:
/// 
/// `i^2 = j^2 = k^2 = ijk = -1`
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
/// # use quaternion_core::{from_axis_angle, point_rotation};
/// # let PI = std::f64::consts::PI;
/// // Generates a quaternion representing the
/// // rotation of π/2[rad] around the y-axis.
/// let q = from_axis_angle([0.0, 1.0, 0.0], PI/2.0);
/// 
/// // Rotate the point.
/// let r = point_rotation(q, [2.0, 2.0, 0.0]);
/// 
/// assert!( (r[0] - 0.0).abs() < 1e-12 );
/// assert!( (r[1] - 2.0).abs() < 1e-12 );
/// assert!( (r[2] + 2.0).abs() < 1e-12 );
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
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{from_axis_angle, to_axis_angle};
/// # let PI = std::f64::consts::PI;
/// let axis_ori = [0.0, 1.0, 2.0];
/// let angle_ori = PI / 2.0;
/// let q = from_axis_angle(axis_ori, angle_ori);
/// 
/// let (axis, angle) = to_axis_angle(q);
/// 
/// assert!( (axis_ori[0] - axis[0]).abs() < 1e-12 );
/// assert!( (axis_ori[0] - axis[0]).abs() < 1e-12 );
/// assert!( (axis_ori[0] - axis[0]).abs() < 1e-12 );
/// assert!( (angle_ori - angle).abs() < 1e-12 );
/// ```
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
    let m22_p_m33 = m[1][1] + m[2][2];
    let m22_m_m33 = m[1][1] - m[2][2];
    let (index, max_num) = max4([
        m[0][0] + m22_p_m33,
        m[0][0] - m22_p_m33,
       -m[0][0] + m22_m_m33,
       -m[0][0] - m22_m_m33,
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
    let coef = (tmp + tmp) / norm(q.1);
    if coef.is_infinite() {
        ZERO_VECTOR()
    } else {
        scale(coef, q.1)
    }
}

/// Product of DCM and Vector3
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, matrix_product};
/// # let PI = std::f64::consts::PI;
/// let theta = PI / 2.0;
/// let rot_x = [
///     [1.0, 0.0, 0.0],
///     [0.0, theta.cos(), -theta.sin()],
///     [0.0, theta.sin(),  theta.cos()]
/// ];
/// let v = [0.0, 1.0, 0.0];
/// 
/// let r = matrix_product(rot_x, v);
/// assert!( (r[0] - 0.0).abs() < 1e-12 );
/// assert!( (r[1] - 0.0).abs() < 1e-12 );
/// assert!( (r[2] - 1.0).abs() < 1e-12 );
/// ```
#[inline]
pub fn matrix_product<T>(m: DCM<T>, v: Vector3<T>) -> Vector3<T>
where T: Float {
    let mut r = [T::zero(); 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i] = mul_add(m[i][j], v[j], r[i]);
        }
    }
    r
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
/// Compared to `dot(a, a).sqrt()`, this function is less likely
/// to cause overflow and underflow.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, Quaternion, sum, dot, norm};
/// // --- Vector3 --- //
/// let v: Vector3<f64> = [1.0, 2.0, 3.0];
/// assert!( (14.0_f64.sqrt() - norm(v)).abs() < 1e-12 );
/// 
/// // --- Quaternion --- //
/// let q: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// assert!( (30.0_f64.sqrt() - norm(q)).abs() < 1e-12 );
/// 
/// // --- Check about overflow --- //
/// let v: Vector3<f32> = [1e15, 2e20, -3e15];
/// assert_eq!( dot(v, v).sqrt(), f32::INFINITY );  // Oh...
/// assert_eq!( norm(v), 2e20 );  // Excellent!
/// ```
#[inline]
pub fn norm<T, U>(a: U) -> T 
where T: Float, U: QuaternionOps<T> {
    a.norm()
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
    let norm_v = norm(q.1);
    let norm_q = pythag(q.0, norm_v);
    let coef = (q.0 / norm_q).acos() / norm_v;
    ( norm_q.ln(), scale(coef, q.1) )
}

// exp(q)の結果がVersorとなる条件は，qのスカラー部が0（つまりqが純虚四元数）．
// 
/// Natural logarithm of Versor.
/// 
/// If the argument `q` is guaranteed to be a Versor,
/// it is less calculation cost than the `ln(...)` function.
/// 
/// Only the vector part is returned since the real part is always zero.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, exp, ln_versor};
/// let v: Vector3<f64> = [0.1, 0.2, 0.3];
/// let r = ln_versor( exp(v) );
/// 
/// assert!( (v[0] - r[0]).abs() < 1e-12 );
/// assert!( (v[1] - r[1]).abs() < 1e-12 );
/// assert!( (v[2] - r[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn ln_versor<T>(q: Quaternion<T>) -> Vector3<T>
where T: Float {
    scale( q.0.acos() / norm(q.1), q.1)
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
    let norm_v = norm(q.1);
    let norm_q = pythag(q.0, norm_v);
    let omega = (q.0 / norm_q).acos();
    let (sin, cos) = (t * omega).sin_cos();
    let coef = norm_q.powf(t);
    ( coef * cos, scale((coef / norm_v) * sin, q.1) )
}

/// Power function of Versor.
/// 
/// If the argument `q` is guaranteed to be a Versor, 
/// it is less calculation cost than the `pow(...)` function.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Quaternion, normalize, mul, inv, pow_versor, sqrt};
/// let q: Quaternion<f64> = normalize( (1.0, [2.0, 3.0, 4.0]) );
/// 
/// let q_q = mul(q, q);
/// let q_pow_2 = pow_versor(q, 2.0);
/// assert!( (q_q.0    - q_pow_2.0).abs() < 1e-12 );
/// assert!( (q_q.1[0] - q_pow_2.1[0]).abs() < 1e-12 );
/// assert!( (q_q.1[1] - q_pow_2.1[1]).abs() < 1e-12 );
/// assert!( (q_q.1[2] - q_pow_2.1[2]).abs() < 1e-12 );
/// 
/// let q_sqrt = sqrt(q);
/// let q_pow_0p5 = pow_versor(q, 0.5);
/// assert!( (q_sqrt.0    - q_pow_0p5.0).abs() < 1e-12 );
/// assert!( (q_sqrt.1[0] - q_pow_0p5.1[0]).abs() < 1e-12 );
/// assert!( (q_sqrt.1[1] - q_pow_0p5.1[1]).abs() < 1e-12 );
/// assert!( (q_sqrt.1[2] - q_pow_0p5.1[2]).abs() < 1e-12 );
/// 
/// let q_inv = inv(q);
/// let q_pow_m1 = pow_versor(q, -1.0);
/// assert!( (q_inv.0    - q_pow_m1.0).abs() < 1e-12 );
/// assert!( (q_inv.1[0] - q_pow_m1.1[0]).abs() < 1e-12 );
/// assert!( (q_inv.1[1] - q_pow_m1.1[1]).abs() < 1e-12 );
/// assert!( (q_inv.1[2] - q_pow_m1.1[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn pow_versor<T>(q: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float {
    let (sin, cos) = (t * q.0.acos()).sin_cos();
    ( cos, scale(sin / norm(q.1), q.1) )
}

/// Square root of Quaternion.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Quaternion, mul, sqrt};
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
    let norm_v = norm(q.1);
    let norm_q = pythag(q.0, norm_v);
    let coef = ((norm_q - q.0) * half).sqrt() / norm_v;
    ( ((norm_q + q.0) * half).sqrt(), scale(coef, q.1) )
}

/// Square root of Versor.
/// 
/// If the argument `q` is guaranteed to be a Versor, 
/// it is less calculation cost than the `sqrt(...)` function.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Quaternion, normalize, mul, sqrt_versor};
/// let q: Quaternion<f64> = normalize( (1.0, [2.0, 3.0, 4.0]) );
/// let q_sqrt = sqrt_versor(q);
/// 
/// let result = mul(q_sqrt, q_sqrt);
/// assert!( (q.0    - result.0).abs() < 1e-12 );
/// assert!( (q.1[0] - result.1[0]).abs() < 1e-12 );
/// assert!( (q.1[1] - result.1[1]).abs() < 1e-12 );
/// assert!( (q.1[2] - result.1[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn sqrt_versor<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    let half = cast(0.5);
    let coef = (half - q.0 * half).sqrt() / norm(q.1);
    ( mul_add(q.0, half, half).sqrt(), scale(coef, q.1) )
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

/// Calculate the versor to rotate from vector `a` to vector `b` (Without singularity!).
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
    let norm_a = dot(a, a).sqrt();
    let norm_b = dot(b, b).sqrt();
    if norm_a == T::zero() || norm_b == T::zero() {
        return IDENTITY();
    }
    let a_u = scale(norm_a.recip(), a);
    let b_u = scale(norm_b.recip(), b);

    // 特異点回避はどこで切り替えてもいいけど，とりあえず90degで切り替える
    if dot(a_u, b_u) >= T::zero() {
        let tmp = add(a_u, b_u);
        let axis = scale(dot(tmp, tmp).sqrt().recip(), tmp);
        (T::zero(), axis)
    } else {
        let tmp = sub(a_u, b_u);
        let axis_a2mb = scale(dot(tmp, tmp).sqrt().recip(), tmp);
        let axis_mb2b = orthogonal_vector(b);
        mul(axis_mb2b, axis_a2mb)
    }
}

/// Calculate the versor to rotate from vector `a` to vector `b` by the shortest path.
/// 
/// The parameter `t` adjusts the amount of movement from `a` to `b`, 
/// so that When `t = 1`, it moves to position `b` completely.
/// 
/// If you enter a zero vector, it returns an identity quaternion.
#[inline]
pub fn rotate_a_to_b_param<T>(a: Vector3<T>, b: Vector3<T>, t: T) -> Quaternion<T>
where T: Float {
    let norm_a_b = (dot(a, a) * dot(b, b)).sqrt();
    if norm_a_b == T::zero() {
        return IDENTITY();
    }

    let axis = cross(a, b);
    let norm_axis_inv = dot(axis, axis).sqrt().recip();
    if norm_axis_inv.is_finite() {
        let cos_theta = dot(a, b) / norm_a_b;
        let (sin, cos) = (t * acos_safe(cos_theta) * cast(0.5)).sin_cos();
        (cos, scale(sin * norm_axis_inv, axis))
    } else {
        if dot(a, b) > T::zero() {
            IDENTITY()
        } else {
            (T::zero(), orthogonal_vector(a))  // theta = πの場合
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
        scale_add(t, sub(b, a), a)  // lerp
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
#[inline]
#[allow(non_snake_case)]
fn IDENTITY<T: Float>() -> Quaternion<T> {
    (T::one(), [T::zero(); 3])
}

#[inline]
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

/// （主に呼び出し側の特異点近傍で）NaNにならないように定義域外の値をカットする．
#[inline(always)]
fn acos_safe<T: Float>(x: T) -> T {
    // FloatConstを使いたくないからこの実装とする．
    (if x.abs() > T::one() { x.signum() } else { x }).acos()
}

/// 配列内の最大値とそのインデックスを返す．
#[inline]
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

/// `a`に直交し，ノルムが1であるベクトルを返す．ただし`norm(a) > 0`であること．
/// 
/// 【理論】
/// aとbが直交する場合には両者の内積がゼロになるので，そこから考えて内積が0に
/// なるようにbの要素を作れば良い．
/// 具体的には，aの要素のうち２つを入れ替えて，残り一つの要素をゼロにする．また，
/// 入れ替える要素のうち片方の符号を反転させる．
/// 例えば，`a = [0.5, 1.0, -1.5]`であれば`b = [0.0, 1.5, 1.0]`とすることで
/// aに直交なベクトルbが求まる．
/// 
/// 注意点として，入れ替える要素は片方が必ずゼロ以外で無ければならない．
/// 例えば`a = [1.0, 0.0, 0.0]`のような場合に2つめと3つめの要素を入れ替えて
/// しまうと直交なベクトルとならない．
/// 
/// この関数では，入れ替える要素が一意に決まるように絶対値が最大のものと中間の
/// ものを入れ替え，最大値の符号を反転させる．例えば，`a = [0.5, -0.8, -1.5]`
/// に対しては`b = [0.0, 1.5, -0.8]`とする．
#[inline]
fn orthogonal_vector<T: Float>(a: Vector3<T>) -> Vector3<T> {
    let mut working_array: Vector3<T> = unsafe {MaybeUninit::uninit().assume_init()};

    // aの絶対値が最大となるインデックスを探す（working_arrayにはaの絶対値を入れる）
    working_array[0] = a[0].abs();
    let mut maximum_index: usize = 0;
    let mut max_val = working_array[0];
    for (i, val) in a.iter().enumerate().skip(1) {
        working_array[i] = val.abs();
        if working_array[i] > max_val {
            max_val = working_array[i];
            maximum_index = i;
        }
    }
    // working_arrayの中央値を探す
    let idx1 = (maximum_index + 1) % 3;
    let idx2 = (maximum_index + 2) % 3;
    let median_index = if working_array[idx1] > working_array[idx2] {
        idx1
    } else {
        idx2
    };

    let norm_inv = pythag(a[median_index], a[maximum_index]).recip();
    working_array[median_index] = -a[maximum_index] * norm_inv;
    working_array[maximum_index] = a[median_index]  * norm_inv;
    working_array[3 - (median_index + maximum_index)] = T::zero();

    working_array
}

// f32,f64のメソッドには同様の機能を提供するhypot()が存在するが，
// このクレートはマイコン上での動作も想定しているため内部で平方根の
// 計算を行わないMoler-Morrison algorithmを採用した．収束が速いから
// ノルムの計算に使用しても丸め誤差が蓄積しにくいというのも理由の一つ．
// あと，単純に実装が美しい．
// 
/// Moler-Morrison algorithmにより，ピタゴラスの定理
/// `c^2 = a^2 + b^2`
/// を満たすcを求める．
/// 
/// `c = (a*a + b*b).sqrt()`のように実装した場合，有害なアンダーフロー
/// やオーバーフローが発生する可能性がある．
/// Moler-Morrison algorithmは他の方法に比べて速度面で若干劣るものの，
/// 堅牢かつ移植性の高い実装で正確な計算結果を得ることができる．
/// cがオーバーフローしない範囲の値である限り，本関数による演算で
/// オーバーフローが起こることはない．
/// 
/// 反復回数は有効数字6桁なら2回, 20桁なら3回, 60桁なら4回．
/// aとbの絶対値が大きく異なる場合（例えば|a| >> |b|）には，
/// より少ない反復回数で結果が求まる．
#[inline]
fn pythag<T: Float>(a: T, b: T) -> T {
    let mut a = a.abs();
    let mut b = b.abs();
    if a < b {
        mem::swap(&mut a, &mut b);
    }
    if b == T::zero() {
        return a;
    }

    let two : T = cast(2.0);
    let four: T = cast(4.0);
    for _ in 0..4 {  // loopにするとinfやNaNが入った際に抜けられなくなる
        let mut s = b / a;
        s = s * s;
        let tmp = four + s;
        if tmp == four {  // 収束判定（これでちゃんとbreakできる）
            break;
        }
        s = s / tmp;
        a = mul_add(two, a * s, a);
        b = s * b;
    }
    a
}


#[cfg(test)]
mod test_private_functions {
    use super::*;

    #[test]
    fn test_max4() {
        let array = [-0.5, 2.0, 3.0, 0.1];
        let (idx, val) = max4(array);
        assert_eq!(idx, 2);
        assert!((val - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_orthogonal_vector() {
        let a = normalize([1.0, -5.0, 3.0]);
        let b = orthogonal_vector(a);
        assert!(dot(a, b).abs() < 1e-12);
        assert!((1.0 - norm(b)).abs() < 1e-12);
    }

    #[test]
    fn test_pythag() {
        let a = 13.0f64.sqrt();
        let b = pythag(-2.0, 3.0);
        assert!((a - b).abs() < 1e-12);
    }
}