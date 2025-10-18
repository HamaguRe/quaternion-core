//! Quaternion library written in Rust. 
//! 
//! This crate provides Quaternion operations and conversions between several rotation 
//! representations. Functions are generic and support `f32` and `f64`.
//! 
//! ## Generics
//! 
//! Functions implementing the `QuaternionOps` trait can take both `Quaternion<T>` 
//! and `Vector3<T>` types as arguments (How convenient...!!).
//! 
//! For example:
//! ```
//! use quaternion_core::{Vector3, Quaternion, add};
//! 
//! // --- Vector3 --- //
//! let v1: Vector3<f32> = [1.0, 2.0, 3.0];
//! let v2: Vector3<f32> = [0.1, 0.2, 0.3];
//! 
//! let added_v = add(v1, v2);  // <--- [1.1, 2.2, 3.3]
//! 
//! // --- Quaternion --- //
//! let q1: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
//! let q2: Quaternion<f64> = (0.1, [0.2, 0.3, 0.4]);
//! 
//! let added_q = add(q1, q2);  // <--- (1.1, [2.2, 3.3, 4.4])
//! ```
//! 
//! ## What is Versor?
//! 
//! Versor denotes a Quaternion representing rotation, with a norm of 1. 
//! This is precisely what is meant by a Unit Quaternion.
//! 
//! The documentation for this crate basically writes Versor instead of Unit Quaternion, 
//! but the difference in usage is not clear.
//! Please think Versor = Unit Quaternion.

#![no_std]
#[cfg(feature = "std")]
extern crate std;

#[cfg(any(feature = "std", feature = "libm"))]
use num_traits::float::{Float, FloatConst};

mod euler;
mod generics;
mod pfs;
pub use generics::QuaternionOps;

/// Three dimensional vector (Pure Quaternion)
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
/// ```ignore
/// [
///     [m11, m12, m13],
///     [m21, m22, m23],
///     [m31, m32, m33]
/// ]
/// ```
pub type DCM<T> = [Vector3<T>; 3];

/// Specifies the rotation type of Euler angles.
/// 
/// Considering a fixed `Reference frame` and a rotating `Body frame`, 
/// `Intrinsic rotation` and `Extrinsic rotation` represent the following rotations:
/// 
/// * `Intrinsic`: Rotate around the axes of the `Body frame`
/// * `Extrinsic`: Rotate around the axes of the `Reference frame`
#[cfg_attr(feature = "serde-serialize", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy)]
pub enum RotationType {
    Intrinsic,
    Extrinsic,
}

/// Specifies the rotation sequence of Euler angles.
/// 
/// Each variant reads from left to right.
/// For example, `RotationSequence::XYZ` represents first a rotation 
/// around the X axis, then around the Y axis, and finally around 
/// the Z axis (X ---> Y ---> Z).
#[cfg_attr(feature = "serde-serialize", derive(serde::Serialize, serde::Deserialize))]
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

/// Generate identity quaternion.
/// 
/// It returns `(1.0, [0.0, 0.0, 0.0])`.
#[inline]
pub fn identity<T>() -> Quaternion<T>
where T: Float {
    (T::one(), [T::zero(); 3])
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
/// // rotation of pi/2[rad] around the y-axis.
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
    let theta = angle % ( T::PI() + T::PI() );  // Range reduction: (-2π, 2π)
    let (sin, cos) = ( theta * pfs::cast(0.5) ).sin_cos();
    let coef = sin / norm(axis);
    if coef.is_infinite() {
        identity()
    } else {
        ( cos, scale(coef, axis) )
    }
}

/// Extracts the rotation `axis` and the rotation `angle` around the `axis` from the Versor.
/// 
/// # Returns
///
/// The function returns a tuple `(axis, angle)`, where:
///
/// * **`axis`**: The rotation axis as a **unit vector** (`Vector3`).
/// * **`angle`**: The rotation angle in **radians**. The range is `(-PI, PI]`.
///
/// ## Special Case: Identity Quaternion
///
/// If the input is the **Identity Quaternion** `(1.0, [0.0, 0.0, 0.0])`, 
/// the function returns an angle of zero and a zero axis vector.
///
/// # Singularity
///
/// Because this function returns a normalized rotation axis, 
/// when the norm of the vector part of Versor is zero, a singularity 
/// occurs and accuracy decreases.
/// 
/// If you want to use the calculated `angle` and `axis` as `scale(angle, axis)`, 
/// it is better to use the `to_rotation_vector` function. 
/// The `to_rotation_vector` function can be calculated without singularities.
///
/// # Examples
/// 
/// ```
/// # use quaternion_core::{from_axis_angle, to_axis_angle, sub, normalize};
/// # let PI = std::f64::consts::PI;
/// let axis_ori = [0.0, 1.0, 2.0];
/// let angle_ori = PI / 2.0;
/// let q = from_axis_angle(axis_ori, angle_ori);
/// 
/// let (axis, angle) = to_axis_angle(q);
/// 
/// let diff = sub( normalize(axis_ori), normalize(axis) );
/// assert!( diff[0].abs() < 1e-12 );
/// assert!( diff[1].abs() < 1e-12 );
/// assert!( diff[2].abs() < 1e-12 );
/// assert!( (angle_ori - angle).abs() < 1e-12 );
/// ```
#[inline]
pub fn to_axis_angle<T>(q: Quaternion<T>) -> (Vector3<T>, T)
where T: Float {
    let norm_v = norm(q.1);
    let norm_inv = norm_v.recip();
    if norm_inv.is_infinite() {
        ( [T::zero(); 3], T::zero() )
    } else {
        let mut angle = pfs::cast::<T>(2.0) * norm_v.min( T::one() ).asin();
        if q.0 < T::zero() {
            angle = -angle;
        }
        ( scale(norm_inv, q.1), angle )
    }
}

/// Converts a **Direction Cosine Matrix (DCM)** into a **Versor**
/// 
/// Crucially, this function assumes the input DCM represents a
/// Point Rotation (Frame Fixed), where a vector `v` is rotated by the operation `q v q*`.
///
/// If your input DCM represents a **Frame Rotation (Point Fixed)**,
/// which corresponds to the rotation `q* v q`, you must take the conjugate
/// of the resulting Versor to get the correct rotation:
///
/// ```
/// # use quaternion_core::{from_dcm, conj, to_dcm};
/// # let dcm = to_dcm(conj((1.0, [0.0; 3])));
/// // q_pr represents the Versor for Point Rotation.
/// let q_pr = from_dcm(dcm); 
/// // q_fr is the Versor for the Frame Rotation.
/// let q_fr = conj(q_pr);
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
    let (index, max_num) = pfs::max4([
        m[0][0] + m22_p_m33,
        m[0][0] - m22_p_m33,
       -m[0][0] + m22_m_m33,
       -m[0][0] - m22_m_m33,
    ]);

    let half: T = pfs::cast(0.5);
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

/// Converts a **Versor** into a **Direction Cosine Matrix (DCM)**.
///
/// **By default, the output DCM represents a Point Rotation (Frame Fixed)**,
/// which rotates a vector `v` by the quaternion operation `q v q*`.
///
/// If you need a DCM that represents a **Frame Rotation (Point Fixed)**
/// (the rotation `q* v q`), you must take the conjugate of the Versor:
///
/// ```
/// # use quaternion_core::{to_dcm, conj};
/// # let q = (1.0, [0.0; 3]);
/// // This DCM represents the Frame Rotation.
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
    let two = pfs::cast(2.0);

    // Compute these value only once.
    let (q0_q0, [q0_q1, q0_q2, q0_q3]) = scale(q.0, q);
    let q1_q2 = q.1[0] * q.1[1];
    let q1_q3 = q.1[0] * q.1[2];
    let q2_q3 = q.1[1] * q.1[2];

    let m11 = pfs::mul_add(pfs::mul_add(q.1[0], q.1[0], q0_q0), two, neg_one);
    let m12 = (q1_q2 - q0_q3) * two;
    let m13 = (q1_q3 + q0_q2) * two;
    let m21 = (q1_q2 + q0_q3) * two;
    let m22 = pfs::mul_add(pfs::mul_add(q.1[1], q.1[1], q0_q0), two, neg_one);
    let m23 = (q2_q3 - q0_q1) * two;
    let m31 = (q1_q3 - q0_q2) * two;
    let m32 = (q2_q3 + q0_q1) * two;
    let m33 = pfs::mul_add(pfs::mul_add(q.1[2], q.1[2], q0_q0), two, neg_one);

    [
        [m11, m12, m13],
        [m21, m22, m23],
        [m31, m32, m33],
    ]
}

/// Converts a **Euler Angles** into a **Versor**.
///
/// This function requires two parameters to fully define the rotation:
///
/// 1.  `RotationType`: Specifies whether the rotation is **Intrinsic** or **Extrinsic**.
/// 2.  `RotationSequence`: Defines the three-axis sequence (e.g., XYZ, ZYX, XZX, ...).
///
/// The input `angles` array must contain the three angles in **radians**,
/// corresponding to the specified rotation sequence: `angles[0]` -> `angles[1]` -> `angles[2]`.
/// (Each angle should be within the range: `[-2*PI, 2*PI]`).
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
where T: Float {
    match rt {
        RotationType::Intrinsic => euler::from_intrinsic_euler_angles(rs, angles),
        RotationType::Extrinsic => euler::from_extrinsic_euler_angles(rs, angles),
    }
}

/// Converts a **Versor** into a **Euler Angles**.
///
/// This function requires two parameters to fully define the rotation:
///
/// 1.  `RotationType`: Specifies whether the rotation is **Intrinsic** or **Extrinsic**.
/// 2.  `RotationSequence`: Defines the three-axis sequence (e.g., XYZ, ZYX, XZX, ...).
///
/// The output `angles` array contains the three angles, corresponding to the sequence:
/// `angles[0]` -> `angles[1]` -> `angles[2]`.
/// Each angle is returned in the range: `(-PI, PI]`.
///
/// # Singularity (Gimbal Lock) 
///
/// ## RotationType::Intrinsic
/// 
/// For Proper Euler Angles (ZXZ, XYX, YZY, ZYZ, XZX, YXY), the singularity is reached 
/// when the sine of the second rotation angle is 0 (angle = 0, ±π, ...), and for 
/// Tait-Bryan angles (XYZ, YZX, ZXY, XZY, ZYX, YXZ), the singularity is reached when 
/// the cosine of the second rotation angle is 0 (angle = ±π/2).
/// 
/// ## RotationType::Extrinsic
/// 
/// As in the case of Intrinsic rotation, for Proper Euler Angles, the singularity occurs 
/// when the sine of the second rotation angle is 0 (angle = 0, ±π, ...), and for 
/// Tait-Bryan angles, the singularity occurs when the cosine of the second rotation angle 
/// is 0 (angle = ±π/2).
/// 
/// ## Resolution at Singularity
/// * For **Intrinsic** rotation, the **third angle** (`angles[2]`) is set to 0 \[rad\].
/// * For **Extrinsic** rotation, the **first angle** (`angles[0]`) is set to 0 \[rad\].
///
/// # Examples
/// 
/// Depending on the rotation angle of each axis, it may not be possible to recover 
/// the same rotation angle as the original. However, they represent the same rotation in 3D space.
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

/// Converts a **Rotation Vector** into a **Versor**.
///
/// A Rotation Vector is a convenient representation where:
/// 
/// 1.  Its **direction** defines the **rotation axis**.
/// 2.  Its **norm** defines the **rotation angle** (in radians).
///
/// There are no particular restrictions on the norm of the input Rotation Vector.
/// Also, even if a zero vector is input, the conversion to a Versor can be performed 
/// without falling into a singularity.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{from_rotation_vector, scale, point_rotation};
/// # let PI = std::f64::consts::PI;
/// let angle = PI / 2.0;
/// let axis = [1.0, 0.0, 0.0];
/// 
/// // This represents a rotation of pi/2 around the x-axis.
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
pub fn from_rotation_vector<T>(r: Vector3<T>) -> Quaternion<T>
where T: Float + FloatConst {
    let half: T = pfs::cast(0.5);

    let theta = norm(r) % ( T::PI() + T::PI() );  // Range reduction: [0, 2π)
    let half_theta = half * theta;
    (half_theta.cos(), scale(half * pfs::sinc(half_theta), r))
}

/// Converts a **Versor** into a **Rotation Vector**.
/// 
/// A Rotation Vector is a convenient representation where:
/// 
/// 1.  Its **direction** defines the **rotation axis**.
/// 2.  Its **norm** defines the **rotation angle** (in radians).
/// 
/// The resulting vector's norm (the rotation angle) is always constrained
/// to the range: `[0, PI]`
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
    // half_thetaをasinで求めるとhalf_theta=±pi/2のときに精度低下するのでatanを使う．
    let half_theta = (norm(q.1) / q.0).atan();  // [-pi/2, pi/2]
    let coef = pfs::cast::<T>(2.0) / pfs::sinc(half_theta);   // [2, pi]
    scale(coef.copysign(q.0), q.1)  // sinc関数は偶関数なので，half_thetaの符号をここで反映する必要がある．
}

/// Product of **Direction Cosine Matrix (DCM)** and **Vector3**
/// 
/// This is the product of a 3x3 matrix and a 3D vector.
/// It is used to apply the rotation represented by the DCM to the vector.
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
    [
        dot(m[0], v),
        dot(m[1], v),
        dot(m[2], v),
    ]
}

/// Sum all elements of a Quaternion or Vector3.
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

/// Addition of two Quaternions or Vector3s: `a + b`
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

/// Subtraction of two Quaternions or Vector3s: `a - b`
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

/// Scaling a Quaternion or Vector3 by a scalar factor: `s * a`
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

/// Scaling and addition in one step: `s * a + b`
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

/// Calculate the element-wise product of two Quaternions or Vector3s.
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

/// Hadamard product and addition in one step.
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

/// Dot product of two Quaternions or Vector3s: `a · b`
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
    a.dot(b)
}

/// Cross product of two Quaternions or Vector3s: `a × b`
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

/// Calculate L2 norm of Vector3s or Quaternions.
/// 
/// Compared to `dot(a, a).sqrt()`, this function is less likely
/// to cause overflow and underflow.
/// 
/// When the `norm-sqrt` feature is enabled, the default 
/// implementation is compiled with `dot(a, a).sqrt()` instead.
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
/// 
/// #[cfg(not(feature = "norm-sqrt"))]
/// assert_eq!( norm(v), 2e20 );  // Excellent!
/// ```
#[inline]
pub fn norm<T, U>(a: U) -> T 
where T: Float, U: QuaternionOps<T> {
    a.norm()
}

/// Normalize a Vector3 or Quaternion.
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
    scale(norm(a).recip(), a)
}

/// Invert the sign of a Vector3s or Quaternions.
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

/// Multiplication of two Quaternions or Pure Quaternoins (Vector3).
/// 
/// This is an operation also known as the Hamilton product.
/// In this operation, the Vector3 is treated as a pure quaternoin (which is a quaternion with zero real part).
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

/// Calculate the conjugate of a Quaternion.
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

/// Calculate the inverse of a Quaternion or Pure Quaternion (Vector3).
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

// acosは[-pi/2, pi/2]の範囲でしか値を返さないので、qのとり方によってはlnで完全に復元できない。
// q == ln(exp(q)) が成り立つのはcos(norm(q.1))が[-pi/2, pi/2]の範囲内にある場合のみ。
//
/// Exponential function of a Quaternion or Pure Quaternion (Vector3).
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

// acosは[-pi/2, pi/2]の範囲でしか値を返さないので、qのとり方によってはlnで完全に復元できない。
// q == ln(exp(q)) が成り立つのはcos(norm(q.1))が[-pi/2, pi/2]の範囲内にある場合のみ。
// 
/// Natural logarithm of a Quaternion.
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
    let norm_q = pfs::norm2(q.0, norm_v);
    let coef = (q.0 / norm_q).acos() / norm_v;
    ( norm_q.ln(), scale(coef, q.1) )
}

// exp(q)の結果がVersorとなる条件は，qのスカラー部が0（つまりqが純虚四元数）．
// 
/// Natural logarithm of a Versor.
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

/// Power function of a Quaternion.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Quaternion, mul, inv, pow, sqrt};
/// let q: Quaternion<f64> = (1.0, [2.0, 3.0, 4.0]);
/// 
/// let q_pow_0 = pow(q, 0.0);
/// assert!( (1.0 - q_pow_0.0).abs() < 1e-12 );
/// assert!( (0.0 - q_pow_0.1[0]).abs() < 1e-12 );
/// assert!( (0.0 - q_pow_0.1[1]).abs() < 1e-12 );
/// assert!( (0.0 - q_pow_0.1[2]).abs() < 1e-12 );
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
    let norm_q = pfs::norm2(q.0, norm_v);
    let omega = (norm_v / q.0).atan();
    let coef = norm_q.powf(t);
    let tmp = t * omega;
    let numer = coef * t * pfs::sinc(tmp);
    let denom = norm_q * pfs::sinc(omega);
    ( coef * tmp.cos(), scale(numer / denom, q.1) )
}

/// Power function of a Versor.
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
/// let q_pow_0 = pow_versor(q, 0.0);
/// assert!( (1.0 - q_pow_0.0).abs() < 1e-12 );
/// assert!( (0.0 - q_pow_0.1[0]).abs() < 1e-12 );
/// assert!( (0.0 - q_pow_0.1[1]).abs() < 1e-12 );
/// assert!( (0.0 - q_pow_0.1[2]).abs() < 1e-12 );
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
    let omega = ( norm(q.1) / q.0 ).atan();
    let tmp = t * omega;
    ( tmp.cos(), scale(t * pfs::sinc(tmp) / pfs::sinc(omega), q.1) )
}

/// Square root of a Quaternion.
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
    let half = pfs::cast(0.5);
    let norm_v = norm(q.1);
    let norm_q = pfs::norm2(q.0, norm_v);
    let omega = (norm_v / q.0).atan();
    let coef = half * pfs::sinc(omega*half) / (norm_q.sqrt() * pfs::sinc(omega));
    ( ((norm_q + q.0) * half).sqrt(), scale(coef, q.1) )
}

/// Square root of a Versor.
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
    let half = pfs::cast(0.5);
    let omega = (norm(q.1) / q.0).atan();
    let coef = half * pfs::sinc(omega*half) / pfs::sinc(omega);
    ( pfs::mul_add(q.0, half, half).sqrt(), scale(coef, q.1) )
}

/// Rotates a point by the Versor (Point Rotation - Frame Fixed)
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
    scale_add(pfs::cast(2.0), cross(q.1, tmp), v)
}

/// Rotates a frame by the Versor (Frame Rotation - Point Fixed)
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
    scale_add(pfs::cast(2.0), cross(tmp, q.1), v)
}

/// Calculate the Versor to rotate from vector `a` to vector `b` (Without singularity!).
/// 
/// This function calculates `q` satisfying `b = point_rotation(q, a)` when `norm(a) = norm(b)`.
/// 
/// # Characteristics
/// 
/// * **Robustness:** This function can accurately calculate the versor regardless of the 
/// combination of vector orientations (however, `norm(a) > 0` and `norm(b) > 0`).
/// * **Axis Ambiguity:** This function provides **no guarantees** regarding the direction 
/// or regularity of the rotation axis. If you require a rotation axis that is **orthogonal** 
/// to both vector `a` and vector `b`, you should use `rotate_a_to_b_shortest` function.
/// 
/// # Returns
/// 
/// Returns `None` if either input vector `a` or `b` is a zero vector,
/// as a rotation cannot be uniquely defined in that case.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, cross, rotate_a_to_b, point_rotation, normalize};
/// let a: Vector3<f64> = [1.5, -0.5, 0.2];
/// let b: Vector3<f64> = [0.1, 0.6, 1.0];
/// 
/// let q = rotate_a_to_b(a, b).unwrap();
/// let b_check = point_rotation(q, a);
/// 
/// let b_u = normalize(b);
/// let b_check_u = normalize(b_check);
/// assert!( (b_u[0] - b_check_u[0]).abs() < 1e-12 );
/// assert!( (b_u[1] - b_check_u[1]).abs() < 1e-12 );
/// assert!( (b_u[2] - b_check_u[2]).abs() < 1e-12 );
/// ```
#[inline]
pub fn rotate_a_to_b<T>(a: Vector3<T>, b: Vector3<T>) -> Option<Quaternion<T>>
where T: Float {
    let norm_inv_a = norm(a).recip();
    let norm_inv_b = norm(b).recip();
    if norm_inv_a.is_infinite() || norm_inv_b.is_infinite() {
        return None;
    }
    let a = scale(norm_inv_a, a);
    let b = scale(norm_inv_b, b);

    let a_add_b = add(a, b);
    let dot_a_add_b = dot(a_add_b, a_add_b);
    if dot_a_add_b >= T::one() {  // arccos(a・b) = 120degで切り替え
        Some( (T::zero(), scale(dot_a_add_b.sqrt().recip(), a_add_b)) )
    } else {
        let norm_a_sub_b = (pfs::cast::<T>(4.0) - dot_a_add_b).sqrt();
        let axis_a2mb = scale(norm_a_sub_b.recip(), sub(a, b));  // a ==> -b
        let axis_mb2b = pfs::orthogonal_vector(b);  // -b ==> b
        Some( mul(axis_mb2b, axis_a2mb) )
    }
}

/// Calculate the versor to rotate from vector `a` to vector `b` by the shortest path.
/// 
/// This function calculates `q` satisfying `b = point_rotation(q, a)` when `norm(a) = norm(b)`.
///
/// # Characteristics
///
/// * **Shortest Path:** This method guarantees the rotation axis is always **orthogonal** to 
/// both vector **a** and vector **b**, representing the geometrically shortest path between 
/// the two directions.
/// * **Rotation Angle:** The rotation angle is in the range `[0, PI]` radians.
/// * **Parallel Vectors:** If vector `a` and `b` are **opposite** and parallel, the 
/// rotation axis is theoretically ambiguous. This function provides a valid rotation, 
/// but the axis direction is not guaranteed (although it will be orthogonal to the `a` and `b`).
///
/// # Performance Consideration
///
/// This function is slightly more computationally intensive than `rotate_a_to_b`, 
/// especially when the angle between the vectors is near PI. If an orthogonal rotation 
/// axis is not strictly required, `rotate_a_to_b` may be preferred for performance.
///
/// # Returns
///
/// Returns `None` if either input vector **a** or **b** is a zero vector,
/// as a rotation cannot be uniquely defined in that case.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, scale, cross, to_rotation_vector, from_rotation_vector, rotate_a_to_b_shortest, point_rotation, normalize};
/// let a: Vector3<f64> = [1.5, -0.5, 0.2];
/// let b: Vector3<f64> = [0.1, 0.6, 1.0];
/// 
/// let q = rotate_a_to_b_shortest(a, b).unwrap();
/// let b_check = point_rotation(q, a);
/// 
/// let b_u = normalize(b);
/// let b_check_u = normalize(b_check);
/// assert!( (b_u[0] - b_check_u[0]).abs() < 1e-12 );
/// assert!( (b_u[1] - b_check_u[1]).abs() < 1e-12 );
/// assert!( (b_u[2] - b_check_u[2]).abs() < 1e-12 );
/// 
/// // --- If the amount of displacement of vector `a` is to be adjusted --- //
/// // The parameter `t` adjusts the amount of movement from `a` to `b`. 
/// // When `t = 1`, `a` moves completely to position `b`.
/// let t = 0.5;
/// let mut q = rotate_a_to_b_shortest(a, b).unwrap();
/// let r = to_rotation_vector(q);  // To avoid singularities, proceed via the rotation vector.
/// q = from_rotation_vector(scale(t, r));
/// ```
#[inline]
pub fn rotate_a_to_b_shortest<T>(a: Vector3<T>, b: Vector3<T>) -> Option<Quaternion<T>>
where T: Float {
    let norm_inv_a = norm(a).recip();
    let norm_inv_b = norm(b).recip();
    if norm_inv_a.is_infinite() || norm_inv_b.is_infinite() {
        return None;
    }
    let a = scale(norm_inv_a, a);
    let b = scale(norm_inv_b, b);

    let half = pfs::cast::<T>(0.5);
    let a_dot_b = dot(a, b);  // == cos(theta)
    if a_dot_b > pfs::cast(-0.94) {  // theta < 160 deg
        let q_s = (half + half * a_dot_b).sqrt();  // == cos(theta/2)
        Some( (q_s, scale(half / q_s, cross(a, b))) )  // a --> b
    } else {
        let b_cross_a = cross(b, a);
        let q_s = (half - half * a_dot_b).sqrt();  // == cos((pi-theta)/2)
        let q_a2mb = (q_s, scale(half / q_s, b_cross_a));  // a --> -b
        let r_o = pfs::orthogonal_vector(b);

        let c = sub(b_cross_a, r_o);
        let r_o_cross_c = cross(r_o, c);

        // rho < pi/2
        let norm_c = norm(c);  // norm_c > 0
        let norm_c_sin_rho = norm(r_o_cross_c).copysign( dot(r_o_cross_c, b) );  // norm(c) * sin(rho)
        let norm_c_cos_rho = (norm_c * norm_c - norm_c_sin_rho * norm_c_sin_rho).sqrt();  // norm(c) * cos(rho)

        // r_o と b_cross_a の間の角度lambdaを求める
        let lambda = norm_c_sin_rho.atan2(T::one() - norm_c_cos_rho);  // (-pi, pi]

        let (sin, cos) = (half * lambda).sin_cos();
        let q_b = (cos, scale(sin, b));
        let r = point_rotation(q_b, r_o);  // -b --> a

        let q_a2b_s = -dot(r, q_a2mb.1);
        let q_a2b_v = scale_add(q_a2mb.0, r, cross(r, q_a2mb.1));
        Some( (q_a2b_s, q_a2b_v) )  // a --> b
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
    // 最短経路で補間する
    if dot(a, b).is_sign_negative() {
        // bの符号を反転
        scale_add(-t, add(a, b), a)
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
    // 最短経路で補間する
    let mut dot = dot(a, b);
    if dot.is_sign_negative() {
        b = negate(b);
        dot = -dot;
    }
    // If the distance between quaternions is close enough, use lerp.
    if dot > pfs::cast(0.9995) {  // Approximation error < 0.017%
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
