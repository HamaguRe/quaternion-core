//! Quaternion Libraly (f32 & f64)

#![no_std]
#[cfg(feature = "std")]
extern crate std;

use num_traits::{Float, FloatConst};

mod simd;
mod euler;
pub use simd::FloatSimd;

/// `[i, j, k]`
pub type Vector3<T> = [T; 3];

/// `(1, [i, j, k])`
pub type Quaternion<T> = (T, Vector3<T>);

/// Direction Cosine Matrix
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
/// * `Intrinsic`: Rotate around the axes of the `Body-frame`
/// * `Extrinsic`: Rotate around the axes of the `Reference-frame`
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
/// The `axis` does not have to be a unit vector.
/// 
/// If you enter a zero vector, it returns an identity quaternion.
#[inline]
pub fn from_axis_angle<T>(axis: Vector3<T>, angle: T) -> Quaternion<T>
where T: Float + FloatConst {
    let theta = angle % ( T::PI() + T::PI() );  // limit to (-2π, 2π)
    let f = ( theta * cast(0.5) ).sin_cos();
    let coef = f.0 / norm_vec(axis);
    if coef.is_infinite() {
        IDENTITY()
    } else {
        ( f.1, scale_vec(coef, axis) )
    }
}

/// Compute the rotation `axis` (unit vector) and the rotation `angle`\[rad\] 
/// around the axis from the versor.
/// 
/// If identity quaternion is entered, `angle` returns zero and 
/// the `axis` returns a zero vector.
/// 
/// Range of `angle`: `(-π, π]`
#[inline]
pub fn to_axis_angle<T>(q: Quaternion<T>) -> (Vector3<T>, T)
where T: Float {
    let norm_q_v = norm_vec(q.1);
    let coef = norm_q_v.recip();
    if coef.is_infinite() {
        ( ZERO_VECTOR(), T::zero() )
    } else {
        // 少しの誤差は見逃す．
        let tmp = norm_q_v.min( T::one() ).asin();
        ( scale_vec(coef, q.1), (tmp + tmp).copysign(q.0) ) // theta = 2*tmp
    }
}

/// Convert a DCM to a Quaternion representing 
/// the `q v q*` rotation (Point Rotation - Frame Fixed).
/// 
/// When convert to a DCM representing `q* v q` rotation
/// (Frame Rotation - Point Fixed) to a Quaternion, do the following:
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
/// # let pi = std::f64::consts::PI;
/// // Make these as you like.
/// let v = [1.0, 0.5, -8.0];
/// let q = from_axis_angle([0.2, 1.0, -2.0], pi/4.0);
/// 
/// // --- Point rotation --- //
/// {
///     let m = to_dcm(q);
///     let q_check = from_dcm(m);
///     
///     assert!( (q.0 - q_check.0).abs() < 1e-12 );
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
///     assert!( (q.0 - q_check.0).abs() < 1e-12 );
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

/// Convert a Quaternion to a DCM representing 
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
/// ```
/// # use quaternion_core::{
/// #     from_axis_angle, to_dcm, conj, 
/// #     matrix_product, point_rotation, frame_rotation
/// # };
/// # let pi = std::f64::consts::PI;
/// // Make these as you like.
/// let v = [1.0, 0.5, -8.0];
/// let q = from_axis_angle([0.2, 1.0, -2.0], pi/4.0);
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
where T: Float + FloatSimd<T> {
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

/// Convert Euler angles to Quaternion.
/// 
/// The type of rotation (Intrinsic or Extrinsic) is specified by `RotationType` enum, 
/// and the rotation sequence (XZX, XYZ, ...) is specified by `RotationSequence` enum.
/// 
/// Each element of `angles` should be specified in the range `[-2π, 2π]`.
/// 
/// Sequences: `angles[0]` ---> `angles[1]` ---> `angles[2]`
/// 
/// # Example
/// 
/// ```
/// # use quaternion_core::{from_axis_angle, mul, from_euler_angles, point_rotation};
/// # let pi = std::f64::consts::PI;
/// use quaternion_core::{RotationType::*, RotationSequence::XYZ};
/// 
/// let angles = [pi/6.0, 1.6*pi, -pi/4.0];
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

/// Convert Quaternion (Unit quaternion) to Euler angles.
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
/// # Example
/// 
/// Depending on the rotation angle of each axis, it may not be possible to recover the 
/// same rotation angle as the original. However, they represent the same rotation in 3D space.
/// 
/// ```
/// # use quaternion_core::{from_euler_angles, to_euler_angles, point_rotation};
/// # let pi = std::f64::consts::PI;
/// use quaternion_core::{RotationType::*, RotationSequence::XYZ};
/// 
/// let angles = [pi/6.0, pi/4.0, pi/3.0];
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
where T: Float + FloatConst + FloatSimd<T> {
    match rt {
        RotationType::Intrinsic => euler::to_intrinsic_euler_angles(rs, q),
        RotationType::Extrinsic => euler::to_extrinsic_euler_angles(rs, q),
    }
}

/// Convert Rotation vector to Quaternion (Unit quaternion).
/// 
/// The Rotation vector itself represents the axis of rotation, 
/// and the norm represents the angle of rotation around the axis.
/// 
/// `angle` range is: `(0, 2π)`
#[inline]
pub fn from_rotation_vector<T>(v: Vector3<T>) -> Quaternion<T>
where T: Float {
    let theta = norm_vec(v);
    let f = ( theta * cast(0.5) ).sin_cos();
    let coef = f.0 / theta;
    if coef.is_infinite() {
        IDENTITY()
    } else {
        ( f.1, scale_vec(coef, v) )
    }
}

/// Convert Versor to rotation vector.
/// 
/// The rotation vector itself represents the axis of rotation, 
/// and the norm represents the angle of rotation around the axis.
/// 
/// `angle` range is: `(0, 2π)`
#[inline]
pub fn to_rotation_vector<T>(q: Quaternion<T>) -> Vector3<T>
where T: Float {
    let tmp = acos_safe(q.0);
    let coef = (tmp + tmp) / norm_vec(q.1);  // 2*tmp
    if coef.is_infinite() {
        ZERO_VECTOR()
    } else {
        scale_vec(coef, q.1)
    }
}

/// Product of matrix and vector
/// 
/// Rotate vectors using a directional cosine matrix.
#[inline]
pub fn matrix_product<T>(m: DCM<T>, v: Vector3<T>) -> Vector3<T>
where T: Float {
    [
        dot_vec(m[0], v),
        dot_vec(m[1], v),
        dot_vec(m[2], v),
    ]
}

/// Calculate the sum of each element of Vector3.
#[inline]
pub fn sum_vec<T>(v: Vector3<T>) -> T
where T: Float {
    v[0] + v[1] + v[2]
}

/// Calculate the sum of each element of Quaternion.
#[inline]
pub fn sum<T>(q: Quaternion<T>) -> T
where T: FloatSimd<T> {
    T::sum(q)
}

/// Calculate `a + b`
#[inline]
pub fn add_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T>
where T: Float {
    [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ]
}

/// Calculate `a + b`
#[inline]
pub fn add<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>
where T: FloatSimd<T> {
    T::add(a, b)
}

/// Calculate `a - b`
#[inline]
pub fn sub_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T>
where T: Float {
    [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ]
}

/// Calculate `a - b`
#[inline]
pub fn sub<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>
where T: FloatSimd<T> {
    T::sub(a, b)
}

/// Calculate `s * v`
/// 
/// Multiplication of scalar and vector.
#[inline]
pub fn scale_vec<T>(s: T, v: Vector3<T>) -> Vector3<T>
where T: Float {
    [ s*v[0], s*v[1], s*v[2] ]
}

/// Calculate `s * q`
/// 
/// Multiplication of scalar and quaternion.
#[inline]
pub fn scale<T>(s: T, q: Quaternion<T>) -> Quaternion<T>
where T: FloatSimd<T> {
    T::scale(s, q)
}

/// Calculate `s*a + b`
/// 
/// If the `fma` feature is enabled, the FMA calculation is performed using the `mul_add` method. 
/// If not enabled, it's computed by unfused multiply-add (s*a + b).
#[inline]
pub fn scale_add_vec<T>(s: T, a: Vector3<T>, b: Vector3<T>) -> Vector3<T>
where T: Float {
    [
        mul_add(s, a[0], b[0]),
        mul_add(s, a[1], b[1]),
        mul_add(s, a[2], b[2]),
    ]
}

/// Calculate `s*a + b`
/// 
/// If the `fma` feature is enabled, the FMA calculation is performed using the `mul_add` method. 
/// If not enabled, it's computed by unfused multiply-add (s*a + b).
#[inline]
pub fn scale_add<T>(s: T, a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>
where T: FloatSimd<T> {
    T::scale_add(s, a, b)
}

/// Hadamard product of Vector.
/// 
/// Calculate `a ∘ b`
#[inline]
pub fn hadamard_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T>
where T: Float {
    [ a[0]*b[0], a[1]*b[1], a[2]*b[2] ]
}

/// Hadamard product of Quaternion.
/// 
/// Calculate `a ∘ b`
#[inline]
pub fn hadamard<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>
where T: FloatSimd<T> {
    T::hadamard(a, b)
}

/// Hadamard product and Addiction of Vector.
/// 
/// Calculate `a ∘ b + c`
/// 
/// If the `fma` feature is enabled, the FMA calculation is performed using the `mul_add` method. 
/// If not enabled, it's computed by unfused multiply-add (s*a + b).
#[inline]
pub fn hadamard_add_vec<T>(a: Vector3<T>, b: Vector3<T>, c: Vector3<T>) -> Vector3<T>
where T: Float {
    [
        mul_add(a[0], b[0], c[0]),
        mul_add(a[1], b[1], c[1]),
        mul_add(a[2], b[2], c[2]),
    ]
}

/// Hadamard product and Addiction of Quaternion.
/// 
/// Calculate `a ∘ b + c`
/// 
/// If the `fma` feature is enabled, the FMA calculation is performed using the `mul_add` method. 
/// If not enabled, it's computed by unfused multiply-add (s*a + b).
#[inline]
pub fn hadamard_add<T>(a: Quaternion<T>, b: Quaternion<T>, c: Quaternion<T>) -> Quaternion<T>
where T: FloatSimd<T> {
    T::hadamard_add(a, b, c)
}

/// Dot product of vector.
#[inline]
pub fn dot_vec<T>(a: Vector3<T>, b: Vector3<T>) -> T
where T: Float {
    sum_vec( hadamard_vec(a, b) )
}

/// Dot product of quaternion.
#[inline]
pub fn dot<T>(a: Quaternion<T>, b: Quaternion<T>) -> T 
where T: FloatSimd<T> {
    T::dot(a, b)
}

/// Cross product.
#[inline]
pub fn cross_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T>
where T: Float {
    [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    ]
}

/// Calculate L2 norm.
#[inline]
pub fn norm_vec<T>(v: Vector3<T>) -> T
where T: Float {
    dot_vec(v, v).sqrt()
}

/// Calculate L2 norm.
#[inline]
pub fn norm<T>(q: Quaternion<T>) -> T 
where T: Float + FloatSimd<T> {
    dot(q, q).sqrt()
}

/// Normalization of vector3.
/// 
/// If you enter a zero vector, it returns a zero vector.
#[inline]
pub fn normalize_vec<T>(v: Vector3<T>) -> Vector3<T>
where T: Float {
    let coef = norm_vec(v).recip();
    if coef.is_infinite() {
        ZERO_VECTOR()
    } else {
        scale_vec(coef, v)
    }
}

/// Normalization of quaternion.
#[inline]
pub fn normalize<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float + FloatSimd<T> {
    scale( norm(q).recip(), q )
}

/// Invert the sign of a vector.
/// 
/// return: `-v`
#[inline]
pub fn negate_vec<T>(v: Vector3<T>) -> Vector3<T>
where T: Float {
    [ -v[0], -v[1], -v[2] ]
}

/// Invert the sign of a quaternion.
/// 
/// return: `-q`
#[inline]
pub fn negate<T>(q: Quaternion<T>) -> Quaternion<T>
where T: FloatSimd<T> {
    T::negate(q)
}

/// Product of pure quaternions.
/// 
/// `ab ≡ -a・b + a×b` (!= ba)
#[inline]
pub fn mul_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Quaternion<T>
where T: Float {
    ( -dot_vec(a, b), cross_vec(a, b) )
}

/// Hamilton product.
/// 
/// The product order is `ab (!= ba)`
#[inline]
pub fn mul<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>
where T: Float + FloatSimd<T> {
    let a0_b = scale(a.0, b);
    (
        a0_b.0 - dot_vec(a.1, b.1),
        add_vec( scale_add_vec(b.0, a.1, a0_b.1), cross_vec(a.1, b.1) )
    )
}

/// Compute the conjugate quaternion.
#[inline]
pub fn conj<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( q.0, negate_vec(q.1) )
}

/// Compute the inverse pure quaternion.
#[inline]
pub fn inv_vec<T>(v: Vector3<T>) -> Vector3<T>
where T: Float {
    scale_vec( dot_vec(v, v).recip(), negate_vec(v) )
}

/// Compute the inverse quaternion.
#[inline]
pub fn inv<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float + FloatSimd<T> {
    scale( dot(q, q).recip(), conj(q) )
}

/// Exponential function of vector3.
#[inline]
pub fn exp_vec<T>(v: Vector3<T>) -> Quaternion<T>
where T: Float {
    let norm_v = norm_vec(v);
    let (sin, cos) = norm_v.sin_cos();
    ( cos, scale_vec(sin / norm_v, v) )
}

/// Exponential function of quaternion.
#[inline]
pub fn exp<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    let norm_q_v = norm_vec(q.1);
    let (sin, cos) = norm_q_v.sin_cos();
    let coef = q.0.exp();
    ( coef * cos, scale_vec((coef * sin) / norm_q_v, q.1) )
}

/// Natural logarithm of quaternion.
#[inline]
pub fn ln<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    let tmp = dot_vec(q.1, q.1);
    let norm_q = (q.0*q.0 + tmp).sqrt();
    let coef = (q.0 / norm_q).acos() / tmp.sqrt();
    ( norm_q.ln(), scale_vec(coef, q.1) )
}

/// Natural logarithm of versor.
/// 
/// If it is guaranteed to be a versor, it is less computationally 
/// expensive than the `ln` function.
/// 
/// Only the vector part is returned since the real part is always zero.
#[inline]
pub fn ln_versor<T>(q: Quaternion<T>) -> Vector3<T>
where T: Float {
    scale_vec( acos_safe(q.0) / norm_vec(q.1), q.1)
}

/// Power function of quaternion.
#[inline]
pub fn pow<T>(q: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float {
    let tmp = dot_vec(q.1, q.1);
    let norm_q = (q.0*q.0 + tmp).sqrt();
    let omega = (q.0 / norm_q).acos();
    let (sin, cos) = (t * omega).sin_cos();
    let coef = norm_q.powf(t);
    ( coef * cos, scale_vec((coef * sin) / tmp.sqrt(), q.1) )
}

/// Power function of versor.
/// 
/// If it is guaranteed to be a versor, it is less computationally 
/// expensive than the `pow` function. 
#[inline]
pub fn pow_versor<T>(q: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float {
    let (sin, cos) = (t * acos_safe(q.0)).sin_cos();
    ( cos, scale_vec(sin / norm_vec(q.1), q.1) )
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
/// # let pi = std::f64::consts::PI;
/// // Make these as you like.
/// let v = [1.0, 0.5, -8.0];
/// let q = from_axis_angle([0.2, 1.0, -2.0], pi);
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
    let tmp = scale_add_vec(q.0, v, cross_vec(q.1, v));
    scale_add_vec(cast(2.0), cross_vec(q.1, tmp), v)
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
/// # let pi = std::f64::consts::PI;
/// // Make these as you like.
/// let v = [1.0, 0.5, -8.0];
/// let q = from_axis_angle([0.2, 1.0, -2.0], pi);
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
    let tmp = scale_add_vec(q.0, v, cross_vec(v, q.1));
    scale_add_vec(cast(2.0), cross_vec(tmp, q.1), v)
}

/// Calculate a versor to rotate from vector `a` to `b`.
/// 
/// If you enter a zero vector, it returns an identity quaternion.
/// 
/// # Examples
/// 
/// ```
/// # use quaternion_core::{Vector3, cross_vec, rotate_a_to_b, point_rotation};
/// 
/// let a: Vector3<f64> = [1.5, -0.5, 0.2];
/// let b: Vector3<f64> = [0.1, 0.6, 1.0];
/// 
/// let q = rotate_a_to_b(a, b);
/// let b_check = point_rotation(q, a);
/// 
/// let cross = cross_vec(b, b_check);
/// assert!( cross[0].abs() < 1e-12 );
/// assert!( cross[1].abs() < 1e-12 );
/// assert!( cross[2].abs() < 1e-12 );
/// ```
#[inline]
pub fn rotate_a_to_b<T>(a: Vector3<T>, b: Vector3<T>) -> Quaternion<T>
where T: Float {
    let half: T = cast(0.5);

    let t = dot_vec(a, b);
    let s_square = dot_vec(a, a) * dot_vec(b, b);
    let e_half = half * (t / s_square.sqrt());
    let v = ((half - e_half) / (s_square - t * t)).sqrt();

    // vがfiniteならeもfiniteである．
    if v.is_finite() {
        ( (half + e_half).sqrt(), scale_vec(v, cross_vec(a, b)) )
    } else {
        IDENTITY()
    }
}

/// Calculate a versor to rotate from vector `a` to `b`.
/// 
/// The parameter `t` adjusts the amount of movement from `a` to `b`, 
/// so that When `t=1`, it moves to position `b` completely.
/// 
/// If you enter a zero vector, it returns an identity quaternion.
/// 
/// If `t=1` at all times, it is less computationally expensive to use `rotate_a_to_b` function.
#[inline]
pub fn rotate_a_to_b_param<T>(a: Vector3<T>, b: Vector3<T>, t: T) -> Quaternion<T>
where T: Float {
    let dot_ab = dot_vec(a, b);
    let norm_ab_square = dot_vec(a, a) * dot_vec(b, b);
    let tmp_acos = dot_ab / norm_ab_square.sqrt();
    if tmp_acos.is_infinite() {
        IDENTITY()
    } else {
        let theta = acos_safe(tmp_acos);
        let (sin, cos) = ( t * theta * cast(0.5) ).sin_cos();
        let coef_v = sin / (norm_ab_square - dot_ab * dot_ab).sqrt();
        if coef_v.is_finite() {
            ( cos, scale_vec(coef_v, cross_vec(a, b)) )
        } else {
            IDENTITY()
        }
    }
}

/// Lerp (Linear interpolation)
/// 
/// Generate a quaternion that interpolate the shortest path from `a` to `b` 
/// (The norm of `a` and `b` must be 1).
/// The argument `t (0 <= t <= 1)` is the interpolation parameter.
/// 
/// Normalization is not performed internally because 
/// it increases the computational complexity.
#[inline]
pub fn lerp<T>(a: Quaternion<T>, b: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float + FloatSimd<T> {
    debug_assert!(
        t >= T::zero() && t <= T::one(), 
        "Parameter `t` must be in the range [0, 1]."
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
/// Generate a quaternion that interpolate the shortest path from `a` to `b`.
/// The argument `t(0 <= t <= 1)` is the interpolation parameter.
/// 
/// The norm of `a` and `b` must be 1 (Versor).
#[inline]
pub fn slerp<T>(a: Quaternion<T>, mut b: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float + FloatSimd<T> {
    debug_assert!(
        t >= T::zero() && t <= T::one(), 
        "Parameter `t` must be in the range [0, 1]."
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
    // たまにacosが抜けると計算時間を把握しにくくなるから，この実装とする．
    ( x.abs().min( T::one() ) * x.signum() ).acos()
}