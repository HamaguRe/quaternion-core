//! Quaternion Libraly (f32 & f64)

#![no_std]
#[cfg(feature = "std")]
extern crate std;

use num_traits::float::{Float, FloatConst};

/// [i, j, k]
pub type Vector3<T> = [T; 3];

/// (1, [i, j, k])
pub type Quaternion<T> = (T, Vector3<T>);

/// Direction Cosine Matrix
pub type DCM<T> = [Vector3<T>; 3];


/// 回転角\[rad\]と軸ベクトルを指定してVersorを生成する．
/// 
/// axisは単位ベクトルでなくても良い．  
/// 零ベクトルを入力した場合は，恒等四元数を返す．
/// 
/// Generate Versor by specifying rotation angle\[rad\] and axis vector.  
/// The "axis" need not be unit vector.  
/// If you enter a zero vector, it returns an identity quaternion.
#[inline]
pub fn from_axis_angle<T>(axis: Vector3<T>, angle: T) -> Quaternion<T>
where T: Float + FloatConst {
    let norm_vec = norm_vec(axis);
    if norm_vec < T::epsilon() {  // ゼロ除算回避
        IDENTITY()
    } else {
        let tmp = angle % ( T::PI() + T::PI() );  // limit to (-2π, 2π)
        let f = ( tmp * cast(0.5) ).sin_cos();
        ( f.1, scale_vec(f.0 / norm_vec, axis) )
    }
}

/// Versorから回転軸（単位ベクトル）と軸回りの回転角\[rad\]を求める．
/// 
/// Compute the rotation axis (unit vector) and the rotation angle\[rad\] 
/// around the axis from the versor.
/// 
/// return: `(axis, angle)`
/// 
/// Range of angle: `-π <= angle <= π`
#[inline]
pub fn to_axis_angle<T>(q: Quaternion<T>) -> (Vector3<T>, T)
where T: Float + FloatConst {
    let norm_vec = norm_vec(q.1);
    if norm_vec < T::epsilon() {
        ( ZERO_VECTOR(), T::zero() )
    } else {
        let axis = scale_vec( norm_vec.recip(), q.1 );
        let tmp = asin_safe(norm_vec);
        let angle = copysign(tmp + tmp, q.0);  // 2*tmp
        (axis, angle)
    }
}

/// 位置ベクトル回転（`q v q*`）を表す方向余弦行列をVersorに変換する．
/// 
/// 座標系回転（`q* v q`）を表す方向余弦行列を変換する場合には，
/// ```ignore
/// let q = conj( from_dcm(dcm) );
/// ```
/// とする．
/// 
/// Generate versor from direction cosine matrix, 
/// representing rotation of position vector.
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
            q1 = half * tmp;
            q0 = (m[2][1] - m[1][2]) * coef;
            q2 = (m[0][1] + m[1][0]) * coef;
            q3 = (m[0][2] + m[2][0]) * coef;
        },
        2 => {
            q2 = half * tmp;
            q0 = (m[0][2] - m[2][0]) * coef;
            q1 = (m[0][1] + m[1][0]) * coef;
            q3 = (m[1][2] + m[2][1]) * coef;
        },
        3 => {
            q3 = half * tmp;
            q0 = (m[1][0] - m[0][1]) * coef;
            q1 = (m[0][2] + m[2][0]) * coef;
            q2 = (m[1][2] + m[2][1]) * coef;
        },
        _ => unreachable!(),
    };

    (q0, [q1, q2, q3])
}

/// 位置ベクトル回転（`q v q*`）を表すVersorを，方向余弦行列に変換する．
/// 
/// 座標系回転（`q* v q`）を表すVersorを変換する場合には，
/// ```ignore
/// let dcm = to_dcm( conj(q) );
/// ```
/// とする．
#[inline]
pub fn to_dcm<T>(q: Quaternion<T>) -> DCM<T>
where T: Float {
    let one = T::one();
    let two = cast(2.0);

    // Compute these value only once.
    let q0_q0 = q.0 * q.0;
    let q0_q1 = q.0 * q.1[0];
    let q0_q2 = q.0 * q.1[1];
    let q0_q3 = q.0 * q.1[2];
    let q1_q2 = q.1[0] * q.1[1];
    let q1_q3 = q.1[0] * q.1[2];
    let q2_q3 = q.1[1] * q.1[2];

    let m11 = mul_add(mul_add(q.1[0], q.1[0], q0_q0), two, -one);
    let m12 = (q1_q2 - q0_q3) * two;
    let m13 = (q1_q3 + q0_q2) * two;
    let m21 = (q1_q2 + q0_q3) * two;
    let m22 = mul_add(mul_add(q.1[1], q.1[1], q0_q0), two, -one);
    let m23 = (q2_q3 - q0_q1) * two;
    let m31 = (q1_q3 - q0_q2) * two;
    let m32 = (q2_q3 + q0_q1) * two;
    let m33 = mul_add(mul_add(q.1[2], q.1[2], q0_q0), two, -one);

    [
        [m11, m12, m13],
        [m21, m22, m23],
        [m31, m32, m33],
    ]
}

/// `z-y-x`系のオイラー角[rad]を四元数に変換する．
/// 
/// ypr: [yaw, pitch, roll]
/// 
/// Aerospace angle/axis sequences is `[yaw, pitch, roll] / [z, y, x]`.  
#[inline]
pub fn from_euler_angles<T>(ypr: Vector3<T>) -> Quaternion<T> 
where T: Float {
    let [alpha, beta, gamma] = scale_vec(cast(0.5), ypr);

    // Compute these value only once
    let (sin_alpha, cos_alpha) = alpha.sin_cos();
    let (sin_beta,  cos_beta)  = beta.sin_cos();
    let (sin_gamma, cos_gamma) = gamma.sin_cos();
    let cos_alpha_cos_beta = cos_alpha * cos_beta;
    let sin_alpha_sin_beta = sin_alpha * sin_beta;

    let q0 = cos_alpha_cos_beta   * cos_gamma + sin_alpha_sin_beta   * sin_gamma;
    let q1 = cos_alpha_cos_beta   * sin_gamma - sin_alpha_sin_beta   * cos_gamma;
    let q2 = cos_alpha * sin_beta * cos_gamma + sin_alpha * cos_beta * sin_gamma;
    let q3 = sin_alpha * cos_beta * cos_gamma - cos_alpha * sin_beta * sin_gamma;

    (q0, [q1, q2, q3])
}

/// 位置ベクトル回転（`q v q*`）を表すVersorを`z-y-x系`のオイラー角\[rad\]に変換する．
/// 
/// 座標系回転（`q* v q`）を表すVersorを変換する場合には，
/// ```ignore
/// let euler_angles = to_euler_angle( conj(q) );
/// ```
/// とする．
/// 
/// Aerospace angle/axis sequences is `[yaw, pitch, roll] / [z, y, x]`.  
/// pitch angle is limited to [-π/2, π/2]．
/// 
/// return: `[yaw, pitch, roll]`
#[inline]
pub fn to_euler_angle<T>(q: Quaternion<T>) -> Vector3<T>
where T: Float + FloatConst {
    let [
        [m11,   _, m13],
        [m21,   _, m23],
        [m31, m32, m33],
    ] = to_dcm(q);

    if m13.abs() < (T::one() - T::epsilon()) {
        [
            (m32 / m33).atan(),  // yaw
            (-m31).asin(),       // pitch
            (m21 / m11).atan(),  // roll
        ]
    } else {  // ジンバルロック
        [
            T::zero(),                       // yaw
            copysign(T::FRAC_PI_2(), -m13),  // pitch
            (m23 / m13).atan(),              // roll
        ]
    }
}

/// 回転ベクトル(rotation vector)をVersorに変換
/// 
/// 回転ベクトルはそれ自体が回転軸を表し，ノルムは軸回りの回転角(0, 2π)を表す．
#[inline]
pub fn from_rotation_vector<T>(v: Vector3<T>) -> Quaternion<T>
where T: Float + FloatConst {
    let theta = norm_vec(v);
    if theta < T::epsilon() {
        IDENTITY()
    } else {
        let f = ( theta * cast(0.5) ).sin_cos();
        ( f.1, scale_vec(f.0 / theta, v) )
    }
}

/// Versorを回転ベクトル(rotation vector)に変換
/// 
/// 回転ベクトルはそれ自体が回転軸を表し，ノルムは軸回りの回転角(0, 2π)を表す．
#[inline]
pub fn to_rotation_vector<T>(q: Quaternion<T>) -> Vector3<T>
where T: Float + FloatConst {
    let norm_vec = norm_vec(q.1);
    if norm_vec < T::epsilon() {
        ZERO_VECTOR()
    } else {
        let tmp = acos_safe(q.0);
        scale_vec((tmp + tmp) / norm_vec, q.1)  // 2*tmp
    }
}

/// 方向余弦行列を用いてベクトルを回転させる．
#[inline]
pub fn matrix_product<T>(m: DCM<T>, v: Vector3<T>) -> Vector3<T>
where T: Float {
    [
        dot_vec(m[0], v),
        dot_vec(m[1], v),
        dot_vec(m[2], v),
    ]
}

/// 右手系と左手系の四元数を変換
/// 
/// x, z軸の向きはそのままで，y軸と全ての軸回りの回転方向を反転
#[inline]
pub fn system_trans<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    // 実部と，ベクトル部の要素どれか一つの符号を反転させれば良い
    ( -q.0, [ q.1[0], -q.1[1], q.1[2] ] )
}

/// Dot product of vector.
#[inline]
pub fn dot_vec<T>(a: Vector3<T>, b: Vector3<T>) -> T
where T: Float {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

/// Dot product of quaternion.
#[inline]
pub fn dot<T>(a: Quaternion<T>, b: Quaternion<T>) -> T 
where T: Float {
    a.0 * b.0 + dot_vec(a.1, b.1)
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

/// Calcurate `a + b`
#[inline]
pub fn add_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T>
where T: Float {
    [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ]
}

/// Calcurate `a + b`
#[inline]
pub fn add<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( a.0 + b.0, add_vec(a.1, b.1) )
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
where T: Float {
    ( a.0 - b.0, sub_vec(a.1, b.1) )
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
where T: Float {
    ( s * q.0, scale_vec(s, q.1) )
}

/// Calculate `s*a + b`
/// 
/// `fma` featureを有効にした場合は`mul_add`メソッドを用いてFMA計算を行い，
/// 有効にしなかった場合は単純な積和（s*a + b）に展開して計算する．
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
/// `fma` featureを有効にした場合は`mul_add`メソッドを用いてFMA計算を行い，
/// 有効にしなかった場合は単純な積和（s*a + b）に展開して計算する．
#[inline]
pub fn scale_add<T>(s: T, a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( mul_add(s, a.0, b.0), scale_add_vec(s, a.1, b.1) )
}

/// Hadamard product of Vector.
#[inline]
pub fn hadamard_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T>
where T: Float {
    [ a[0]*b[0], a[1]*b[1], a[2]*b[2] ]
}

/// Hadamard product of Quaternion.
#[inline]
pub fn hadamard<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( a.0 * b.0, hadamard_vec(a.1, b.1) )
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
where T: Float {
    dot(q, q).sqrt()
}

/// Normalization of vector.
/// 
/// 零ベクトルを入力した場合は零ベクトルを返す．
/// 
/// If you enter a zero vector, it returns a zero vector.
#[inline]
pub fn normalize_vec<T>(v: Vector3<T>) -> Vector3<T>
where T: Float + FloatConst {
    let norm_vec = norm_vec(v);
    if norm_vec < T::epsilon() {
        ZERO_VECTOR()
    } else {
        scale_vec( norm_vec.recip(), v )
    }
}

/// Normalization of quaternion.
#[inline]
pub fn normalize<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    scale( norm(q).recip(), q )
}

/// ベクトルの符号を反転．
/// 
/// return: `-v`
#[inline]
pub fn negate_vec<T>(v: Vector3<T>) -> Vector3<T>
where T: Float {
    [ -v[0], -v[1], -v[2] ]
}

/// 四元数の符号を反転．
/// 
/// return: `-q`
#[inline]
pub fn negate<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( -q.0, negate_vec(q.1) )
}

/// 純虚四元数同士の積．
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
where T: Float {
    let v0 = scale_vec(a.0, b.1);
    let v1 = scale_vec(b.0, a.1);
    let s = a.0 * b.0 - dot_vec(a.1, b.1);
    let v = add_vec( add_vec(v0, v1), cross_vec(a.1, b.1) );
    (s, v)
}

/// 共役四元数を求める．
/// 
/// Compute the conjugate quaternion.
#[inline]
pub fn conj<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( q.0, negate_vec(q.1) )
}

/// 逆純虚四元数を求める．
/// 
/// 零ベクトルを入力した場合は零ベクトルを返す．
/// 
/// Compute the inverse pure quaternion.
#[inline]
pub fn inv_vec<T>(v: Vector3<T>) -> Vector3<T>
where T: Float + FloatConst {
    let dot_vec = dot_vec(v, v);
    if dot_vec < T::epsilon() {
        ZERO_VECTOR()
    } else {
        scale_vec( dot_vec.recip(), negate_vec(v) )
    }
}

/// 逆四元数を求める．
/// 
/// Compute the inverse quaternion.
#[inline]
pub fn inv<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    scale( dot(q, q).recip(), conj(q) )
}

/// Exponential function of vector.
#[inline]
pub fn exp_vec<T>(v: Vector3<T>) -> Quaternion<T>
where T: Float + FloatConst {
    let norm_vec = norm_vec(v);
    let (sin, cos) = norm_vec.sin_cos();
    ( cos, scale_vec(sin / norm_vec, v) )
}

/// Exponential function of quaternion.
#[inline]
pub fn exp<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float + FloatConst {
    let norm_vec = norm_vec(q.1);
    let (sin, cos) = norm_vec.sin_cos();
    let coef = q.0.exp();
    ( coef * cos, scale_vec((coef * sin) / norm_vec, q.1) )
}

/// Natural logarithm of quaternion.
#[inline]
pub fn ln<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float + FloatConst {
    let tmp = dot_vec(q.1, q.1);
    let norm = (q.0*q.0 + tmp).sqrt();
    let coef = acos_safe(q.0 / norm) / tmp.sqrt();
    ( norm.ln(), scale_vec(coef, q.1) )
}

/// Natural logarithm of versor.
/// 
/// Versorであることが保証されている場合にはln関数よりも計算量を減らせる．
/// 実部は必ず0になるので省略．
#[inline]
pub fn ln_versor<T>(q: Quaternion<T>) -> Vector3<T>
where T: Float + FloatConst {
    scale_vec(q.0.acos() / norm_vec(q.1), q.1)
}

/// Power function of quaternion.
#[inline]
pub fn pow<T>(q: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float + FloatConst {
    let tmp = dot_vec(q.1, q.1);
    let norm = (q.0*q.0 + tmp).sqrt();
    let omega = acos_safe(q.0 / norm);
    let (sin, cos) = (t * omega).sin_cos();
    let coef = norm.powf(t);
    ( coef * cos, scale_vec((coef * sin) / tmp.sqrt(), q.1) )
}

/// Power function of versor.
/// 
/// Versorであることが保証されている場合にはpow関数よりも計算量を減らせる．
#[inline]
pub fn pow_versor<T>(q: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float + FloatConst {
    let (sin, cos) = ( t * q.0.acos() ).sin_cos();
    ( cos, scale_vec(sin / norm_vec(q.1), q.1) )
}

/// 位置ベクトルの回転
/// 
/// `q v q*  (||q|| = 1)`
#[inline]
pub fn vector_rotation<T>(q: Quaternion<T>, v: Vector3<T>) -> Vector3<T>
where T: Float {
    let tmp = scale_add_vec(q.0, v, cross_vec(q.1, v));
    scale_add_vec(cast(2.0), cross_vec(q.1, tmp), v)
}

/// 座標系の回転
/// 
/// `q* v q  (||q|| = 1)`
#[inline]
pub fn frame_rotation<T>(q: Quaternion<T>, v: Vector3<T>) -> Vector3<T>
where T: Float {
    let tmp = scale_add_vec(q.0, v, cross_vec(v, q.1));
    scale_add_vec(cast(2.0), cross_vec(tmp, q.1), v)
}

/// 位置ベクトル`a`を 位置ベクトル`b`と同じ場所へ最短距離で回転させるVersorを求める．
/// 
/// 零ベクトルを入力した場合は，恒等四元数を返す．
/// 
/// Find a versor to rotate from vector `a` to `b`.
/// If you enter a zero vector, it returns an identity quaternion.
/// 
/// Range of rotation parameter: `0 <= t <= 1`
#[inline]
pub fn rotate_a_to_b<T>(a: Vector3<T>, b: Vector3<T>, t: T) -> Quaternion<T>
where T: Float + FloatConst {
    debug_assert!(
        T::zero() <= t && t <= T::one(),
        "rotate_a_to_b(): The rotate parameter `t` is out of range."
    );

    let dot_ab = dot_vec(a, b);
    let norm_ab_square = dot_vec(a, a) * dot_vec(b, b);
    let tmp = norm_ab_square - dot_ab * dot_ab;

    if tmp < T::epsilon() {
        IDENTITY()
    } else {
        let angle = acos_safe( dot_ab / norm_ab_square.sqrt() );
        let (sin, cos) = ( t * angle * cast(0.5) ).sin_cos();
        ( cos, scale_vec(sin / tmp.sqrt(), cross_vec(a, b)) )
    }
}

/// Lerp (Linear interpolation)
/// 
/// Generate a quaternion that interpolate the route from `a` to `b`.
/// The argument `t (0 <= t <= 1)` is the interpolation parameter.
/// 
/// The norm of `a` and `b` must be 1 (Versor).
#[inline]
pub fn lerp<T>(a: Quaternion<T>, b: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float {
    debug_assert!(
        T::zero() <= t && t <= T::one(),
        "lerp(): The interpolation parameter `t` is out of range."
    );

    normalize( scale_add(t, sub(b, a), a) )
}

/// Slerp (Spherical linear interpolation)
/// 
/// Generate a quaternion that interpolate the route from `a` to `b`.
/// The argument `t(0 <= t <= 1)` is the interpolation parameter.
/// 
/// The norm of `a` and `b` must be 1 (Versor).
#[inline]
pub fn slerp<T>(a: Quaternion<T>, mut b: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float {
    debug_assert!(
        T::zero() <= t && t <= T::one(),
        "slerp(): The interpolation parameter `t` is out of range."
    );

    // 最短経路で補間
    let mut dot = dot(a, b);
    if dot.is_sign_negative() {
        b = negate(b);
        dot = -dot;
    }
    // If the distance between quaternions is close enough, use lerp.
    if dot > cast(0.9995) {  // Approximation error < 0.017%
        lerp(a, b, t)
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

/// 定数呼び出し以外に使わないのでエラー処理を省略
#[inline(always)]
fn cast<T: Float>(x: f64) -> T {
    num_traits::cast::<f64, T>(x).unwrap()
}

/// xの絶対値を持ち，かつsignの符号を持つ値を返す．
#[inline(always)]
fn copysign<T: Float>(x: T, sign: T) -> T {
    if sign.is_sign_positive() {
         x.abs()
    } else {
        -x.abs()
    }
}

/// `fma` featureを有効にした場合は`mul_add`メソッドとして展開され，
/// 有効にしなかった場合は単純な積和（`s*a + b`）に展開してコンパイルされる．
#[inline(always)]
fn mul_add<T: Float>(s: T, a: T, b: T) -> T {
    if cfg!(feature = "fma") {
        s.mul_add(a, b)
    } else {
        s * a + b
    }
}

/// 配列内の最大値を探して，その位置と値を返す．
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

/// 定義域外の値をカットして未定義動作を防ぐ
#[inline(always)]
fn asin_safe<T: Float + FloatConst>(x: T) -> T {
    if x.abs() < T::one() {
        x.asin()
    } else {
        copysign(T::FRAC_PI_2(), x)
    }
}

/// 定義域外の値をカットして未定義動作を防ぐ
#[inline(always)]
fn acos_safe<T: Float + FloatConst>(x: T) -> T {
    if x.abs() < T::one() {
        x.acos()
    } else {
        if x.is_sign_positive() { T::zero() } else { T::PI() }
    }
}