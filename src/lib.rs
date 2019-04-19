// Quaternion Libraly
// f32 or f64

extern crate num_traits;
use num_traits::float::{Float, FloatConst};

pub type Vector3<T> = [T; 3];
pub type Quaternion<T> = (T, Vector3<T>);
pub type DirectionCosines<T> = [Vector3<T>; 3];


/// 恒等四元数を生成．
/// Generate identity quaternion.
#[inline(always)]
pub fn id<T>() -> Quaternion<T>
where T: Float {
    let one  = T::one();
    let zero = T::zero();
    
    (one, [zero; 3])
}

/// 回転角と軸ベクトルを指定して四元数を生成．
/// axisは単位ベクトルでなくても良い．
/// Generate quaternion by specifying rotation angle and axis vector.
/// The "axes" need not be unit vectors.
/// angle[rad]
#[inline(always)]
pub fn from_axis_angle<T>(axis: Vector3<T>, angle: T) -> Quaternion<T> 
where T: Float {
    let two = T::one() + T::one();

    let n = normalize_vec(axis);
    let q_s = (angle / two).cos();
    let s   = (angle / two).sin();
    let q_v = scale_vec(s, n);
    (q_s, q_v)
}

/// 方向余弦行列から四元数を生成．
#[inline(always)]
pub fn from_direction_cosines<T>(m: DirectionCosines<T>) -> Quaternion<T>
where T: Float {
    let one = T::one();
    let two = one + one;
    let four = two + two;

    let q0 = (one / two) * (m[0][0] + m[1][1] + m[2][2] + one).sqrt();
    let q1 = (m[1][2] - m[2][1]) / (four * q0);
    let q2 = (m[2][0] - m[0][2]) / (four * q0);
    let q3 = (m[0][1] - m[1][0]) / (four * q0);
    (q0, [q1, q2, q3])
}

/// オイラー角[rad]から四元数を生成．
#[inline(always)]
pub fn from_euler_angles<T>(roll: T, pitch: T, yaw: T) -> Quaternion<T> 
where T: Float {
    let two = T::one() + T::one();

    let alpha = yaw   / two;
    let beta  = pitch / two;
    let gamma = roll  / two;
    // Compute these value only once
    let sin_alpha = alpha.sin();
    let cos_alpha = alpha.cos();
    let sin_beta  = beta.sin();
    let cos_beta  = beta.cos();
    let sin_gamma = gamma.sin();
    let cos_gamma = gamma.cos();

    let q0 = cos_alpha * cos_beta * cos_gamma + sin_alpha * sin_beta * sin_gamma;
    let q1 = cos_alpha * cos_beta * sin_gamma - sin_alpha * sin_beta * cos_gamma;
    let q2 = cos_alpha * sin_beta * cos_gamma + sin_alpha * cos_beta * sin_gamma;
    let q3 = sin_alpha * cos_beta * cos_gamma - cos_alpha * sin_beta * sin_gamma;
    (q0, [q1, q2, q3])
}

/// 四元数を方向余弦行列（回転行列）に変換
#[inline(always)]
pub fn to_direction_cosines<T>(q: Quaternion<T>) -> DirectionCosines<T>
where T: Float {
    let one = T::one();
    let two = one + one;

    let q0 = q.0;
    let q1 = q.1[0];
    let q2 = q.1[1];
    let q3 = q.1[2];

    let m11 = two * (q0*q0 + q1*q1) - one;
    let m12 = two * (q1*q2 + q0*q3);
    let m13 = two * (q1*q3 - q0*q2);
    let m21 = two * (q1*q2 - q0*q3);
    let m22 = two * (q0*q0 + q2*q2) - one;
    let m23 = two * (q2*q3 + q0*q1);
    let m31 = two * (q1*q3 + q0*q2);
    let m32 = two * (q2*q3 - q0*q1);
    let m33 = two * (q0*q0 + q3*q3) - one;

    [
        [m11, m12, m13],
        [m21, m22, m23],
        [m31, m32, m33]
    ]
}

/// 四元数をオイラー角[rad]に変換
/// Quaternion --> [roll, pitch, yaw]
#[inline(always)]
pub fn to_euler_angles<T>(q: Quaternion<T>) -> Vector3<T> 
where T: Float + FloatConst {
    let one = T::one();
    let two = one + one;

    let q0 = q.0;
    let q1 = q.1[0];
    let q2 = q.1[1];
    let q3 = q.1[2];
    
    let m11 = two * (q0*q0 + q1*q1) - one;
    let m12 = two * (q1*q2 + q0*q3);
    let m13 = two * (q1*q3 - q0*q2);
    let m23 = two * (q2*q3 + q0*q1);
    let m33 = two * (q0*q0 + q3*q3) - one;

    let roll  = (m23 / m33).atan();
    let pitch = asin_safe(-m13);
    let yaw   = (m12 / m11).atan();

    [roll, pitch, yaw]
}

/// ベクトルの内積
/// Dot product of vector
#[inline(always)]
pub fn dot_vec<T>(a: Vector3<T>, b: Vector3<T>) -> T 
where T: Float {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

/// 四元数の内積
/// Dot product of quaternion
#[inline(always)]
pub fn dot<T>(a: Quaternion<T>, b: Quaternion<T>) -> T 
where T: Float {
    a.0 * b.0 + dot_vec(a.1, b.1)
}

/// 外積
/// Cross product
#[inline(always)]
pub fn cross_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T> 
where T: Float {
    let vec0 = a[1]*b[2] - a[2]*b[1];
    let vec1 = a[2]*b[0] - a[0]*b[2];
    let vec2 = a[0]*b[1] - a[1]*b[0];
    [vec0, vec1, vec2]
}

/// 二つのベクトルを加算する
/// Add two vectors
#[inline(always)]
pub fn add_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T> 
where T: Float {
    [a[0]+b[0], a[1]+b[1], a[2]+b[2]]
}

/// 二つの四元数を加算する
/// Add two quaternions.
#[inline(always)]
pub fn add<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    ( a.0 + b.0, add_vec(a.1, b.1) )
}

/// ハミルトン積のために定義した，特別なベクトル同士の積
/// ab ≡ -a・b + a×b
#[inline(always)]
pub fn mul_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Quaternion<T> 
where T: Float {
    ( -dot_vec(a, b), cross_vec(a, b) )
}

/// ハミルトン積
/// 積の順序は "ab"(!=ba)
/// Hamilton product
/// The product order is "ab"(!= ba)
#[inline(always)]
pub fn mul<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    let tmp0 = mul_vec(a.1, b.1);
    let vec0 = scale_vec(a.0, b.1);
    let vec1 = scale_vec(b.0, a.1);
    let tmp1 = ( a.0 * b.0, add_vec(vec0, vec1) );
    add(tmp0, tmp1)
}

/// スカラーとベクトルの積
/// Multiplication of scalar and vector.
#[inline(always)]
pub fn scale_vec<T>(s: T, v: Vector3<T>) -> Vector3<T> 
where T: Float {
    [s*v[0], s*v[1], s*v[2]]
}

/// スカラーと四元数の積
/// Multiplication of scalar and quaternion.
#[inline(always)]
pub fn scale<T>(s: T, a: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    ( s * a.0, scale_vec(s, a.1) )
}

/// L2ノルムを計算
/// Calculate L2 norm
#[inline(always)]
pub fn norm<T>(a: Quaternion<T>) -> T 
where T: Float {
    dot(a, a).sqrt()
}

/// L2ノルムを計算
/// Calculate L2 norm
#[inline(always)]
pub fn norm_vec<T>(r: Vector3<T>) -> T 
where T: Float {
    dot_vec(r, r).sqrt()
}

/// ノルムが1になるように正規化
/// Normalized so that norm is 1
#[inline(always)]
pub fn normalize<T>(a: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    let one = T::one();

    scale(one / norm(a), a)
}

/// 正規化
/// Normalization
#[inline(always)]
pub fn normalize_vec<T>(r: Vector3<T>) -> Vector3<T> 
where T: Float {
    let zero = T::zero();
    let one  = T::one();

    let norm = norm_vec(r);
    if norm == zero {
        return [zero; 3];  // ゼロ除算回避
    }
    scale_vec(one / norm, r)
}

/// 符号反転
#[inline(always)]
pub fn sign_inversion_vec<T>(r: Vector3<T>) -> Vector3<T> 
where T: Float {
    [ -r[0], -r[1], -r[2] ]
}

/// 符号反転
#[inline(always)]
pub fn sign_inversion<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( -q.0, sign_inversion_vec(q.1) )
}

/// 共役四元数を求める
/// Compute conjugated quaternion
#[inline(always)]
pub fn conj<T>(a: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    ( a.0, sign_inversion_vec(a.1) )
}

/// 逆四元数を求める
/// Compute inverse quaternion
#[inline(always)]
pub fn inverse<T>(a: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    let one = T::one();

    let norm_square = dot(a, a);
    scale( one / norm_square, conj(a) )
}

/// ネイピア数eの四元数冪
/// Exponential of Quaternion<T>.
#[inline(always)]
pub fn exp<T>(a: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    let zero = T::zero();

    let coef = a.0.exp();  // coefficient（係数）
    let vec_norm = norm_vec(a.1);

    // An if statement to avoid unnecessary calculation.
    if vec_norm == zero {
        return (coef, [zero; 3]);
    }
    
    let q_s = vec_norm.cos();
    let n   = normalize_vec(a.1);
    let q_v = scale_vec(vec_norm.sin(), n);
    scale( coef, (q_s, q_v) )
}

/// 四元数の冪乗
/// The power of quaternion.
#[inline(always)]
pub fn power<T>(a: Quaternion<T>, t: T) -> Quaternion<T> 
where T: Float + FloatConst {
    let coef = norm(a).powf(t);
    let omega = acos_safe(a.0);
    let n = normalize_vec(a.1);
    let q_s = (t * omega).cos();
    let q_v = scale_vec( (t * omega).sin(), n );
    scale( coef, (q_s, q_v) )
}

/// 四元数の自然対数
/// Natural logarithm of quaternion
#[inline(always)]
pub fn ln<T>(a: Quaternion<T>) -> Quaternion<T> 
where T: Float + FloatConst {
    let vec_norm = norm_vec(a.1);
    let q_s = vec_norm.ln();
    let s = acos_safe(a.0 / vec_norm);
    let n = normalize_vec(a.1);
    let q_v = scale_vec(s, n);
    (q_s, q_v)
}

/// 位置ベクトルの回転
/// r' = q r q*
#[inline(always)]
pub fn vector_rotation<T>(q: Quaternion<T>, r: Vector3<T>) -> Vector3<T> 
where T: Float {
    let two = T::one() + T::one();

    let q = normalize(q);
    let term1 = scale_vec( q.0 * q.0, r );
    let term2 = scale_vec( two * q.0, cross_vec(q.1, r) );
    let term3 = scale_vec( dot_vec(r, q.1), q.1 );
    let term4 = cross_vec( q.1, cross_vec(q.1, r) );
    add_vec( add_vec(term1, term2), add_vec(term3, term4) )
}

/// 座標系の回転
/// r' = q* r q
#[inline(always)]
pub fn coordinate_rotation<T>(q: Quaternion<T>, r: Vector3<T>) -> Vector3<T> 
where T: Float {
    let two = T::one() + T::one();

    let q = normalize(q);
    let term1 = scale_vec( q.0 * q.0, r );
    let term2 = scale_vec( two * q.0, cross_vec(r, q.1) );
    let term3 = scale_vec( dot_vec(r, q.1), q.1 );
    let term4 = cross_vec( q.1, cross_vec(q.1, r) );
    add_vec( add_vec(term1, term2), add_vec(term3, term4) )
}

/// ベクトル "a" を ベクトル "b" へ最短距離で回転させる四元数を求める．
/// Find a quaternion to rotate from vector "a" to "b".
/// 0 <= t <= 1
#[inline(always)]
pub fn rotation_a_to_b<T>(a: Vector3<T>, b: Vector3<T>, t: T) -> Quaternion<T> 
where T: Float + FloatConst {
    let axis = cross_vec(a, b);
    let s = norm_vec(a) * norm_vec(b);
    let theta = acos_safe( dot_vec(a, b) / s );
    from_axis_angle(axis, theta * t)
}

/// The integrate of angular velocity.
/// 機体角速度を積分して，引数に渡した四元数を更新する．
/// omega[rad/sec]
/// dt[sec]
#[inline(always)]
pub fn integration<T>(q: Quaternion<T>, omega: Vector3<T>, dt: T) -> Quaternion<T> 
where T: Float {
    let zero = T::zero();
    let two  = T::one() + T::one();
    
    let arg = scale_vec(dt / two, omega);
    let dq  = exp( (zero, arg) );
    mul(q, dq)
}

/// オイラー法
/// The integration of body angular velocity.
/// 機体角速度を積分して，引数に渡した四元数を更新する．
/// omega[rad/sec]
/// dt[sec]
#[inline(always)]
pub fn integration_euler<T>(q: Quaternion<T>, omega: Vector3<T>, dt: T) -> Quaternion<T> 
where T: Float {
    let zero = T::zero();
    let two  = T::one() + T::one();

    let dq = mul( q, (zero, omega) );  // 空間角速度を用いる
    let dq = scale(dt / two, dq);
    normalize( add(q, dq) )
}

/// 線形補間
/// 引数"a"から"b"への経路を補完する四元数を生成する．
/// 引数t(0 <= t <= 1) は補間パラメータ．
/// Linear interpolation
/// Generate a quaternion that interpolate the route from "a" to "b".
/// The argument t(0 <= t <= 1) is the interpolation parameter.
#[inline(always)]
pub fn lerp<T>(a: Quaternion<T>, b: Quaternion<T>, t: T) -> Quaternion<T> 
where T: Float {
    let one = T::one();

    let a = normalize(a);
    let b = normalize(b);
    let q_1 = scale(one - t, a);
    let q_2 = scale(t, b);
    normalize( add(q_1, q_2) )
}

/// 球状線形補間
/// "a"から"b"への経路を補完する四元数を生成する．
/// 引数t(0 <= t <= 1) は補間パラメータ．
/// Spherical linear interpolation
/// Generate a quaternion that interpolate the route from "a" to "b".
/// The argument t(0 <= t <= 1) is the interpolation parameter.
#[inline(always)]
pub fn slerp<T>(a: Quaternion<T>, b: Quaternion<T>, t: T) -> Quaternion<T> 
where T: Float + FloatConst {
    let zero = T::zero();
    let one  = T::one();
    let threshold: T = num_traits::cast(0.9995).unwrap();

    // Normalize to avoid undefined behavior.
    let a = normalize(a);
    let mut b = normalize(b);
    // 最短経路で補間する．
    let mut dot = dot(a, b);
    if dot < zero {
        b = sign_inversion(b);
        dot = -dot;
    }
    // If the inputs are too close for comfort, linearly interpolate.
    if dot > threshold {
        return lerp(a, b, t);
    }
    // selrp
    let omega = acos_safe(dot);  // Angle between the two quaternion
    let sin_omega = omega.sin();
    let s_1 = ( (one - t) * omega ).sin() / sin_omega;
    let q_1 = scale(s_1, a);
    let s_2 = (t * omega).sin() / sin_omega;
    let q_2 = scale(s_2, b);
    add(q_1, q_2)
}

/// 四元数の冪乗を用いたSlerp
/// Sherical linear interpolation. 
/// Use quaternion's exponential.
/// "a" --> "b".
#[inline(always)]
pub fn slerp_1<T>(a: Quaternion<T>, b: Quaternion<T>, t: T) -> Quaternion<T> 
where T: Float + FloatConst {
    let zero = T::zero();
    let threshold: T = num_traits::cast(0.9995).unwrap();

    let a = normalize(a);
    let mut b = normalize(b);
    // 最短経路で補完
    let mut dot = dot(a, b);
    if dot < zero {
        b = sign_inversion(b);
        dot = -dot;
    }
    // lerp
    if dot > threshold {
        return lerp(a, b, t);
    }
    // slerp
    let tmp = mul( conj(a), b );
    let tmp = power(tmp, t);
    mul(a, tmp)
}

/// 定義域外の値をカットして未定義動作を防ぐ
#[inline(always)]
fn asin_safe<T>(s: T) -> T 
where T: Float + FloatConst {
    let one = T::one();
    let two = one + one;
    let pi  = T::PI();

    let result;
    if s >= one {  // Avoid undefined behavior
        result = pi / two;
    } else if s <= -one {
        result = -pi / two;
    } else {
        result = s.asin();
    }
    result
}

/// 定義域外の値をカットして未定義動作を防ぐ
#[inline(always)]
fn acos_safe<T>(s: T) -> T 
where T: Float + FloatConst {
    let zero = T::zero();
    let one  = T::one();
    let pi   = T::PI();

    let result;
    if s >= one {  // Avoid undefined behavior
        result = zero;
    } else if s <= -one{
        result = pi;
    } else {
        result = s.acos();
    }
    result
}