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
/// The "axis" need not be unit vector.
/// angle[rad]
#[inline(always)]
pub fn from_axis_angle<T>(axis: Vector3<T>, angle: T) -> Quaternion<T> 
where T: Float {
    let half: T = num_traits::cast(0.5).unwrap();

    let n = normalize_vec(axis);
    let q_s = (angle * half).cos();
    let s   = (angle * half).sin();
    let q_v = scale_vec(s, n);
    (q_s, q_v)
}

/// 方向余弦行列から四元数を生成．
#[inline(always)]
pub fn from_direction_cosines<T>(m: DirectionCosines<T>) -> Quaternion<T>
where T: Float {
    let one = T::one();
    let half: T = num_traits::cast(0.5).unwrap();
    let four: T = num_traits::cast(4.0).unwrap();

    let q0 = (m[0][0] + m[1][1] + m[2][2] + one).sqrt() * half;
    let tmp = (four * q0).recip();  // reciprocal
    let q1 = (m[1][2] - m[2][1]) * tmp;
    let q2 = (m[2][0] - m[0][2]) * tmp;
    let q3 = (m[0][1] - m[1][0]) * tmp;
    (q0, [q1, q2, q3])
}

/// オイラー角[rad]から四元数を生成．
#[inline(always)]
pub fn from_euler_angles<T>(roll: T, pitch: T, yaw: T) -> Quaternion<T> 
where T: Float {
    let half: T = num_traits::cast(0.5).unwrap();

    let alpha = yaw   * half;
    let beta  = pitch * half;
    let gamma = roll  * half;
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
pub fn to_direction_cosines<T>((q0, [q1, q2, q3]): Quaternion<T>) -> DirectionCosines<T>
where T: Float {
    let one = T::one();
    let two = T::one() + T::one();

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
pub fn to_euler_angles<T>((q0, [q1, q2, q3]): Quaternion<T>) -> Vector3<T> 
where T: Float + FloatConst {
    let one = T::one();
    let two = T::one() + T::one();

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

/// 回転を表す四元数から，単位ベクトル（回転軸）を取り出す．
/// normalize_vec関数よりも計算量が少ないが，
/// norm(q)==1であることを前提とする．
/// 引数に渡す四元数のノルムが保証できない場合には
/// normalize_vec関数を用いるべき．
#[inline(always)]
pub fn get_unit_vector<T>(q: Quaternion<T>) -> Vector3<T>
where T: Float {
    let zero = T::zero();
    let one  = T::one();

    if q.0 == one {
        return [zero; 3];  // ゼロ除算回避
    }
    scale_vec( (one - q.0 * q.0).sqrt().recip(), q.1)
}

/// 回転を表す四元数から，軸周りの回転角[rad]を取り出す．
#[inline(always)]
pub fn get_angle<T>(q: Quaternion<T>) -> T
where T: Float + FloatConst {
    let two = T::one() + T::one();

    let q = normalize(q);
    two * acos_safe(q.0)
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
/// a・b
#[inline(always)]
pub fn dot<T>(a: Quaternion<T>, b: Quaternion<T>) -> T 
where T: Float {
    a.0 * b.0 + dot_vec(a.1, b.1)
}

/// 外積
/// Cross product
/// a×b
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
    [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ]
}

/// 二つの四元数を加算する
/// Add two quaternions.
#[inline(always)]
pub fn add<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    ( a.0 + b.0, add_vec(a.1, b.1) )
}

// ベクトルの減算
// Vector subtraction
// Calculate "a - b"
#[inline(always)]
pub fn sub_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T>
where T: Float {
    [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ]
}

// 四元数の減算
// Quaternion subtraction
// Calculate "a - b"
#[inline(always)]
pub fn sub<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( a.0 - b.0, sub_vec(a.1, b.1) )
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
    [ s*v[0], s*v[1], s*v[2] ]
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
pub fn norm_vec<T>(r: Vector3<T>) -> T 
where T: Float {
    dot_vec(r, r).sqrt()
}

/// L2ノルムを計算
/// Calculate L2 norm
#[inline(always)]
pub fn norm<T>(a: Quaternion<T>) -> T 
where T: Float {
    dot(a, a).sqrt()
}

/// 正規化
/// Normalization
#[inline(always)]
pub fn normalize_vec<T>(r: Vector3<T>) -> Vector3<T> 
where T: Float {
    let zero = T::zero();

    let norm = norm_vec(r);
    if norm == zero {
        return [zero; 3];  // ゼロ除算回避
    }
    scale_vec( norm.recip(), r )
}

/// ノルムが1になるように正規化
/// Normalized so that norm is 1
#[inline(always)]
pub fn normalize<T>(a: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    scale( norm(a).recip(), a )
}

/// 符号反転
/// return "-r"
#[inline(always)]
pub fn negate_vec<T>(r: Vector3<T>) -> Vector3<T> 
where T: Float {
    [ -r[0], -r[1], -r[2] ]
}

/// 符号反転
/// return "-q"
#[inline(always)]
pub fn negate<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( -q.0, negate_vec(q.1) )
}

/// 共役四元数を求める
/// Compute conjugated quaternion
#[inline(always)]
pub fn conj<T>(a: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    ( a.0, negate_vec(a.1) )
}

/// 逆四元数を求める
/// Compute inverse quaternion
#[inline(always)]
pub fn inverse<T>(a: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    scale( dot(a, a).recip(), conj(a) )
}

/// ネイピア数eのベクトル冪
/// Exponential of Vector3.
#[inline(always)]
pub fn exp_vec<T>(a: Vector3<T>) -> Quaternion<T> 
where T: Float {
    let norm = norm_vec(a);
    let n = normalize_vec(a);
    let q_v = scale_vec(norm.sin(), n);
    (norm.cos(), q_v)
}

/// ネイピア数eの四元数冪
/// Exponential of Quaternion.
#[inline(always)]
pub fn exp<T>(a: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    scale( a.0.exp(), exp_vec(a.1) )
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
    let norm = norm(a);
    let n = normalize_vec(a.1);
    let tmp = acos_safe(a.0 / norm);
    let q_v = scale_vec(tmp, n);
    ( norm.ln(), q_v )
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
    let norm_a_b = norm_vec(a) * norm_vec(b);
    let theta = acos_safe( dot_vec(a, b) / norm_a_b );
    from_axis_angle(axis, theta * t)
}

/// 機体座標系の姿勢変化を積分して，引数に渡した四元数を更新する．
/// Integrate attitude change of body coordinate system, 
/// and update Quaternion passed to the argument.
/// omega[rad/sec]
/// dt[sec]
#[inline(always)]
pub fn integration<T>(q: Quaternion<T>, omega: Vector3<T>, dt: T) -> Quaternion<T> 
where T: Float {
    let half: T = num_traits::cast(0.5).unwrap();
    
    let tmp = scale_vec(dt * half, omega);
    let dq  = exp_vec(tmp);
    mul(dq, q)
}

/// オイラー法
/// Euler method
/// 機体座標系の姿勢変化を積分して，引数に渡した四元数を更新する．
/// Integrate attitude change of body coordinate system, 
/// and update Quaternion passed to the argument.
/// omega[rad/sec]
/// dt[sec]
#[inline(always)]
pub fn integration_euler<T>(q: Quaternion<T>, omega: Vector3<T>, dt: T) -> Quaternion<T> 
where T: Float {
    let zero = T::zero();
    let half: T = num_traits::cast(0.5).unwrap();

    let tmp = mul( (zero, omega), q );
    let dq = scale(dt * half, tmp);
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
    let one  = T::one();
    let threshold: T = num_traits::cast(0.9995).unwrap();

    // Normalize to avoid undefined behavior.
    let a = normalize(a);
    let mut b = normalize(b);
    // 最短経路で補間する．
    let mut dot = dot(a, b);
    if dot.is_sign_negative() == true {
        b = negate(b);
        dot = -dot;
    }
    // If the inputs are too close for comfort, linearly interpolate.
    if dot > threshold {
        return lerp(a, b, t);
    }
    // selrp
    let omega = acos_safe(dot);  // Angle between the two quaternion
    let sin_omega_recip = omega.sin().recip();
    let s_1 = ( (one - t) * omega ).sin() * sin_omega_recip;
    let q_1 = scale(s_1, a);
    let s_2 = (t * omega).sin() * sin_omega_recip;
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
    let threshold: T = num_traits::cast(0.9995).unwrap();

    let a = normalize(a);
    let mut b = normalize(b);
    // 最短経路で補完
    let mut dot = dot(a, b);
    if dot.is_sign_negative() == true {
        b = negate(b);
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
    let frac_pi_2 = T::FRAC_PI_2(); // π/2

    if s >= one {  // Avoid undefined behavior
        return  frac_pi_2;
    } else if s <= -one {
        return -frac_pi_2;
    }
    s.asin()
}

/// 定義域外の値をカットして未定義動作を防ぐ
#[inline(always)]
fn acos_safe<T>(s: T) -> T 
where T: Float + FloatConst {
    let zero = T::zero();
    let one  = T::one();
    let pi   = T::PI();

    if s >= one {  // Avoid undefined behavior
        return zero;
    } else if s <= -one {
        return pi;
    }
    s.acos()
}