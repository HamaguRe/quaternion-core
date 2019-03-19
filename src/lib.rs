// Quaternion Libraly
// f32 or f64

extern crate num_traits;
use num_traits::float::{Float, FloatConst};

pub type Vector3<T> = [T; 3];
pub type Quaternion<T> = (T, Vector3<T>);


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
    let q_v = mul_scalar_vec(s, n);
    (q_s, q_v)
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

/// クォータニオンをオイラー角[rad]に変換
/// Return --> [roll, pitch, yaw]
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
    let m13 = two * (q0*q2 - q1*q3);
    let m23 = two * (q2*q3 + q0*q1);
    let m33 = two * (q0*q0 + q3*q3) - one;

    let roll = (m23 / m33).atan();
    let pitch = asin_safe(m13);
    let yaw = (m12 / m11).atan();

    [roll, pitch, yaw]
}

/// ベクトルの内積
/// Dot product of vector
#[inline(always)]
fn dot_vec<T>(a: Vector3<T>, b: Vector3<T>) -> T 
where T: Float {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

/// クォータニオンの内積
/// Dot product of quaternion
#[inline(always)]
pub fn dot<T>(a: Quaternion<T>, b: Quaternion<T>) -> T 
where T: Float {
    a.0 * b.0 + dot_vec(a.1, b.1)
}

/// 外積
/// Cross product
#[inline(always)]
fn cross_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T> 
where T: Float {
    let vec0 = a[1]*b[2] - a[2]*b[1];
    let vec1 = a[2]*b[0] - a[0]*b[2];
    let vec2 = a[0]*b[1] - a[1]*b[0];
    [vec0, vec1, vec2]
}

/// 二つのベクトルを加算する
/// Add two vectors
#[inline(always)]
fn add_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T> 
where T: Float {
    [a[0]+b[0], a[1]+b[1], a[2]+b[2]]
}

/// 二つのクォータニオンを加算する
/// Add two quaternions.
#[inline(always)]
pub fn add<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    ( a.0 + b.0, add_vec(a.1, b.1) )
}

/// クォータニオン積
/// 積の順序は "ab"(!=ba)
/// Multiplication of quaternion.
/// The product order is "ab"(!= ba)
#[inline(always)]
pub fn mul<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    let q_s = a.0 * b.0 - dot_vec(a.1, b.1);
    let vec_1 = mul_scalar_vec(a.0, b.1);
    let vec_2 = mul_scalar_vec(b.0, a.1);
    let vec_3 = cross_vec(a.1, b.1);
    let q_v = add_vec( vec_1, add_vec(vec_2, vec_3) );
    (q_s, q_v)
}

/// スカラーとベクトルの積
/// Multiplication of scalar and vector.
#[inline(always)]
pub fn mul_scalar_vec<T>(s: T, v: Vector3<T>) -> Vector3<T> 
where T: Float {
    [s*v[0], s*v[1], s*v[2]]
}

/// スカラーとクォータニオンの積
/// Multiplication of scalar and quaternion.
#[inline(always)]
pub fn mul_scalar_quat<T>(s: T, a: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    ( s * a.0, mul_scalar_vec(s, a.1) )
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
    mul_scalar_quat(one / norm(a), a)
}

/// 正規化
/// Normalization
#[inline(always)]
pub fn normalize_vec<T>(r: Vector3<T>) -> Vector3<T> 
where T: Float {
    let zero = T::zero();
    let one  = T::one();

    let norm = norm_vec(r);
    if norm == zero {  // ゼロ除算回避
        return [zero; 3];
    }
    mul_scalar_vec(one / norm, r)
}

/// 共役クォータニオンを求める
/// Compute conjugated quaternion
#[inline(always)]
pub fn conj<T>(a: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    ( a.0, [-a.1[0], -a.1[1], -a.1[2]] )
}

/// 逆クォータニオンを求める
/// Compute inverse quaternion
#[inline(always)]
pub fn inverse<T>(a: Quaternion<T>) -> Quaternion<T> 
where T: Float {
    let one = T::one();
    let norm_square = dot(a, a);
    mul_scalar_quat( one / norm_square, conj(a) )
}

/// ネイピア数eのクォータニオン冪
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
    let q_v = mul_scalar_vec(vec_norm.sin(), n);
    mul_scalar_quat( coef, (q_s, q_v) )
}

/// クォータニオンの冪乗
/// The power of quaternion.
#[inline(always)]
pub fn power<T>(a: Quaternion<T>, t: T) -> Quaternion<T> 
where T: Float + FloatConst {
    let coef = norm(a).powf(t);
    let omega = acos_safe(a.0);
    let n = normalize_vec(a.1);
    let q_s = (t * omega).cos();
    let q_v = mul_scalar_vec( (t * omega).sin(), n );
    mul_scalar_quat( coef, (q_s, q_v) )
}

/// クォータニオンの自然対数
/// Natural logarithm of quaternion
#[inline(always)]
pub fn ln<T>(a: Quaternion<T>) -> Quaternion<T> 
where T: Float + FloatConst {
    let vec_norm = norm_vec(a.1);
    let q_s = vec_norm.ln();
    let s = acos_safe(a.0 / vec_norm);
    let n = normalize_vec(a.1);
    let q_v = mul_scalar_vec(s, n);
    (q_s, q_v)
}

/// 位置ベクトルの回転
/// Coordinate rotate. (r' = a r a*)
#[inline(always)]
pub fn vector_rotation<T>(a: Quaternion<T>, r: Vector3<T>) -> Vector3<T> 
where T: Float {
    let zero = T::zero();
    let a = normalize(a);
    // ベクトルを，スカラー部0のクォータニオンとして計算する．
    let result = mul( mul(a, (zero, r)), conj(a) );
    result.1
}

/// 座標系の回転
/// Vector rotate. (r' = a* r a)
#[inline(always)]
pub fn coordinate_rotation<T>(a: Quaternion<T>, r: Vector3<T>) -> Vector3<T> 
where T: Float {
    let zero = T::zero();
    let a = normalize(a);
    let result = mul( conj(a), mul((zero, r), a) );
    result.1
}

/// ベクトル "a" から "b" への回転を行うクォータニオンを求める．
/// Find a quaternion to rotate from vector "a" to "b".
#[inline(always)]
pub fn rotation_a_to_b<T>(a: Vector3<T>, b: Vector3<T>) -> Quaternion<T> 
where T: Float + FloatConst {
    let axis = cross_vec(a, b);
    let s = norm_vec(a) * norm_vec(b);
    let theta = acos_safe( dot_vec(a, b) / s );
    from_axis_angle(axis, theta)
}

/// The integrate of angular velocity.
/// dt間のクォータニオンの変化量を返す．
/// omega[rad/sec]
/// dt[sec]
#[inline(always)]
pub fn integration<T>(omega: Vector3<T>, dt: T) -> Quaternion<T> 
where T: Float {
    let zero = T::zero();
    let two = T::one() + T::one();
    let arg = mul_scalar_vec(dt / two, omega);
    exp( (zero, arg) )
}

/// ボディ角速度を積分して，引数に渡したクォータニオンを更新する．
/// Update the quaternion "q" passed to the argument.
#[inline(always)]
pub fn vector_integration<T>(q: Quaternion<T>, omega: Vector3<T>, dt: T) -> Quaternion<T> 
where T: Float {
    let dq = integration(omega, dt);
    mul(dq, q)
}

/// 空間角速度を積分して，引数に渡したクォータニオンを積分する．
#[inline(always)]
pub fn coordinate_integration<T>(q: Quaternion<T>, omega: Vector3<T>, dt: T) -> Quaternion<T> 
where T: Float {
    let dq = integration(omega, dt);
    mul(q, dq)
}

/// オイラー法
/// 機体角速度を積分して，引数に渡したクォータニオンを積分する．
/// omega[rad/sec]
/// dt[sec]
#[inline(always)]
pub fn vector_integration_euler<T>(q: Quaternion<T>, omega: Vector3<T>, dt: T) -> Quaternion<T> 
where T: Float {
    let zero = T::zero();
    let two = T::one() + T::one();

    let dq = mul( (zero, omega), q);  // 機体角速度を用いる
    let dq = mul_scalar_quat(dt / two, dq);
    normalize( add(q, dq) )
}

/// オイラー法
/// 空間角速度を積分して，引数に渡したクォータニオンを更新する．
/// omega[rad/sec]
/// dt[sec]
#[inline(always)]
pub fn coordinate_integration_euler<T>(q: Quaternion<T>, omega: Vector3<T>, dt: T) -> Quaternion<T> 
where T: Float {
    let zero = T::zero();
    let two = T::one() + T::one();

    let dq = mul( q, (zero, omega) );  // 空間角速度を用いる
    let dq = mul_scalar_quat(dt / two, dq);
    normalize( add(q, dq) )
}

/// 線形補間
/// 引数"a"から"b"への経路を補完するクォータニオンを生成する．
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
    let q_1 = mul_scalar_quat(one - t, a);
    let q_2 = mul_scalar_quat(t, b);
    let result = add(q_1, q_2);
    normalize(result)
}

/// 球状線形補間
/// "a"から"b"への経路を補完するクォータニオンを生成する．
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
        b = mul_scalar_quat(-one, b);
        dot = -dot;
    }
    // If the inputs are too close for comfort, linearly interpolate.
    if dot > threshold {
        return lerp(a, b, t);
    }
    // selrp
    let omega = acos_safe(dot);  // Angle between the two quaternion
    let sin_omega = omega.sin();
    let s_1 = ((one - t)*omega).sin() / sin_omega;
    let q_1 = mul_scalar_quat(s_1, a);
    let s_2 = (t * omega).sin() / sin_omega;
    let q_2 = mul_scalar_quat(s_2, b);
    add(q_1, q_2)
}

/// クォータニオンの冪乗を用いたSlerp
/// Sherical linear interpolation. 
/// Use quaternion's exponential.
/// "a" --> "b".
#[inline(always)]
pub fn slerp_1<T>(a: Quaternion<T>, b: Quaternion<T>, t: T) -> Quaternion<T> 
where T: Float + FloatConst {
    let zero = T::zero();
    let one  = T::one();
    let threshold: T = num_traits::cast(0.9995).unwrap();

    let a = normalize(a);
    let mut b = normalize(b);
    // 最短経路で補完
    let mut dot = dot(a, b);
    if dot < zero {
        b = mul_scalar_quat(-one, b);
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

/// 定義域外の値はカットして未定義動作を防ぐ
fn asin_safe<T>(s: T) -> T 
where T: Float + FloatConst {
    let one = T::one();
    let two = one + one;
    let pi = T::PI();

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
fn acos_safe<T>(s: T) -> T 
where T: Float + FloatConst {
    let zero = T::zero();
    let one = T::one();
    let pi = T::PI();

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