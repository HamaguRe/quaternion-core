pub type Vector3 = [f64; 3];
pub type Quaternion = (f64, Vector3);


/// 恒等四元数を生成．
/// Generate identity quaternion.
#[inline(always)]
pub fn id() -> Quaternion {
    (1.0, [0.0; 3])
}

/// 回転角と軸ベクトルを指定して四元数を生成．
/// axisは単位ベクトルでなくても良い．
/// Generate quaternion by specifying rotation angle and axis vector.
/// The "axes" need not be unit vectors.
/// angle[rad]
#[inline(always)]
pub fn axis_angle(axis: Vector3, angle: f64) -> Quaternion {
    let n = normalize_vec(axis);
    let q_s = (angle / 2.0).cos();
    let s = (angle / 2.0).sin();
    let q_v = mul_scalar_vec(s, n);
    (q_s, q_v)
}

/// ベクトルの内積
/// Dot product of vector
#[inline(always)]
pub fn dot_vec(a: Vector3, b: Vector3) -> f64 {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

/// クォータニオンの内積
/// Dot product of quaternion
#[inline(always)]
pub fn dot(a: Quaternion, b: Quaternion) -> f64 {
    a.0 * b.0 + dot_vec(a.1, b.1)
}

/// 外積
/// Cross product
#[inline(always)]
pub fn cross_vec(a: Vector3, b: Vector3) -> Vector3 {
    let vec0 = a[1]*b[2] - a[2]*b[1];
    let vec1 = a[2]*b[0] - a[0]*b[2];
    let vec2 = a[0]*b[1] - a[1]*b[0];
    [vec0, vec1, vec2]
}

/// 二つのベクトルを加算する
/// Add two vectors
#[inline(always)]
pub fn add_vec(a: Vector3, b: Vector3) -> Vector3 {
    [a[0]+b[0], a[1]+b[1], a[2]+b[2]]
}

/// 二つのクォータニオンを加算する
/// Add two quaternions.
#[inline(always)]
pub fn add(a: Quaternion, b: Quaternion) -> Quaternion {
    ( a.0 + b.0, add_vec(a.1, b.1) )
}

/// クォータニオン積
/// 積の順序は "ab"(!=ba)
/// Multiplication of quaternion.
/// The product order is "ab"(!= ba)
#[inline(always)]
pub fn mul(a: Quaternion, b: Quaternion) -> Quaternion {
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
pub fn mul_scalar_vec(s: f64, v: Vector3) -> Vector3 {
    [s*v[0], s*v[1], s*v[2]]
}

/// スカラーとクォータニオンの積
/// Multiplication of scalar and quaternion.
#[inline(always)]
pub fn mul_scalar_quat(s: f64, a: Quaternion) -> Quaternion {
    ( s * a.0, mul_scalar_vec(s, a.1) )
}

/// L2ノルムを計算
/// Calculate L2 norm
#[inline(always)]
pub fn norm(a: Quaternion) -> f64 {
    dot(a, a).sqrt()
}

/// L2ノルムを計算
/// Calculate L2 norm
#[inline(always)]
pub fn norm_vec(r: Vector3) -> f64 {
    dot_vec(r, r).sqrt()
}

/// ノルムが1になるように正規化
/// Normalized so that norm is 1
#[inline(always)]
pub fn normalize(a: Quaternion) -> Quaternion {
    mul_scalar_quat(1.0 / norm(a), a)
}

/// 正規化
/// Normalization
#[inline(always)]
pub fn normalize_vec(r: Vector3) -> Vector3 {
    let norm = norm_vec(r);
    if norm == 0.0 {
        return [0.0; 3];
    }
    mul_scalar_vec(1.0 / norm, r)
}

/// 共役クォータニオンを求める
/// Compute conjugated quaternion
#[inline(always)]
pub fn conj(a: Quaternion) -> Quaternion {
    ( a.0, [-a.1[0], -a.1[1], -a.1[2]] )
}

/// 逆クォータニオンを求める
/// Compute inverse quaternion
#[inline(always)]
pub fn inverse(a: Quaternion) -> Quaternion {
    let conj = conj(a);
    let norm_square = dot(a, a);
    mul_scalar_quat(1.0 / norm_square, conj)
}

/// ネイピア数eのクォータニオン冪
/// Exponential of Quaternion.
#[inline(always)]
pub fn exp(a: Quaternion) -> Quaternion {
    let coef = a.0.exp();  // coefficient（係数）
    let vec_norm = norm_vec(a.1);
    // An if statement to avoid unnecessary calculation.
    if vec_norm == 0.0 {
        return (coef, [0.0; 3]);
    }
    let q_s = vec_norm.cos();
    let n   = normalize_vec(a.1);
    let q_v = mul_scalar_vec(vec_norm.sin(), n);
    mul_scalar_quat( coef, (q_s, q_v) )
}

/// クォータニオンの冪乗
/// The power of quaternion.
#[inline(always)]
pub fn power(a: Quaternion, t: f64) -> Quaternion {
    let coef = norm(a).powf(t);
    let omega = a.0.acos();
    let n = normalize_vec(a.1);
    let q_s = (t * omega).cos();
    let q_v = mul_scalar_vec( (t * omega).sin(), n );
    mul_scalar_quat( coef, (q_s, q_v) )
}

/// クォータニオンの自然対数
/// Natural logarithm of quaternion
#[inline(always)]
pub fn ln(a: Quaternion) -> Quaternion {
    let vec_norm = norm_vec(a.1);
    let q_s = vec_norm.ln();
    let s = (a.0 / vec_norm).acos();
    let n = normalize_vec(a.1);
    let q_v = mul_scalar_vec(s, n);
    (q_s, q_v)
}

/// 位置ベクトルの回転
/// Coordinate rotate. (r <-- a r a*)
#[inline(always)]
pub fn vector_rotation(a: Quaternion, r: Vector3) -> Vector3 {
    let a = normalize(a);
    let a_conj = conj(a);
    // ベクトルを，スカラー部0のクォータニオンとして計算する．
    let result = mul( a, mul((0.0, r), a_conj) );
    result.1
}

/// 座標系の回転
/// Vector rotate. (r <-- a* r a)
#[inline(always)]
pub fn coordinate_rotation(a: Quaternion, r: Vector3) -> Vector3 {
    let a = normalize(a);
    let a_conj = conj(a);
    let result = mul( a_conj, mul((0.0, r), a) );
    result.1
}

/// ベクトル "a" から "b" への回転を行うクォータニオンを求める．
/// Find a quaternion to rotate from vector "a" to "b".
/// t (0 <= t <= 1)
#[inline(always)]
pub fn rotation_a_to_b(a: Vector3, b: Vector3, t: f64) -> Quaternion {
    let axis = cross_vec(a, b);
    let s = norm_vec(a) * norm_vec(b);
    let theta = (dot_vec(a, b) / s).acos();
    axis_angle(axis, theta * t)
}

/// The integrate of angular velocity.
/// 角速度を積分して，引数に渡したクォータニオン"q"を更新する．
/// Update the quaternion "q" passed to the argument.
/// 
/// omega[rad/sec]
/// dt[sec]
#[inline(always)]
pub fn integration(omega: Vector3, q: Quaternion, dt: f64) -> Quaternion {
    let arg = mul_scalar_vec(dt / 2.0, omega);
    let dq = exp( (0.0, arg) );
    normalize( mul(dq, q) )
}

/// The integrate of angular velocity.
/// 角速度を積分して，引数に渡したクォータニオン"q"を更新する．
/// 近似式を用いるため，"integration()"関数よりも計算量が少ない．
/// "dt"を大きくしすぎると真値との誤差が大きくなる．
/// Update the quaternion "q" passed to the argument.
/// Since it uses an approximate expression,
/// calculation amount is smaller than "integration()" function.
/// If "dt" is made too large, the error becomes large.
/// 
/// omega[rad/sec]
/// dt[sec]
#[inline(always)]
pub fn integration_1(omega: Vector3, q: Quaternion, dt: f64) -> Quaternion {
    let dq = mul((0.0, omega), q);
    let dq = mul_scalar_quat(dt / 2.0, dq);
    normalize( add(q, dq) )
}

/// 線形補間
/// "a"から"b"への経路を補完するクォータニオンを生成する．
/// 引数t(0 <= t <= 1) は補間パラメータ．
/// Linear interpolation
/// Generate a quaternion that interpolate the route from "a" to "b".
/// The argument t(0 <= t <= 1) is the interpolation parameter.
#[inline(always)]
pub fn lerp(a: Quaternion, b: Quaternion, t: f64) -> Quaternion {
    let a = normalize(a);
    let b = normalize(b);
    let q_1 = mul_scalar_quat(1.0 - t, a);
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
pub fn slerp(a: Quaternion, b: Quaternion, t: f64) -> Quaternion {
    // Normalize to avoid undefined behavior.
    let a = normalize(a);
    let mut b = normalize(b);

    // 最短経路で補間する．
    let mut dot = dot(a, b);
    if dot < 0.0 {
        b = mul_scalar_quat(-1.0, b);
        dot = -dot;
    }

    // If the inputs are too close for comfort, linearly interpolate.
    const DOT_THRESHOLD: f64 = 0.9995;
    if dot > DOT_THRESHOLD {
        return lerp(a, b, t);
    }

    // selrp
    let omega = dot.acos();  // Angle between the two quaternion
    let sin_omega = omega.sin();
    let s_1 = ((1.0 - t)*omega).sin() / sin_omega;
    let q_1 = mul_scalar_quat(s_1, a);
    let s_2 = (t * omega).sin() / sin_omega;
    let q_2 = mul_scalar_quat(s_2, b);
    add(q_1, q_2)
}

/// Sherical linear interpolation. 
/// Use quaternion's exponential.
/// "a" --> "b".
#[inline(always)]
pub fn slerp_1(a: Quaternion, b: Quaternion, t: f64) -> Quaternion {
    let a = normalize(a);
    let mut b = normalize(b);

    let dot = dot(a, b);
    if dot < 0.0 {
        b = mul_scalar_quat(-1.0, b);
    }

    // a.0 が1より大きくなった時の適切な処理がわからない...
    // Set it to 1 to avoid undefined behavior.
    // The domain of the scalar part is [-1 <= a.0 <= 1].
    let mut arg = mul( conj(a), b );
    if arg.0.abs() > 1.0 {
        arg.0 = arg.0.round();
    }
    let tmp = power( arg, t);
    mul(a, tmp)
}
