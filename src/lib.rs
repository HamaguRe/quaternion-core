// Quaternion Libraly
// f64 only
//
// get_unit_vector() は計算に特異点を持つ．

pub type Vector3<T> = [T; 3];
pub type Quaternion<T> = (T, Vector3<T>);
pub type DirectionCosines<T> = [Vector3<T>; 3];

const PI: f64 = std::f64::consts::PI;
const THRESHOLD: f64 = 0.9995;
const ZERO_VECTOR: Vector3<f64> = [0.0; 3];
pub const IDENTITY: Quaternion<f64> = (1.0, ZERO_VECTOR);  // Identity Quaternion


/// 回転角[rad]と軸ベクトルを指定して四元数を生成．
/// axisは単位ベクトルでなくても良い．
/// 零ベクトルを入力した場合は，恒等四元数を返す．
/// Generate quaternion by specifying rotation angle[rad] and axis vector.
/// The "axis" need not be unit vector.
/// If you enter a zero vector, it returns an identity quaternion.
#[inline(always)]
pub fn from_axis_angle(axis: Vector3<f64>, angle: f64) -> Quaternion<f64> {
    let norm = norm_vec(axis);
    if (norm == 0.0) | (angle == 0.0) {
        return IDENTITY;
    }
    let f = (angle * 0.5).sin_cos();
    ( f.1, scale_vec( f.0 / norm, axis ) )
}

/// 位置ベクトルの回転を表す方向余弦行列から四元数を生成．
#[inline(always)]
pub fn from_direction_cosines_vector(m: DirectionCosines<f64>) -> Quaternion<f64> {
    let q0 = (m[0][0] + m[1][1] + m[2][2] + 1.0).sqrt() * 0.5;
    let q1 = m[2][1] - m[1][2];
    let q2 = m[0][2] - m[2][0];
    let q3 = m[1][0] - m[0][1];
    let coef = (4.0 * q0).recip();
    normalize( ( q0, scale_vec(coef, [q1, q2, q3]) ) )
}

/// 座標系回転を表す方向余弦行列から四元数を生成．
/// Generate quaternion from direction cosine matrix.
#[inline(always)]
pub fn from_direction_cosines_frame(m: DirectionCosines<f64>) -> Quaternion<f64> {
    let q0 = (m[0][0] + m[1][1] + m[2][2] + 1.0).sqrt() * 0.5;
    let q1 = m[1][2] - m[2][1];
    let q2 = m[2][0] - m[0][2];
    let q3 = m[0][1] - m[1][0];
    let coef = (4.0 * q0).recip();  // reciprocal
    normalize( ( q0, scale_vec(coef, [q1, q2, q3]) ) )
}

/// 回転を表す四元数から，回転軸（単位ベクトル）と軸回りの回転角[rad]を取り出す．
/// Compute the rotation axis (unit vector) and the rotation angle[rad] 
/// around the axis from the quaternion representing the rotation.
/// return "(axis, angle)"
#[inline(always)]
pub fn to_axis_angle(q: Quaternion<f64>) -> (Vector3<f64>, f64) {
    let q = normalize(q);
    let axis = get_unit_vector(q);
    let angle = 2.0 * acos_safe(q.0);
    (axis, angle)
}

/// 位置ベクトルの回転を表す四元数を，方向余弦行列（回転行列）に変換．
/// q r q* と同じ回転を表す．
#[inline(always)]
pub fn to_direction_cosines_vector(q: Quaternion<f64>) -> DirectionCosines<f64> {
    let (q0, [q1, q2, q3]) = normalize(q);
    // Compute these value only once.
    let q0_q0 = q0 * q0;
    let q0_q1 = q0 * q1;
    let q0_q2 = q0 * q2;
    let q0_q3 = q0 * q3;
    let q1_q2 = q1 * q2;
    let q1_q3 = q1 * q3;
    let q2_q3 = q2 * q3;

    let m11 = (q0_q0 + q1*q1).mul_add(2.0, -1.0);
    let m12 = (q1_q2 - q0_q3) * 2.0;
    let m13 = (q1_q3 + q0_q2) * 2.0;
    let m21 = (q1_q2 + q0_q3) * 2.0;
    let m22 = (q0_q0 + q2*q2).mul_add(2.0, -1.0);
    let m23 = (q2_q3 - q0_q1) * 2.0;
    let m31 = (q1_q3 - q0_q2) * 2.0;
    let m32 = (q2_q3 + q0_q1) * 2.0;
    let m33 = (q0_q0 + q3*q3).mul_add(2.0, -1.0);

    [
        [m11, m12, m13],
        [m21, m22, m23],
        [m31, m32, m33]
    ]
}

/// 座標系の回転を表す四元数を，方向余弦行列（回転行列）に変換．
/// q* r q と同じ回転を表す．
#[inline(always)]
pub fn to_direction_cosines_frame(q: Quaternion<f64>) -> DirectionCosines<f64> {
    let (q0, [q1, q2, q3]) = normalize(q);
    // Compute these value only once.
    let q0_q0 = q0 * q0;
    let q0_q1 = q0 * q1;
    let q0_q2 = q0 * q2;
    let q0_q3 = q0 * q3;
    let q1_q2 = q1 * q2;
    let q1_q3 = q1 * q3;
    let q2_q3 = q2 * q3;

    let m11 = (q0_q0 + q1*q1).mul_add(2.0, -1.0);
    let m12 = (q1_q2 + q0_q3) * 2.0;
    let m13 = (q1_q3 - q0_q2) * 2.0;
    let m21 = (q1_q2 - q0_q3) * 2.0;
    let m22 = (q0_q0 + q2*q2).mul_add(2.0, -1.0);
    let m23 = (q2_q3 + q0_q1) * 2.0;
    let m31 = (q1_q3 + q0_q2) * 2.0;
    let m32 = (q2_q3 - q0_q1) * 2.0;
    let m33 = (q0_q0 + q3*q3).mul_add(2.0, -1.0);

    [
        [m11, m12, m13],
        [m21, m22, m23],
        [m31, m32, m33]
    ]
}

/// 方向余弦行列（回転行列）を用いてベクトルを回転させる．
/// 四元数を用いる場合よりも計算量が少ないため，多くのベクトルに対して回転を適用する場合には
/// 一度方向余弦行列に変換してからベクトルを回転させたほうが速度面で有利になるかもしれない．
#[inline(always)]
pub fn matrix_product(m: DirectionCosines<f64>, r: Vector3<f64>) -> Vector3<f64> {
    [
        dot_vec(m[0], r),
        dot_vec(m[1], r),
        dot_vec(m[2], r)
    ]
}

/// 回転を表す四元数から，単位ベクトル（回転軸）を取り出す．
/// normalize_vec関数よりも計算量が少ないが，
/// norm(q)==1であることを前提とする．
/// 引数に渡す四元数のノルムが保証できない場合には
/// normalize_vec関数を用いるべき．
/// 単位四元数を入力した場合には零ベクトルを返す．
#[inline(always)]
pub fn get_unit_vector(q: Quaternion<f64>) -> Vector3<f64> {
    let coef = ( q.0.mul_add(-q.0, 1.0) ).sqrt();
    if coef == 0.0 {
        return ZERO_VECTOR;  // ゼロ除算回避
    }
    scale_vec( coef.recip(), q.1 )
}

/// ベクトルのスカラー積（内積）
/// Dot product of vector
#[inline(always)]
pub fn dot_vec(a: Vector3<f64>, b: Vector3<f64>) -> f64 {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

/// 四元数のスカラー積（内積）
/// Dot product of quaternion
/// a・b
#[inline(always)]
pub fn dot(a: Quaternion<f64>, b: Quaternion<f64>) -> f64 {
    a.0 * b.0 + dot_vec(a.1, b.1)
}

/// ベクトル積（外積）
/// Cross product
/// a×b
#[inline(always)]
pub fn cross_vec(a: Vector3<f64>, b: Vector3<f64>) -> Vector3<f64> {
    [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    ]
}

/// 二つのベクトルを加算する
/// Add two vectors
#[inline(always)]
pub fn add_vec(a: Vector3<f64>, b: Vector3<f64>) -> Vector3<f64> {
    [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ]
}

/// 二つの四元数を加算する
/// Add two quaternions.
#[inline(always)]
pub fn add(a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
    ( a.0 + b.0, add_vec(a.1, b.1) )
}

/// ベクトルの減算
/// Vector subtraction
/// Calculate "a - b"
#[inline(always)]
pub fn sub_vec(a: Vector3<f64>, b: Vector3<f64>) -> Vector3<f64> {
    [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ]
}

/// 四元数の減算
/// Quaternion subtraction
/// Calculate "a - b"
#[inline(always)]
pub fn sub(a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
    ( a.0 - b.0, sub_vec(a.1, b.1) )
}

/// スカラーとベクトルの積
/// Multiplication of scalar and vector.
#[inline(always)]
pub fn scale_vec(s: f64, v: Vector3<f64>) -> Vector3<f64> {
    [ s*v[0], s*v[1], s*v[2] ]
}

/// スカラーと四元数の積
/// Multiplication of scalar and quaternion.
#[inline(always)]
pub fn scale(s: f64, a: Quaternion<f64>) -> Quaternion<f64> {
    ( s * a.0, scale_vec(s, a.1) )
}

/// compute "s*a + b"
/// CPUがFMA(Fused multiply–add)命令をサポートしている場合，
/// "add_vec(scale_vec(s, a), b)" よりも高速に計算出来る．
/// If the CPU supports FMA(Fused multiply–add) instructions,
/// this is faster than "add_vec(scale_vec(s, a), b)".
#[inline(always)]
pub fn scale_add_vec(s: f64, a: Vector3<f64>, b: Vector3<f64>) -> Vector3<f64> {
    [
        s.mul_add(a[0], b[0]),
        s.mul_add(a[1], b[1]),
        s.mul_add(a[2], b[2])
    ]
}

/// compute "s*a + b"
/// CPUがFMA(Fused multiply–add)命令をサポートしている場合，
/// "add(scale(s, a), b)" よりも高速に計算出来る．
/// If the CPU supports FMA(Fused multiply–add) instructions,
/// this is faster than "add(scale(s, a), b)".
#[inline(always)]
pub fn scale_add(s: f64, a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
    ( s.mul_add(a.0, b.0), scale_add_vec(s, a.1, b.1) )
}

/// L2ノルムを計算
/// Calculate L2 norm
#[inline(always)]
pub fn norm_vec(r: Vector3<f64>) -> f64 {
    dot_vec(r, r).sqrt()
}

/// L2ノルムを計算
/// Calculate L2 norm
#[inline(always)]
pub fn norm(a: Quaternion<f64>) -> f64 {
    dot(a, a).sqrt()
}

/// ノルムが1になるように正規化
/// 零ベクトルを入力した場合は零ベクトルを返す．
/// Normalization
/// If you enter a zero vector, it returns a zero vector.
#[inline(always)]
pub fn normalize_vec(r: Vector3<f64>) -> Vector3<f64> {
    let norm = norm_vec(r);
    if norm == 0.0 {
        return ZERO_VECTOR;  // ゼロ除算回避
    }
    scale_vec( norm.recip(), r )
}

/// ノルムが1になるように正規化
/// Normalized so that norm is 1
#[inline(always)]
pub fn normalize(a: Quaternion<f64>) -> Quaternion<f64> {
    scale( norm(a).recip(), a )
}

/// 符号反転
/// return "-r"
#[inline(always)]
pub fn negate_vec(r: Vector3<f64>) -> Vector3<f64> {
    [ -r[0], -r[1], -r[2] ]
}

/// 符号反転
/// return "-q"
#[inline(always)]
pub fn negate(q: Quaternion<f64>) -> Quaternion<f64> {
    ( -q.0, negate_vec(q.1) )
}

/// ハミルトン積のために定義した，特別なベクトル同士の積（非可換）．
/// この積は，二つのベクトルを純虚四元数と見なした場合のハミルトン積と等しい．
/// "ab ≡ -a・b + a×b" (!= ba)
#[inline(always)]
pub fn mul_vec(a: Vector3<f64>, b: Vector3<f64>) -> Quaternion<f64> {
    ( -dot_vec(a, b), cross_vec(a, b) )
}

/// ハミルトン積
/// 積の順序は "ab"(!=ba)
/// Hamilton product
/// The product order is "ab"(!= ba)
#[inline(always)]
pub fn mul(a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
    let tmp0 = scale_vec(a.0, b.1);
    let tmp1 = scale_vec(b.0, a.1);
    let term1 = ( a.0 * b.0, add_vec(tmp0, tmp1) );
    let term2 = mul_vec(a.1, b.1);
    add(term1, term2)
}

/// 共役四元数を求める
/// Compute the conjugate quaternion
#[inline(always)]
pub fn conj(a: Quaternion<f64>) -> Quaternion<f64> {
    ( a.0, negate_vec(a.1) )
}

/// 逆四元数を求める
/// Compute the inverse quaternion
#[inline(always)]
pub fn inverse(a: Quaternion<f64>) -> Quaternion<f64> {
    scale( dot(a, a).recip(), conj(a) )
}

/// ネイピア数eのベクトル冪
/// Exponential of Vector3.
#[inline(always)]
pub fn exp_vec(a: Vector3<f64>) -> Quaternion<f64> {
    let norm = norm_vec(a);
    if norm == 0.0 {
        return IDENTITY;
    }
    let f = norm.sin_cos();
    ( f.1, scale_vec( f.0 / norm, a ) )
}

/// ネイピア数eの四元数冪
/// Exponential of Quaternion.
#[inline(always)]
pub fn exp(a: Quaternion<f64>) -> Quaternion<f64> {
    scale( a.0.exp(), exp_vec(a.1) )
}

/// 四元数の自然対数
/// Natural logarithm of quaternion.
#[inline(always)]
pub fn ln(a: Quaternion<f64>) -> Quaternion<f64> {
    let norm = norm(a);
    let norm_vec = norm_vec(a.1);
    if norm_vec == 0.0 {
        return (norm.ln(), ZERO_VECTOR);
    }
    let tmp = acos_safe(a.0 / norm);
    ( norm.ln(), scale_vec(tmp / norm_vec, a.1) )
}

/// 四元数の冪乗
/// The power of quaternion.
#[inline(always)]
pub fn power(a: Quaternion<f64>, t: f64) -> Quaternion<f64> {
    let coef = norm(a).powf(t);
    let f = ( t * acos_safe(a.0) ).sin_cos();
    let norm_vec = norm_vec(a.1);
    if norm_vec == 0.0 {
        return (coef * f.1, ZERO_VECTOR);
    }
    let q_v = scale_vec(f.0 / norm_vec, a.1);
    scale( coef, (f.1, q_v) )
}

/// 位置ベクトルの回転
/// r' = q r q*  (||q|| = 1)
/// 引数"q"のノルムは1でなければならない．
/// The norm of argument "q" must be 1.
#[inline(always)]
pub fn vector_rotation(q: Quaternion<f64>, r: Vector3<f64>) -> Vector3<f64> {
    let dot = dot_vec(q.1, r);
    let cross = cross_vec(q.1, r);
    let term1 = add_vec( scale_vec(q.0, r),   scale_vec(2.0, cross) );
    let term2 = add_vec( scale_vec(dot, q.1), cross_vec(q.1, cross) );
    scale_add_vec(q.0, term1, term2)
}

/// 座標系の回転
/// r' = q* r q  (||q|| = 1)
/// 引数"q"のノルムは1でなければならない．
/// The norm of argument "q" must be 1.
#[inline(always)]
pub fn frame_rotation(q: Quaternion<f64>, r: Vector3<f64>) -> Vector3<f64> {
    let dot = dot_vec(q.1, r);
    let cross = cross_vec(q.1, r);
    let term1 = add_vec( scale_vec(q.0, r),  scale_vec(-2.0, cross) );
    let term2 = add_vec( scale_vec(dot, q.1), cross_vec(q.1, cross) );
    scale_add_vec(q.0, term1, term2)
}

/// ベクトル "a" を ベクトル "b" へ最短距離で回転させる四元数を求める．
/// 零ベクトルを入力した場合は，恒等四元数を返す．
/// Find a quaternion to rotate from vector "a" to "b".
/// If you enter a zero vector, it returns an identity quaternion.
/// 0 <= t <= 1
#[inline(always)]
pub fn rotate_a_to_b(a: Vector3<f64>, b: Vector3<f64>, t: f64) -> Quaternion<f64> {
    let a = normalize_vec(a);
    let b = normalize_vec(b);
    let axis = cross_vec(a, b);
    let angle = acos_safe( dot_vec(a, b) );
    from_axis_angle(axis, angle * t)
}

/// 機体座標系の姿勢変化を積分して，引数に渡した四元数を更新する．
/// Integrate attitude change of body coordinate frame, 
/// and update Quaternion passed to the argument.
/// accel[rad/sec]
/// dt[sec]
#[inline(always)]
pub fn integration(q: Quaternion<f64>, accel: Vector3<f64>, dt: f64) -> Quaternion<f64> {    
    let tmp = scale_vec(dt * 0.5, accel);
    let dq  = exp_vec(tmp);
    mul(dq, q)
}

/// オイラー法
/// Euler method
/// 機体座標系の姿勢変化を積分して，引数に渡した四元数を更新する．
/// Integrate attitude change of body coordinate frame, 
/// and update Quaternion passed to the argument.
/// accel[rad/sec]
/// dt[sec]
#[inline(always)]
pub fn integration_euler(q: Quaternion<f64>, accel: Vector3<f64>, dt: f64) -> Quaternion<f64> {
    let tmp = mul( (0.0, accel), q );
    normalize( scale_add(dt*0.5, tmp, q) )
}

/// 線形補間
/// 引数"a"から"b"への経路を補完する四元数を生成する．
/// 引数t(0 <= t <= 1) は補間パラメータ．
/// "a"と"b"のノルムは必ず1でなければならない．
/// Linear interpolation
/// Generate a quaternion that interpolate the route from "a" to "b".
/// The argument t(0 <= t <= 1) is the interpolation parameter.
/// The norm of "a" and "b" must be 1.
#[inline(always)]
pub fn lerp(a: Quaternion<f64>, b: Quaternion<f64>, t: f64) -> Quaternion<f64> {
    let term1 = scale(1.0 - t, a);
    let term2 = scale(t, b);
    normalize( add(term1, term2) )
}

/// 球状線形補間
/// "a"から"b"への経路を補完する四元数を生成する．
/// 引数t(0 <= t <= 1) は補間パラメータ．
/// "a"と"b"のノルムは必ず1でなければならない．
/// Spherical linear interpolation
/// Generate a quaternion that interpolate the route from "a" to "b".
/// The norm of "a" and "b" must be 1.
/// The argument t(0 <= t <= 1) is the interpolation parameter.
#[inline(always)]
pub fn slerp(a: Quaternion<f64>, mut b: Quaternion<f64>, t: f64) -> Quaternion<f64> {
    // 最短経路で補間
    let mut dot = dot(a, b);
    if dot.is_sign_negative() == true {
        b = negate(b);
        dot = -dot;
    }
    // If the inputs are too close for comfort, linearly interpolate.
    if dot > THRESHOLD {
        return lerp(a, b, t);
    }
    // selrp
    let omega = dot.acos();  // Angle between the two quaternion
    let s1 = ( (1.0 - t) * omega ).sin();
    let s2 = (t * omega).sin();
    let term1 = scale(s1, a);
    let term2 = scale(s2, b);
    scale( omega.sin().recip(), add(term1, term2) )
}

/// 四元数の冪乗を用いたSlerpアルゴリズム．
/// 実装方法が違うだけで，計算結果はslerp関数と同じ．
/// "a"と"b"のノルムは必ず1でなければならない．
/// Slerp algorithm using quaternion powers.
/// The calculation result is the same as the slerp() function.
/// The norm of "a" and "b" must be 1.
/// "a" --> "b"
#[inline(always)]
pub fn slerp_1(a: Quaternion<f64>, mut b: Quaternion<f64>, t: f64) -> Quaternion<f64> {
    // 最短経路で補間
    let mut dot = dot(a, b);
    if dot.is_sign_negative() == true {
        b = negate(b);
        dot = -dot;
    }
    // lerp
    if dot > THRESHOLD {
        return lerp(a, b, t);
    }
    // slerp
    let tmp = power( mul( b, conj(a) ), t );
    mul(tmp, a)
}

/// 定義域外の値をカットして未定義動作を防ぐ
#[inline(always)]
fn acos_safe(s: f64) -> f64 {
    if s >= 1.0 {  // Avoid undefined behavior
        return 0.0;
    } else if s <= -1.0 {
        return PI;
    }
    s.acos()
}
