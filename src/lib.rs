// Quaternion Libraly
// f64 only

pub type Vector3<T> = [T; 3];
pub type Quaternion<T> = (T, Vector3<T>);
pub type DirectionCosines<T> = [Vector3<T>; 3];

const PI: f64 = std::f64::consts::PI;
const FRAC_PI_2: f64 = std::f64::consts::FRAC_PI_2;  // π/2
const THRESHOLD: f64 = 0.9995;


/// 恒等四元数を生成．
/// Generate identity quaternion.
#[inline(always)]
pub fn id() -> Quaternion<f64> {
    (1.0, [0.0; 3])
}

/// 回転角と軸ベクトルを指定して四元数を生成．
/// axisは単位ベクトルでなくても良い．
/// 零ベクトルを入力した場合は，恒等四元数を返す．
/// Generate quaternion by specifying rotation angle and axis vector.
/// The "axis" need not be unit vector.
/// If you enter a zero vector, it returns an identity quaternion.
/// angle[rad]
#[inline(always)]
pub fn from_axis_angle(axis: Vector3<f64>, angle: f64) -> Quaternion<f64> {
    let norm = norm_vec(axis);
    if (norm == 0.0) || (angle == 0.0) {
        return id();
    }
    let n = scale_vec( norm.recip(), axis );
    let f = (angle * 0.5).sin_cos();
    ( f.1, scale_vec(f.0, n) )
}

/// 方向余弦行列から四元数を生成．
#[inline(always)]
pub fn from_direction_cosines(m: DirectionCosines<f64>) -> Quaternion<f64> {
    let q0 = (m[0][0] + m[1][1] + m[2][2] + 1.0).sqrt() * 0.5;
    let tmp = (4.0 * q0).recip();  // reciprocal
    let q1 = (m[1][2] - m[2][1]) * tmp;
    let q2 = (m[2][0] - m[0][2]) * tmp;
    let q3 = (m[0][1] - m[1][0]) * tmp;
    (q0, [q1, q2, q3])
}

/// オイラー角[rad]から四元数を生成．
#[inline(always)]
pub fn from_euler_angles(roll: f64, pitch: f64, yaw: f64) -> Quaternion<f64> {
    let alpha = yaw   * 0.5;
    let beta  = pitch * 0.5;
    let gamma = roll  * 0.5;
    // Compute these value only once
    let (sin_alpha, cos_alpha) = alpha.sin_cos();
    let (sin_beta,  cos_beta)  = beta.sin_cos();
    let (sin_gamma, cos_gamma) = gamma.sin_cos();

    let q0 = cos_alpha * cos_beta * cos_gamma + sin_alpha * sin_beta * sin_gamma;
    let q1 = cos_alpha * cos_beta * sin_gamma - sin_alpha * sin_beta * cos_gamma;
    let q2 = cos_alpha * sin_beta * cos_gamma + sin_alpha * cos_beta * sin_gamma;
    let q3 = sin_alpha * cos_beta * cos_gamma - cos_alpha * sin_beta * sin_gamma;
    (q0, [q1, q2, q3])
}

/// 回転を表す四元数から，回転軸（単位ベクトル）と軸回りの回転角[rad]を取り出す．
/// return "(axis, angle)"
#[inline(always)]
pub fn to_axis_angle(q: Quaternion<f64>) -> (Vector3<f64>, f64) {
    let q = normalize(q);
    let axis = get_unit_vector(q);
    let angle = 2.0 * acos_safe(q.0);
    (axis, angle)
}

/// 四元数を方向余弦行列（回転行列）に変換
#[inline(always)]
pub fn to_direction_cosines((q0, [q1, q2, q3]): Quaternion<f64>) -> DirectionCosines<f64> {
    let m11 = (q0*q0 + q1*q1).mul_add(2.0, -1.0);
    let m12 = (q1*q2 + q0*q3) * 2.0;
    let m13 = (q1*q3 - q0*q2) * 2.0;
    let m21 = (q1*q2 - q0*q3) * 2.0;
    let m22 = (q0*q0 + q2*q2).mul_add(2.0, -1.0);
    let m23 = (q2*q3 + q0*q1) * 2.0;
    let m31 = (q1*q3 + q0*q2) * 2.0;
    let m32 = (q2*q3 - q0*q1) * 2.0;
    let m33 = (q0*q0 + q3*q3).mul_add(2.0, -1.0);

    [
        [m11, m12, m13],
        [m21, m22, m23],
        [m31, m32, m33]
    ]
}

/// 四元数をオイラー角[rad]に変換
/// Quaternion --> [roll, pitch, yaw]
#[inline(always)]
pub fn to_euler_angles((q0, [q1, q2, q3]): Quaternion<f64>) -> Vector3<f64> {
    let m11 = (q0*q0 + q1*q1).mul_add(2.0, -1.0);
    let m12 = (q1*q2 + q0*q3) * 2.0;
    let m13 = (q1*q3 - q0*q2) * 2.0;
    let m23 = (q2*q3 + q0*q1) * 2.0;
    let m33 = (q0*q0 + q3*q3).mul_add(2.0, -1.0);

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
pub fn get_unit_vector(q: Quaternion<f64>) -> Vector3<f64> {
    if q.0 == 1.0 {
        return [0.0; 3];  // ゼロ除算回避
    }
    let tmp = q.0.mul_add(-q.0, 1.0);
    scale_vec( tmp.sqrt().recip(), q.1 )
}

/// ベクトルの内積
/// Dot product of vector
#[inline(always)]
pub fn dot_vec(a: Vector3<f64>, b: Vector3<f64>) -> f64 {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

/// 四元数の内積
/// Dot product of quaternion
/// a・b
#[inline(always)]
pub fn dot(a: Quaternion<f64>, b: Quaternion<f64>) -> f64 {
    a.0.mul_add( b.0, dot_vec(a.1, b.1) )
}

/// ベクトル積（外積）
/// Cross product
/// a×b
#[inline(always)]
pub fn cross_vec(a: Vector3<f64>, b: Vector3<f64>) -> Vector3<f64> {
    let vec0 = a[1]*b[2] - a[2]*b[1];
    let vec1 = a[2]*b[0] - a[0]*b[2];
    let vec2 = a[0]*b[1] - a[1]*b[0];
    [vec0, vec1, vec2]
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

// ベクトルの減算
// Vector subtraction
// Calculate "a - b"
#[inline(always)]
pub fn sub_vec(a: Vector3<f64>, b: Vector3<f64>) -> Vector3<f64> {
    [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ]
}

// 四元数の減算
// Quaternion subtraction
// Calculate "a - b"
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
/// CPUがFMA命令をサポートしている場合，
/// "add_vec(scale_vec(s, a), b)" よりも高速に計算出来る．
/// If the CPU supports FMA instructions,
/// this is faster than "add_vec(scale_vec(s, a), b)".
#[inline(always)]
pub fn scale_add_vec(s: f64, a: Vector3<f64>, b: Vector3<f64>) -> Vector3<f64> {
    let vec0 = s.mul_add(a[0], b[0]);
    let vec1 = s.mul_add(a[1], b[1]);
    let vec2 = s.mul_add(a[2], b[2]);
    [vec0, vec1, vec2]
}

/// compute "s*a + b"
/// CPUがFMA命令をサポートしている場合，
/// "add(scale(s, a), b)" よりも高速に計算出来る．
/// If the CPU supports FMA instructions,
/// this is faster than "add_vec(scale_vec(s, a), b)".
#[inline(always)]
pub fn scale_add(s: f64, a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
    let q0 = s.mul_add(a.0, b.0);
    let vec = scale_add_vec(s, a.1, b.1);
    (q0, vec)
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
/// Normalization
#[inline(always)]
pub fn normalize_vec(r: Vector3<f64>) -> Vector3<f64> {
    let norm = norm_vec(r);
    if norm == 0.0 {
        return [0.0; 3];  // ゼロ除算回避
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

/// ハミルトン積のために定義した，特別なベクトル同士の積
/// ab ≡ -a・b + a×b
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
    let q_s = a.0 * b.0 - dot_vec(a.1, b.1);
    let tmp = scale_add_vec( b.0, a.1, cross_vec(a.1, b.1) );
    let q_v = scale_add_vec( a.0, b.1, tmp );
    (q_s, q_v)
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
        return id();
    }
    let n = scale_vec( norm.recip(), a );
    let f = norm.sin_cos();
    ( f.1, scale_vec(f.0, n) )
}

/// ネイピア数eの四元数冪
/// Exponential of Quaternion.
#[inline(always)]
pub fn exp(a: Quaternion<f64>) -> Quaternion<f64> {
    scale( a.0.exp(), exp_vec(a.1) )
}

/// 四元数の自然対数
/// Natural logarithm of quaternion
#[inline(always)]
pub fn ln(a: Quaternion<f64>) -> Quaternion<f64> {
    let norm = norm(a);
    let n = normalize_vec(a.1);
    let tmp = acos_safe(a.0 / norm);
    ( norm.ln(), scale_vec(tmp, n) )
}

/// 四元数の冪乗
/// The power of quaternion.
#[inline(always)]
pub fn power(a: Quaternion<f64>, t: f64) -> Quaternion<f64> {
    let coef = norm(a).powf(t);
    let n = normalize_vec(a.1);
    let f = ( t * acos_safe(a.0) ).sin_cos();
    let q_v = scale_vec(f.0, n);
    scale( coef, (f.1, q_v) )
}

/// 位置ベクトルの回転
/// r' = q r q*
#[inline(always)]
pub fn vector_rotation(q: Quaternion<f64>, r: Vector3<f64>) -> Vector3<f64> {
    let q = normalize(q);
    let tmp1  = scale_vec( 2.0, cross_vec(q.1, r) );
    let term1 = scale_add_vec(q.0, r, tmp1);
    let tmp2  = cross_vec( q.1, cross_vec(q.1, r) );
    let term2 = scale_add_vec( dot_vec(r, q.1), q.1, tmp2 );
    scale_add_vec(q.0, term1, term2)
}

/// 座標系の回転
/// r' = q* r q
#[inline(always)]
pub fn coordinate_rotation(q: Quaternion<f64>, r: Vector3<f64>) -> Vector3<f64> {
    let q = normalize(q);
    let tmp1  = scale_vec( 2.0, cross_vec(r, q.1) );
    let term1 = scale_add_vec(q.0, r, tmp1);
    let tmp2  = cross_vec( q.1, cross_vec(q.1, r) );
    let term2 = scale_add_vec( dot_vec(r, q.1), q.1, tmp2 );
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
    let theta = acos_safe( dot_vec(a, b) );
    from_axis_angle(axis, theta * t)
}

/// 機体座標系の姿勢変化を積分して，引数に渡した四元数を更新する．
/// Integrate attitude change of body coordinate system, 
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
/// Integrate attitude change of body coordinate system, 
/// and update Quaternion passed to the argument.
/// accel[rad/sec]
/// dt[sec]
#[inline(always)]
pub fn integration_euler(q: Quaternion<f64>, accel: Vector3<f64>, dt: f64) -> Quaternion<f64> {
    let tmp = mul( (0.0, accel), q );
    let q_int = scale_add(dt*0.5, tmp, q);
    normalize( q_int )
}

/// 線形補間
/// 引数"a"から"b"への経路を補完する四元数を生成する．
/// 引数t(0 <= t <= 1) は補間パラメータ．
/// Linear interpolation
/// Generate a quaternion that interpolate the route from "a" to "b".
/// The argument t(0 <= t <= 1) is the interpolation parameter.
#[inline(always)]
pub fn lerp(a: Quaternion<f64>, b: Quaternion<f64>, t: f64) -> Quaternion<f64> {
    let q = scale_add( 1.0 - t, a, scale(t, b) );
    normalize(q)
}

/// 球状線形補間
/// "a"から"b"への経路を補完する四元数を生成する．
/// 引数t(0 <= t <= 1) は補間パラメータ．
/// Spherical linear interpolation
/// Generate a quaternion that interpolate the route from "a" to "b".
/// The argument t(0 <= t <= 1) is the interpolation parameter.
#[inline(always)]
pub fn slerp(a: Quaternion<f64>, b: Quaternion<f64>, t: f64) -> Quaternion<f64> {
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
    if dot > THRESHOLD {
        return lerp(a, b, t);
    }
    // selrp
    let omega = acos_safe(dot);  // Angle between the two quaternion
    let s1 = ( (1.0 - t) * omega ).sin();
    let s2 = (t * omega).sin();
    let tmp = scale_add( s1, a, scale(s2, b) );
    scale( omega.sin().recip(), tmp )
}

/// 四元数の冪乗を用いたSlerp
/// Sherical linear interpolation. 
/// Use quaternion's exponential.
/// "a" --> "b".
#[inline(always)]
pub fn slerp_1(a: Quaternion<f64>, b: Quaternion<f64>, t: f64) -> Quaternion<f64> {
    let a = normalize(a);
    let mut b = normalize(b);
    // 最短経路で補完
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
fn asin_safe(s: f64) -> f64 {
    if s >= 1.0 {  // Avoid undefined behavior
        return  FRAC_PI_2;
    } else if s <= -1.0 {
        return -FRAC_PI_2;
    }
    s.asin()
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
