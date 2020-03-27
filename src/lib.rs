// Quaternion Libraly
// f64 only
//
// Versorは「回転子」を意味し，ノルムは1に制限される．
// Versor means the "rotator", the norm is limited to 1.

pub type Vector3<T> = [T; 3];
pub type Quaternion<T> = (T, Vector3<T>);  // (1, [i, j, k])
pub type DCM<T> = [Vector3<T>; 3];  // Direction Cosines Matrix

const EPSILON: f64 = std::f64::EPSILON;
const PI: f64 = std::f64::consts::PI;
const FRAC_PI_2: f64 = std::f64::consts::FRAC_PI_2;  // π/2
const TWO_PI: f64 = 2.0 * PI;
const THRESHOLD: f64 = 0.9995;  // Used in slerp (Approximation error < 0.017%)
const ZERO_VECTOR: Vector3<f64> = [0.0; 3];
pub const IDENTITY: Quaternion<f64> = (1.0, [0.0; 3]);  // Identity Quaternion

// SIMD実装
mod to_fast;
pub use to_fast::*;


/// 回転角[rad]と軸ベクトルを指定して四元数を生成する．
/// axisは単位ベクトルでなくても良い．
/// 零ベクトルを入力した場合は，恒等四元数を返す．
/// Generate quaternion by specifying rotation angle[rad] and axis vector.
/// The "axis" need not be unit vector.
/// If you enter a zero vector, it returns an identity quaternion.
#[inline(always)]
pub fn from_axis_angle(axis: Vector3<f64>, angle: f64) -> Quaternion<f64> {
    let norm = norm_vec(axis);
    if norm < EPSILON {  // ゼロ除算回避
        IDENTITY
    } else {
        let tmp = (angle % TWO_PI).copysign(angle);  // limit to (-2π, 2π)
        let f = (tmp * 0.5).sin_cos();
        ( f.1, scale_vec(f.0 / norm, axis) )
    }
}

/// Versorから回転軸（単位ベクトル）と軸回りの回転角[rad]を求める．
/// Compute the rotation axis (unit vector) and the rotation angle[rad] 
/// around the axis from the versor.
/// return "(axis, angle)"
/// -π <= angle <= π
#[inline(always)]
pub fn to_axis_angle(q: Quaternion<f64>) -> (Vector3<f64>, f64) {
    let norm_vec = norm_vec(q.1);
    if norm_vec < EPSILON {
        (ZERO_VECTOR, 0.0)
    } else {
        let axis = scale_vec( norm_vec.recip(), q.1 );
        let angle = ( 2.0 * asin_safe(norm_vec) ).copysign(q.0);
        (axis, angle)
    }
}

/// 位置ベクトルの回転を表す方向余弦行列からversorを生成．
/// Generate versor from direction cosine matrix 
/// representing rotation of position vector.
#[inline(always)]
pub fn from_dcm_vector(m: DCM<f64>) -> Quaternion<f64> {
    // ゼロ除算が発生しないように，4通りの式で求めたうちの最大値を係数として使う．
    let tmp_list = [
         m[0][0] + m[1][1] + m[2][2],
         m[0][0] - m[1][1] - m[2][2],
        -m[0][0] + m[1][1] - m[2][2],
        -m[0][0] - m[1][1] + m[2][2],
    ];
    let index = max4_index(tmp_list);

    let tmp = (tmp_list[index] + 1.0).sqrt();
    let coef = (2.0 * tmp).recip();
    let q = match index {
        0 => {
            let q0 = 0.5 * tmp;
            let q1 = (m[2][1] - m[1][2]) * coef;
            let q2 = (m[0][2] - m[2][0]) * coef;
            let q3 = (m[1][0] - m[0][1]) * coef;
            (q0, [q1, q2, q3])
        },
        1 => {
            let q1 = 0.5 * tmp;
            let q0 = (m[2][1] - m[1][2]) * coef;
            let q2 = (m[0][1] + m[1][0]) * coef;
            let q3 = (m[0][2] + m[2][0]) * coef;
            (q0, [q1, q2, q3])
        },
        2 => {
            let q2 = 0.5 * tmp;
            let q0 = (m[0][2] - m[2][0]) * coef;
            let q1 = (m[0][1] + m[1][0]) * coef;
            let q3 = (m[1][2] + m[2][1]) * coef;
            (q0, [q1, q2, q3])
        },
        3 => {
            let q3 = 0.5 * tmp;
            let q0 = (m[1][0] - m[0][1]) * coef;
            let q1 = (m[0][2] + m[2][0]) * coef;
            let q2 = (m[1][2] + m[2][1]) * coef;
            (q0, [q1, q2, q3])
        },
        _ => IDENTITY,  // ここに来ることは無い
    };
    normalize(q)
}

/// 位置ベクトル回転を表すVersorを，方向余弦行列（回転行列）に変換．
/// q v q* と同じ回転を表す．
#[inline(always)]
pub fn to_dcm_vector(q: Quaternion<f64>) -> DCM<f64> {
    let (q0, [q1, q2, q3]) = q;
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

/// 座標系回転を表す方向余弦行列からVersorを生成．
/// Generate versor from direction cosine matrix 
/// representing coordinate system rotation.
#[inline(always)]
pub fn from_dcm_frame(m: DCM<f64>) -> Quaternion<f64> {
    conj( from_dcm_vector(m) )
}

/// 座標系回転を表すVersorを，方向余弦行列（回転行列）に変換．
/// q* v q と同じ回転を表す．
#[inline(always)]
pub fn to_dcm_frame(q: Quaternion<f64>) -> DCM<f64> {
    to_dcm_vector( conj(q) )
}

/// 方向余弦行列を用いてベクトルを回転させる．
#[inline(always)]
pub fn matrix_product(m: DCM<f64>, v: Vector3<f64>) -> Vector3<f64> {
    [
        dot_vec(m[0], v),
        dot_vec(m[1], v),
        dot_vec(m[2], v)
    ]
}

/// 右手系と左手系の四元数を変換
/// x, z軸はそのままで，y軸と全ての軸回りの回転方向を反転
#[inline(always)]
pub fn system_trans(q: Quaternion<f64>) -> Quaternion<f64> {
    // 実部と，ベクトル部の要素どれか一つの符号を反転させれば良い
    ( -q.0, [ q.1[0], -q.1[1], q.1[2] ] )
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
#[cfg(not(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "avx")))]
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

/// Calcurate "a + b"
#[inline(always)]
pub fn add_vec(a: Vector3<f64>, b: Vector3<f64>) -> Vector3<f64> {
    [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ]
}

/// Calcurate "a + b"
#[inline(always)]
#[cfg(not(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "avx")))]
pub fn add(a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
    ( a.0 + b.0, add_vec(a.1, b.1) )
}

/// Calculate "a - b"
#[inline(always)]
pub fn sub_vec(a: Vector3<f64>, b: Vector3<f64>) -> Vector3<f64> {
    [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ]
}

/// Calculate "a - b"
#[inline(always)]
#[cfg(not(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "avx")))]
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
#[cfg(not(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "avx")))]
pub fn scale(s: f64, q: Quaternion<f64>) -> Quaternion<f64> {
    ( s * q.0, scale_vec(s, q.1) )
}

/// Compute "s*a + b"
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

/// Compute "s*a + b"
/// CPUがFMA(Fused multiply–add)命令をサポートしている場合，
/// "add(scale(s, a), b)" よりも高速に計算出来る．
/// If the CPU supports FMA(Fused multiply–add) instructions,
/// this is faster than "add(scale(s, a), b)".
#[inline(always)]
#[cfg(not(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "fma")))]
pub fn scale_add(s: f64, a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
    ( s.mul_add(a.0, b.0), scale_add_vec(s, a.1, b.1) )
}

/// L2ノルムを計算
/// Calculate L2 norm
#[inline(always)]
pub fn norm_vec(v: Vector3<f64>) -> f64 {
    dot_vec(v, v).sqrt()
}

/// L2ノルムを計算
/// Calculate L2 norm
#[inline(always)]
pub fn norm(q: Quaternion<f64>) -> f64 {
    dot(q, q).sqrt()
}

/// ノルムが1になるように正規化
/// 零ベクトルを入力した場合は零ベクトルを返す．
/// Normalization
/// If you enter a zero vector, it returns a zero vector.
#[inline(always)]
pub fn normalize_vec(v: Vector3<f64>) -> Vector3<f64> {
    let norm = norm_vec(v);
    if norm < EPSILON {
        ZERO_VECTOR
    } else {
        scale_vec( norm.recip(), v )
    }
}

/// ノルムが1になるように正規化
/// Normalized so that norm is 1
#[inline(always)]
pub fn normalize(q: Quaternion<f64>) -> Quaternion<f64> {
    scale( norm(q).recip(), q )
}

/// 符号反転
/// return "-v"
#[inline(always)]
pub fn negate_vec(v: Vector3<f64>) -> Vector3<f64> {
    [ -v[0], -v[1], -v[2] ]
}

/// 符号反転
/// return "-q"
#[inline(always)]
pub fn negate(q: Quaternion<f64>) -> Quaternion<f64> {
    ( -q.0, negate_vec(q.1) )
}

/// 純虚四元数同士のハミルトン積
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
pub fn conj(q: Quaternion<f64>) -> Quaternion<f64> {
    ( q.0, negate_vec(q.1) )
}

/// 逆純虚四元数を求める
/// 零ベクトルを入力した場合は零ベクトルを返す
/// Compute the inverse pure quaternion
pub fn inv_vec(v: Vector3<f64>) -> Vector3<f64> {
    let dot_vec = dot_vec(v, v);
    if dot_vec < EPSILON {
        ZERO_VECTOR
    } else {
        scale_vec( dot_vec.recip() , negate_vec(v) )
    }
}

/// 逆四元数を求める
/// Compute the inverse quaternion
#[inline(always)]
pub fn inv(q: Quaternion<f64>) -> Quaternion<f64> {
    scale( dot(q, q).recip(), conj(q) )
}

/// ネイピア数を底とする指数函数
/// Exponential function of Vector3
#[inline(always)]
pub fn exp_vec(v: Vector3<f64>) -> Quaternion<f64> {
    let norm = norm_vec(v);
    if norm < EPSILON {
        IDENTITY
    } else {
        let f = norm.sin_cos();
        ( f.1, scale_vec(f.0 / norm, v) )
    }
}

/// ネイピア数を底とする指数函数
/// Exponential of Quaternion.
#[inline(always)]
pub fn exp(q: Quaternion<f64>) -> Quaternion<f64> {
    scale( q.0.exp(), exp_vec(q.1) )
}

/// 四元数の自然対数
/// Natural logarithm of quaternion.
#[inline(always)]
pub fn ln(q: Quaternion<f64>) -> Quaternion<f64> {
    let norm = norm(q);
    let norm_vec = norm_vec(q.1);
    if norm_vec < EPSILON {
        ( norm.ln(), ZERO_VECTOR )
    } else {
        let coef = acos_safe(q.0 / norm) / norm_vec;
        ( norm.ln(), scale_vec(coef, q.1) )
    }
}

/// Versor(単位四元数)の自然対数
/// Versorであることが保証されている場合にはln関数よりも計算量を減らせる．
/// 実部は必ず0になるので省略．
#[inline(always)]
pub fn ln_versor(q: Quaternion<f64>) -> Vector3<f64> {
    let norm_vec = norm_vec(q.1);
    if norm_vec < EPSILON {
        ZERO_VECTOR
    } else {
        let coef = acos_safe(q.0) / norm_vec;
        scale_vec(coef, q.1)
    }
}

/// 四元数の冪函数
/// The power function of quaternion.
#[inline(always)]
pub fn pow(q: Quaternion<f64>, t: f64) -> Quaternion<f64> {
    let norm = norm(q);
    let norm_vec = norm_vec(q.1);
    let omega = acos_safe(q.0 / norm);
    let f = (t * omega).sin_cos();
    let coef = norm.powf(t);
    if norm_vec < EPSILON {
        (coef * f.1, ZERO_VECTOR)
    } else {
        let n = scale_vec(f.0 / norm_vec, q.1);
        scale( coef, (f.1, n) )
    }
}

/// Versor(単位四元数)の冪函数
/// Versorであることが保証されている場合にはpow関数よりも計算量を減らせる．
/// The power function of Versor.
#[inline(always)]
pub fn pow_versor(q: Quaternion<f64>, t: f64) -> Quaternion<f64> {
    let norm_vec = norm_vec(q.1);
    let f = ( t * acos_safe(q.0) ).sin_cos();
    if norm_vec < EPSILON {
        ( f.1, ZERO_VECTOR )
    } else {
        ( f.1, scale_vec(f.0 / norm_vec, q.1) )
    }
}

/// 位置ベクトルの回転
/// q v q*  (||q|| = 1)
/// 引数qは必ずVersor(単位四元数)でなければならない．
/// Argument must be Versor(unit quaternion).
#[inline(always)]
pub fn vector_rotation(q: Quaternion<f64>, v: Vector3<f64>) -> Vector3<f64> {
    let dot = dot_vec(q.1, v);
    let cross = cross_vec(q.1, v);
    let term1 = add_vec( scale_vec(q.0, v),   scale_vec(2.0, cross) );
    let term2 = add_vec( scale_vec(dot, q.1), cross_vec(q.1, cross) );
    scale_add_vec(q.0, term1, term2)
}

/// 座標系の回転
/// q* v q  (||q|| = 1)
/// 引数qは必ずVersor(単位四元数)でなければならない．
/// Argument q must be Versor(unit quaternion).
#[inline(always)]
pub fn frame_rotation(q: Quaternion<f64>, v: Vector3<f64>) -> Vector3<f64> {
    let dot = dot_vec(q.1, v);
    let cross = cross_vec(q.1, v);
    let term1 = sub_vec( scale_vec(q.0, v),   scale_vec(2.0, cross) );
    let term2 = add_vec( scale_vec(dot, q.1), cross_vec(q.1, cross) );
    scale_add_vec(q.0, term1, term2)
}

/// 位置ベクトル "a" を 位置ベクトル "b" と同じ場所へ最短距離で回転させるVersorを求める．
/// 零ベクトルを入力した場合は，恒等四元数を返す．
/// Find a versor to rotate from vector "a" to "b".
/// If you enter a zero vector, it returns an identity quaternion.
/// 0 <= t <= 1
#[inline(always)]
pub fn rotate_a_to_b(mut a: Vector3<f64>, mut b: Vector3<f64>, t: f64) -> Quaternion<f64> {
    a = normalize_vec(a);
    b = normalize_vec(b);
    let axis = cross_vec(a, b);
    let angle = acos_safe( dot_vec(a, b) );
    from_axis_angle(axis, angle * t)
}

/// 機体座標系の姿勢変化を積分して，引数に渡した四元数を更新する．
/// Integrate attitude change of body coordinate frame, 
/// and update Quaternion passed to the argument.
/// omega[rad/sec]: Angular velocity
/// dt[sec]: Time step width
#[inline(always)]
pub fn integration(q: Quaternion<f64>, omega: Vector3<f64>, dt: f64) -> Quaternion<f64> {    
    let dq = exp_vec( scale_vec(dt*0.5, omega) );
    mul(dq, q)
}

/// 近似式(Approximate expression)
/// 機体座標系の姿勢変化を積分して，引数に渡した四元数を更新する．
/// Integrate attitude change of body coordinate frame, 
/// and update quaternion passed to the argument.
/// omega[rad/sec]: Angular velocity
/// dt[sec]: Time step width
#[inline(always)]
pub fn integration_approx(q: Quaternion<f64>, omega: Vector3<f64>, dt: f64) -> Quaternion<f64> {
    let f = mul_vec(omega, q.1);
    let omega_q = ( f.0, scale_add_vec(q.0, omega, f.1) );
    normalize( scale_add(dt*0.5, omega_q, q) )
}

/// 線形補間(Linear interpolation)
/// 引数"a"から"b"への経路を補完する四元数を生成する．
/// 引数t(0 <= t <= 1) は補間パラメータ．
/// "a"と"b"のノルムは必ず1でなければならない．
/// Generate a quaternion that interpolate the route from "a" to "b".
/// The argument t(0 <= t <= 1) is the interpolation parameter.
/// The norm of "a" and "b" must be 1.
#[inline(always)]
pub fn lerp(a: Quaternion<f64>, b: Quaternion<f64>, t: f64) -> Quaternion<f64> {
    let term1 = scale(1.0 - t, a);
    let term2 = scale(t, b);
    normalize( add(term1, term2) )
}

/// 球状線形補間(Spherical linear interpolation)
/// "a"から"b"への経路を補完する四元数を生成する．
/// 引数t(0 <= t <= 1) は補間パラメータ．
/// "a"と"b"のノルムは必ず1でなければならない．
/// Generate a quaternion that interpolate the route from "a" to "b".
/// The norm of "a" and "b" must be 1.
/// The argument t(0 <= t <= 1) is the interpolation parameter.
#[inline(always)]
pub fn slerp(a: Quaternion<f64>, mut b: Quaternion<f64>, t: f64) -> Quaternion<f64> {
    // 最短経路で補間
    let mut dot = dot(a, b);
    if dot.is_sign_negative() {
        b = negate(b);
        dot = -dot;
    }
    // If the distance between quaternions is close enough, use lerp.
    if dot > THRESHOLD {
        lerp(a, b, t)
    } else {
        let omega = dot.acos();  // Angle between the two quaternions.
        let tmp = t * omega;
        let s1 = (omega - tmp).sin();
        let s2 = tmp.sin();
        let term1 = scale(s1, a);
        let term2 = scale(s2, b);
        let coef = (1.0 - dot*dot).sqrt().recip();
        scale( coef, add(term1, term2) )
    }
}

/// 配列内の最大値を探して，その位置を返す
/// return: index of max_num
#[inline(always)]
fn max4_index(nums: [f64; 4]) -> usize {
    let mut max_num = nums[0];
    let mut index = 0;
    for i in 1..4 {
        if nums[i] > max_num {
            max_num = nums[i];
            index = i;
        }
    }
    index
}

/// 定義域外の値をカットして未定義動作を防ぐ
#[inline(always)]
fn asin_safe(x: f64) -> f64 {
    if x.abs() >= 1.0 {
        FRAC_PI_2.copysign(x)
    } else {
        x.asin()
    }
}

/// 定義域外の値をカットして未定義動作を防ぐ
#[inline(always)]
fn acos_safe(x: f64) -> f64 {
    if x.abs() >= 1.0 {
        if x.is_sign_positive() {0.0} else {PI}
    } else {
        x.acos()
    }
}