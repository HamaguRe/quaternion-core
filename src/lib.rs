// Quaternion Libraly
// f32 & f64
//
// Versorは「回転子」を意味し，ノルムは1に制限される．
// Versor means the "rotator", the norm is limited to 1.

pub type Vector3<T> = [T; 3];  // [i, j, k]
pub type Quaternion<T> = (T, Vector3<T>);  // (1, [i, j, k])
pub type DCM<T> = [Vector3<T>; 3];  // Direction Cosine Matrix

use num_traits::float::{Float, FloatConst};


/// 回転角[rad]と軸ベクトルを指定してVersorを生成する．
/// axisは単位ベクトルでなくても良い．
/// 零ベクトルを入力した場合は，恒等四元数を返す．
/// Generate Versor by specifying rotation angle[rad] and axis vector.
/// The "axis" need not be unit vector.
/// If you enter a zero vector, it returns an identity quaternion.
#[inline(always)]
pub fn from_axis_angle<T>(axis: Vector3<T>, angle: T) -> Quaternion<T>
where T: Float + FloatConst {
    #[allow(non_snake_case)]
    let TWO_PI  = cast::<T>(2.0) * T::PI();

    let norm_vec = norm_vec(axis);
    if norm_vec < T::epsilon() {  // ゼロ除算回避
        IDENTITY()
    } else {
        let tmp = angle % TWO_PI;  // limit to (-2π, 2π)
        let f = ( tmp * cast(0.5) ).sin_cos();
        ( f.1, scale_vec(f.0 / norm_vec, axis) )
    }
}

/// Versorから回転軸（単位ベクトル）と軸回りの回転角[rad]を求める．
/// Compute the rotation axis (unit vector) and the rotation angle[rad] 
/// around the axis from the versor.
/// return "(axis, angle)"
/// -π <= angle <= π
#[inline(always)]
pub fn to_axis_angle<T>(q: Quaternion<T>) -> (Vector3<T>, T)
where T: Float + FloatConst {
    let norm_vec = norm_vec(q.1);
    if norm_vec < T::epsilon() {
        ( ZERO_VECTOR(), T::zero() )
    } else {
        let axis = scale_vec( norm_vec.recip(), q.1 );
        let angle = copysign( cast::<T>(2.0) * asin_safe(norm_vec), q.0 );
        (axis, angle)
    }
}

/// 位置ベクトル回転(q v q*)を表す方向余弦行列をVersorに変換．
/// 座標系回転(q* v q)を表す方向余弦行列を変換する場合には，
/// let q = conj( from_dcm(dcm) );
/// とする．
/// Generate versor from direction cosine matrix 
/// representing rotation of position vector.
#[inline(always)]
pub fn from_dcm<T>(m: DCM<T>) -> Quaternion<T>
where T: Float {
    let half: T = cast(0.5);

    // ゼロ除算を避けるために，4通りの式で求めたうちの最大値を係数として使う．
    let (index, max_num) = max4([
        m[0][0] + m[1][1] + m[2][2],
        m[0][0] - m[1][1] - m[2][2],
       -m[0][0] + m[1][1] - m[2][2],
       -m[0][0] - m[1][1] + m[2][2],
    ]);

    let tmp = ( max_num + T::one() ).sqrt();
    let coef = (cast::<T>(2.0) * tmp).recip();
    match index {
        0 => {
            let q0 = half * tmp;
            let q1 = (m[2][1] - m[1][2]) * coef;
            let q2 = (m[0][2] - m[2][0]) * coef;
            let q3 = (m[1][0] - m[0][1]) * coef;
            (q0, [q1, q2, q3])
        },
        1 => {
            let q1 = half * tmp;
            let q0 = (m[2][1] - m[1][2]) * coef;
            let q2 = (m[0][1] + m[1][0]) * coef;
            let q3 = (m[0][2] + m[2][0]) * coef;
            (q0, [q1, q2, q3])
        },
        2 => {
            let q2 = half * tmp;
            let q0 = (m[0][2] - m[2][0]) * coef;
            let q1 = (m[0][1] + m[1][0]) * coef;
            let q3 = (m[1][2] + m[2][1]) * coef;
            (q0, [q1, q2, q3])
        },
        3 => {
            let q3 = half * tmp;
            let q0 = (m[1][0] - m[0][1]) * coef;
            let q1 = (m[0][2] + m[2][0]) * coef;
            let q2 = (m[1][2] + m[2][1]) * coef;
            (q0, [q1, q2, q3])
        },
        _ => unreachable!(),
    }
}

/// 位置ベクトル回転(q v q*)を表すVersorを，方向余弦行列に変換．
/// 座標系回転(q* v q)を表すVersorを変換する場合には，
/// let dcm = to_dcm( conj(q) );
/// とする．
#[inline(always)]
pub fn to_dcm<T>(q: Quaternion<T>) -> DCM<T>
where T: Float {
    let one = T::one();
    let two = cast(2.0);

    let (q0, [q1, q2, q3]) = q;
    // Compute these value only once.
    let q0_q0 = q0 * q0;
    let q0_q1 = q0 * q1;
    let q0_q2 = q0 * q2;
    let q0_q3 = q0 * q3;
    let q1_q2 = q1 * q2;
    let q1_q3 = q1 * q3;
    let q2_q3 = q2 * q3;

    let m11 = (q0_q0 + q1*q1).mul_add(two, -one);
    let m12 = (q1_q2 - q0_q3) * two;
    let m13 = (q1_q3 + q0_q2) * two;
    let m21 = (q1_q2 + q0_q3) * two;
    let m22 = (q0_q0 + q2*q2).mul_add(two, -one);
    let m23 = (q2_q3 - q0_q1) * two;
    let m31 = (q1_q3 - q0_q2) * two;
    let m32 = (q2_q3 + q0_q1) * two;
    let m33 = (q0_q0 + q3*q3).mul_add(two, -one);

    [
        [m11, m12, m13],
        [m21, m22, m23],
        [m31, m32, m33],
    ]
}

/// 位置ベクトル回転を表すVersorをz-y-x系のオイラー角[rad]に変換する．
/// 座標系回転を表すVersorを変換する場合には，
/// let euler_angles = to_euler_angle( conj(q) );
/// とする．
/// Aerospace angle/axis sequences is [yaw, pitch, roll] / [z, y, x].
/// pitch angle is limited to [-π/2, π/2]．
/// return: [yaw, pitch, roll]
#[inline(always)]
pub fn to_euler_angle<T>(q: Quaternion<T>) -> Vector3<T>
where T: Float + FloatConst {
    let [
        [m11, m12, m13],
        [  _, m22, m23],
        [  _, m32, m33],
    ] = to_dcm(q);

    let yaw;
    let pitch;
    let roll;
    if m13.abs() < (T::one() - T::epsilon()) {  // ジンバルロックが起きていない場合
        pitch = m13.asin();
        yaw  = (-m12 / m11).atan();
        roll = (-m23 / m33).atan();
    } else {
        pitch = copysign(T::FRAC_PI_2(), m13);
        yaw  = T::zero();
        roll = (m32 / m22).atan();
    }

    [yaw, pitch, roll]
}

/// 方向余弦行列を用いてベクトルを回転させる．
#[inline(always)]
pub fn matrix_product<T>(m: DCM<T>, v: Vector3<T>) -> Vector3<T>
where T: Float {
    [
        dot_vec(m[0], v),
        dot_vec(m[1], v),
        dot_vec(m[2], v),
    ]
}

/// 右手系と左手系の四元数を変換
/// x, z軸の向きはそのままで，y軸と全ての軸回りの回転方向を反転
#[inline(always)]
pub fn system_trans<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    // 実部と，ベクトル部の要素どれか一つの符号を反転させれば良い
    ( -q.0, [ q.1[0], -q.1[1], q.1[2] ] )
}

/// Dot product of vector
#[inline(always)]
pub fn dot_vec<T>(a: Vector3<T>, b: Vector3<T>) -> T
where T: Float {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

/// Dot product of quaternion
#[inline(always)]
pub fn dot<T>(a: Quaternion<T>, b: Quaternion<T>) -> T 
where T: Float {
    a.0 * b.0 + dot_vec(a.1, b.1)
}

/// Cross product
#[inline(always)]
pub fn cross_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T>
where T: Float {
    [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    ]
}

/// Calcurate "a + b"
#[inline(always)]
pub fn add_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T>
where T: Float {
    [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ]
}

/// Calcurate "a + b"
#[inline(always)]
pub fn add<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( a.0 + b.0, add_vec(a.1, b.1) )
}

/// Calculate "a - b"
#[inline(always)]
pub fn sub_vec<T>(a: Vector3<T>, b: Vector3<T>) -> Vector3<T>
where T: Float {
    [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ]
}

/// Calculate "a - b"
#[inline(always)]
pub fn sub<T>(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( a.0 - b.0, sub_vec(a.1, b.1) )
}

/// Multiplication of scalar and vector.
#[inline(always)]
pub fn scale_vec<T>(s: T, v: Vector3<T>) -> Vector3<T>
where T: Float {
    [ s*v[0], s*v[1], s*v[2] ]
}

/// Multiplication of scalar and quaternion.
#[inline(always)]
pub fn scale<T>(s: T, q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( s * q.0, scale_vec(s, q.1) )
}

/// Calculate "s*a + b"
/// CPUがFMA(Fused multiply–add)命令をサポートしている場合，
/// "add_vec(scale_vec(s, a), b)" よりも高速に計算出来る．
/// If the CPU supports FMA(Fused multiply–add) instructions,
/// this is faster than "add_vec(scale_vec(s, a), b)".
#[inline(always)]
pub fn scale_add_vec<T>(s: T, a: Vector3<T>, b: Vector3<T>) -> Vector3<T>
where T: Float {
    [
        s.mul_add(a[0], b[0]),
        s.mul_add(a[1], b[1]),
        s.mul_add(a[2], b[2]),
    ]
}

/// Calculate "s*a + b"
/// CPUがFMA(Fused multiply–add)命令をサポートしている場合，
/// "add(scale(s, a), b)" よりも高速に計算出来る．
/// If the CPU supports FMA(Fused multiply–add) instructions,
/// this is faster than "add(scale(s, a), b)".
#[inline(always)]
pub fn scale_add<T>(s: T, a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( s.mul_add(a.0, b.0), scale_add_vec(s, a.1, b.1) )
}

/// Calculate L2 norm
#[inline(always)]
pub fn norm_vec<T>(v: Vector3<T>) -> T
where T: Float {
    dot_vec(v, v).sqrt()
}

/// Calculate L2 norm
#[inline(always)]
pub fn norm<T>(q: Quaternion<T>) -> T 
where T: Float {
    dot(q, q).sqrt()
}

/// ノルムが1になるように正規化
/// 零ベクトルを入力した場合は零ベクトルを返す．
/// Normalization
/// If you enter a zero vector, it returns a zero vector.
#[inline(always)]
pub fn normalize_vec<T>(v: Vector3<T>) -> Vector3<T>
where T: Float + FloatConst {
    let norm_vec = norm_vec(v);
    if norm_vec < T::epsilon() {
        ZERO_VECTOR()
    } else {
        scale_vec( norm_vec.recip(), v )
    }
}

/// ノルムが1になるように正規化
/// Normalized so that norm is 1
#[inline(always)]
pub fn normalize<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    scale( norm(q).recip(), q )
}

/// 符号反転
/// return "-v"
#[inline(always)]
pub fn negate_vec<T>(v: Vector3<T>) -> Vector3<T>
where T: Float {
    [ -v[0], -v[1], -v[2] ]
}

/// 符号反転
/// return "-q"
#[inline(always)]
pub fn negate<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( -q.0, negate_vec(q.1) )
}

/// 純虚四元数同士のハミルトン積
/// "ab ≡ -a・b + a×b" (!= ba)
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
    let tmp0 = scale_vec(a.0, b.1);
    let tmp1 = scale_vec(b.0, a.1);
    let term1 = ( a.0 * b.0, add_vec(tmp0, tmp1) );
    let term2 = mul_vec(a.1, b.1);
    add(term1, term2)
}

/// 共役四元数を求める
/// Compute the conjugate quaternion
#[inline(always)]
pub fn conj<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    ( q.0, negate_vec(q.1) )
}

/// 逆純虚四元数を求める
/// 零ベクトルを入力した場合は零ベクトルを返す
/// Compute the inverse pure quaternion
pub fn inv_vec<T>(v: Vector3<T>) -> Vector3<T>
where T: Float + FloatConst {
    let dot_vec = dot_vec(v, v);
    if dot_vec < T::epsilon() {
        ZERO_VECTOR()
    } else {
        scale_vec( dot_vec.recip() , negate_vec(v) )
    }
}

/// 逆四元数を求める
/// Compute the inverse quaternion
#[inline(always)]
pub fn inv<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float {
    scale( dot(q, q).recip(), conj(q) )
}

/// ネイピア数を底とする指数函数
/// Exponential function of Vector3
#[inline(always)]
pub fn exp_vec<T>(v: Vector3<T>) -> Quaternion<T>
where T: Float + FloatConst {
    let norm_vec = norm_vec(v);
    if norm_vec < T::epsilon() {
        IDENTITY()
    } else {
        let f = norm_vec.sin_cos();
        ( f.1, scale_vec(f.0 / norm_vec, v) )
    }
}

/// ネイピア数を底とする指数函数
/// Exponential of Quaternion.
#[inline(always)]
pub fn exp<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float + FloatConst {
    scale( q.0.exp(), exp_vec(q.1) )
}

/// 四元数の自然対数
/// Natural logarithm of quaternion.
#[inline(always)]
pub fn ln<T>(q: Quaternion<T>) -> Quaternion<T>
where T: Float + FloatConst {
    let tmp = dot_vec(q.1, q.1);
    let norm = (q.0*q.0 + tmp).sqrt();
    let norm_vec = tmp.sqrt();
    if norm_vec < T::epsilon() {
        ( norm.ln(), ZERO_VECTOR() )
    } else {
        let coef = acos_safe(q.0 / norm) / norm_vec;
        ( norm.ln(), scale_vec(coef, q.1) )
    }
}

/// Versor(単位四元数)の自然対数
/// Versorであることが保証されている場合にはln関数よりも計算量を減らせる．
/// 実部は必ず0になるので省略．
#[inline(always)]
pub fn ln_versor<T>(q: Quaternion<T>) -> Vector3<T>
where T: Float + FloatConst {
    let norm_vec = norm_vec(q.1);
    if norm_vec < T::epsilon() {
        ZERO_VECTOR()
    } else {
        let coef = acos_safe(q.0) / norm_vec;
        scale_vec(coef, q.1)
    }
}

/// 四元数の冪函数
/// The power function of quaternion.
#[inline(always)]
pub fn pow<T>(q: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float + FloatConst {
    let tmp = dot_vec(q.1, q.1);
    let norm = (q.0*q.0 + tmp).sqrt();
    let norm_vec = tmp.sqrt();
    let omega = acos_safe(q.0 / norm);
    let f = (t * omega).sin_cos();
    let coef = norm.powf(t);
    if norm_vec < T::epsilon() {
        ( coef * f.1, ZERO_VECTOR() )
    } else {
        let n = scale_vec(f.0 / norm_vec, q.1);
        scale( coef, (f.1, n) )
    }
}

/// Versor(単位四元数)の冪函数
/// Versorであることが保証されている場合にはpow関数よりも計算量を減らせる．
/// The power function of Versor.
#[inline(always)]
pub fn pow_versor<T>(q: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float + FloatConst {
    let norm_vec = norm_vec(q.1);
    let f = ( t * acos_safe(q.0) ).sin_cos();
    if norm_vec < T::epsilon() {
        ( f.1, ZERO_VECTOR() )
    } else {
        ( f.1, scale_vec(f.0 / norm_vec, q.1) )
    }
}

/// 位置ベクトルの回転
/// q v q*  (||q|| = 1)
/// 引数qは必ずVersor(単位四元数)でなければならない．
/// Argument must be Versor(unit quaternion).
#[inline(always)]
pub fn vector_rotation<T>(q: Quaternion<T>, v: Vector3<T>) -> Vector3<T>
where T: Float {
    let dot = dot_vec(q.1, v);
    let cross = cross_vec(q.1, v);
    let term1 = add_vec( scale_vec(q.0, v), scale_vec(cast(2.0), cross) );
    let term2 = scale_add_vec( dot, q.1, cross_vec(q.1, cross) );
    scale_add_vec(q.0, term1, term2)
}

/// 座標系の回転
/// q* v q  (||q|| = 1)
/// 引数qは必ずVersor(単位四元数)でなければならない．
/// Argument q must be Versor(unit quaternion).
#[inline(always)]
pub fn frame_rotation<T>(q: Quaternion<T>, v: Vector3<T>) -> Vector3<T>
where T: Float {
    let dot = dot_vec(q.1, v);
    let cross = cross_vec(q.1, v);
    let term1 = sub_vec( scale_vec(q.0, v), scale_vec(cast(2.0), cross) );
    let term2 = scale_add_vec( dot, q.1, cross_vec(q.1, cross) );
    scale_add_vec(q.0, term1, term2)
}

/// 位置ベクトル "a" を 位置ベクトル "b" と同じ場所へ最短距離で回転させるVersorを求める．
/// 零ベクトルを入力した場合は，恒等四元数を返す．
/// Find a versor to rotate from vector "a" to "b".
/// If you enter a zero vector, it returns an identity quaternion.
/// 0 <= t <= 1
#[inline(always)]
pub fn rotate_a_to_b<T>(a: Vector3<T>, b: Vector3<T>, t: T) -> Quaternion<T>
where T: Float + FloatConst {
    let axis = cross_vec(a, b);
    let norm_axis = norm_vec(axis);

    if norm_axis < T::epsilon() {
        IDENTITY()
    } else {
        let norm_ab = ( dot_vec(a, a) * dot_vec(b, b) ).sqrt();
        let angle = acos_safe( dot_vec(a, b) / norm_ab );
        let f = ( t * angle * cast(0.5) ).sin_cos();
        ( f.1, scale_vec(f.0 / norm_axis, axis) )
    }
}

/// 線形補間(Linear interpolation)
/// 引数"a"から"b"への経路を補完する四元数を生成する．
/// 引数t(0 <= t <= 1) は補間パラメータ．
/// "a"と"b"のノルムは必ず1でなければならない．
/// Generate a quaternion that interpolate the route from "a" to "b".
/// The argument t(0 <= t <= 1) is the interpolation parameter.
/// The norm of "a" and "b" must be 1.
#[inline(always)]
pub fn lerp<T>(a: Quaternion<T>, b: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float {
    normalize( scale_add(t, sub(b, a), a) )
}

/// 球状線形補間(Spherical linear interpolation)
/// "a"から"b"への経路を補完する四元数を生成する．
/// 引数t(0 <= t <= 1) は補間パラメータ．
/// "a"と"b"のノルムは必ず1でなければならない．
/// Generate a quaternion that interpolate the route from "a" to "b".
/// The norm of "a" and "b" must be 1.
/// The argument t(0 <= t <= 1) is the interpolation parameter.
#[inline(always)]
pub fn slerp<T>(a: Quaternion<T>, mut b: Quaternion<T>, t: T) -> Quaternion<T>
where T: Float {
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
        let term1 = scale(s1, a);
        let term2 = scale(s2, b);
        let coef = (T::one() - dot*dot).sqrt().recip();
        scale( coef, add(term1, term2) )
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

// 定数呼び出し以外に使わないのでエラー処理を省略
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

/// 配列内の最大値を探して，その位置と値を返す．
#[inline(always)]
fn max4<T: Float>(nums: [T; 4]) -> (usize, T) {
    let mut max_num = nums[0];
    let mut index = 0;
    for (i, num) in nums.iter().enumerate().skip(1) {
        if num > &max_num {
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
        copysign( T::FRAC_PI_2(), x )
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