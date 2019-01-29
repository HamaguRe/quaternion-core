pub type Vector3 = [f64; 3];
pub type Quaternion = (f64, Vector3);

const E: f64 = std::f64::consts::E;  // Napier's constant


/// 恒等四元数(identity quaternion)を作成
#[inline(always)]
pub fn id() -> Quaternion {
    (1.0, [0.0; 3])
}

/// theta[rad]
#[inline(always)]
pub fn axis_angle(theta: f64, axis: Vector3) -> Quaternion {
    let axis = normalize_vec(axis);
    let q_s = (theta / 2.0).cos();
    let s = (theta / 2.0).sin();
    let q_v = mul_scalar_vec(s, axis);
    (q_s, q_v)
}

#[inline(always)]
pub fn dot_vec(a: Vector3, b: Vector3) -> f64 {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

#[inline(always)]
pub fn dot(a: Quaternion, b: Quaternion) -> f64 {
    a.0 * b.0 + dot_vec(a.1, b.1)
}

#[inline(always)]
pub fn cross_vec(a: Vector3, b: Vector3) -> Vector3 {
    let vec0 = a[1]*b[2] - a[2]*b[1];
    let vec1 = a[2]*b[0] - a[0]*b[2];
    let vec2 = a[0]*b[1] - a[1]*b[0];
    [vec0, vec1, vec2]
}

#[inline(always)]
pub fn add_vec(a: Vector3, b: Vector3) -> Vector3 {
    [a[0]+b[0], a[1]+b[1], a[2]+b[2]]
}

/// Add two quaternions.
#[inline(always)]
pub fn add(a: Quaternion, b: Quaternion) -> Quaternion {
    ( a.0 + b.0, add_vec(a.1, b.1) )
}

#[inline(always)]
pub fn conj(a: Quaternion) -> Quaternion {
    ( a.0, [-a.1[0], -a.1[1], -a.1[2]] )
}

/// Quatrnion multiply. 
/// This function is calc "ab"(!= ba).
#[inline(always)]
pub fn mul(a: Quaternion, b: Quaternion) -> Quaternion {
    let scalar = a.0 * b.0 - dot_vec(a.1, b.1);
    let vec_1 = mul_scalar_vec(a.0, b.1);
    let vec_2 = mul_scalar_vec(b.0, a.1);
    let vec_3 = cross_vec(a.1, b.1);
    let vector = add_vec( vec_1, add_vec(vec_2, vec_3) );
    (scalar, vector)
}

/// Calc multiply scalar and vector.
#[inline(always)]
pub fn mul_scalar_vec(s: f64, v: Vector3) -> Vector3 {
    [s*v[0], s*v[1], s*v[2]]
}

/// Calc multiply scalar and quaternion.
#[inline(always)]
pub fn mul_scalar_quat(s: f64, a: Quaternion) -> Quaternion {
    ( s * a.0, mul_scalar_vec(s, a.1) )
}

/// Calculate "L2 norm"
#[inline(always)]
pub fn norm(a: Quaternion) -> f64 {
    dot(a, a).sqrt()
}

#[inline(always)]
pub fn norm_vec(r: Vector3) -> f64 {
    dot_vec(r, r).sqrt()
}

#[inline(always)]
pub fn normalize(a: Quaternion) -> Quaternion {
    mul_scalar_quat(1.0 / norm(a), a)
}

#[inline(always)]
pub fn normalize_vec(r: Vector3) -> Vector3 {
    let norm = norm_vec(r);
    if norm == 0.0 {
        return [0.0; 3];
    }
    mul_scalar_vec(1.0 / norm, r)
}

/// 逆クォータニオンを求める．
#[inline(always)]
pub fn inverse(a: Quaternion) -> Quaternion {
    let conj = conj(a);
    let norm_square = dot(a, a);
    mul_scalar_quat(1.0 / norm_square, conj)
}

/// Exponential of Quaternion.
#[inline(always)]
pub fn exp(a: Quaternion) -> Quaternion {
    let coef = E.powf(a.0);  // coefficient（係数）
    let vec_norm = norm_vec(a.1);

    // ゼロ除算対策はnorm_vec()側でやってるからこのif文は無くても良い．
    // 軽量化を求めるか，シンプルさを求めるか...
    if vec_norm == 0.0 {
        return (coef, [0.0; 3]);
    }

    let q_s = vec_norm.cos();
    let n   = normalize_vec(a.1);
    let q_v = mul_scalar_vec(vec_norm.sin(), n);
    mul_scalar_quat( coef, (q_s, q_v) )
}

/// Natural logarithm.
#[inline(always)]
pub fn ln(a: Quaternion) -> Quaternion {
    let vec_norm = norm_vec(a.1);
    let q_s = vec_norm.ln();
    let s = (a.0 / vec_norm).acos();
    let n = normalize_vec(a.1);
    let q_v = mul_scalar_vec(s, n);
    (q_s, q_v)
}

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

/// Coordinate rotate. (r <-- a r a*)
#[inline(always)]
pub fn vector_rotation(a: Quaternion, r: Vector3) -> Vector3 {
    let a = normalize(a);
    let a_conj = conj(a);
    // ベクトルを，スカラー部0のクォータニオンとして計算する．
    let result = mul( a, mul((0.0, r), a_conj) );
    result.1
}

/// Vector rotate. (r <-- a* r a)
#[inline(always)]
pub fn coordinate_rotation(a: Quaternion, r: Vector3) -> Vector3 {
    let a = normalize(a);
    let a_conj = conj(a);
    let result = mul( a_conj, mul((0.0, r), a) );
    result.1
}

/// 角速度で積分して，微小変化を表すクォータニオンを返す．
/// The integration of angular velocity. 
/// omega[rad/s]
/// dt[s]
#[inline(always)]
pub fn integration(omega: Vector3, dt: f64) -> Quaternion {
    let arg = mul_scalar_vec(dt / 2.0, omega);
    exp( (0.0, arg) )
}

/// linear interpolation.
/// Calculate "a" to "b".
#[inline(always)]
pub fn lerp(a: Quaternion, b: Quaternion, t: f64) -> Quaternion {
    let a = normalize(a);
    let b = normalize(b);
    let q_1 = mul_scalar_quat(1.0 - t, a);
    let q_2 = mul_scalar_quat(t, b);
    let result = add(q_1, q_2);
    normalize(result)
}

/// spherical linear interpolation.
/// Calclate "a" to "b".
#[inline(always)]
pub fn slerp(a: Quaternion, b: Quaternion, t: f64) -> Quaternion {
    // Normalize to avoid undefined behavior.
    let a = normalize(a);
    let mut b = normalize(b);

    // Dot production of quaternion.
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
/// Calculate "a" to "b".
#[inline(always)]
pub fn slerp_1(a: Quaternion, b: Quaternion, t: f64) -> Quaternion {
    let a = normalize(a);
    let mut b = normalize(b);
    // Dot production of quaternion.
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
    // slerp
    let a_conj = conj(a);
    let tmp = power( mul(a_conj, b), t);
    mul(a, tmp)
}


/// Tests
#[cfg(test)]
mod test {
    use super::*;
    const PI: f64 = std::f64::consts::PI;
    const EPSILON: f64 = 0.00000001;


    #[test]
    fn test_add() {
        let a = (0.5, [1.0, 1.0, 1.0]);
        let b = (0.5, [1.0, 1.0, 1.0]);
        assert_eq!(add(a, b), (1.0, [2.0, 2.0, 2.0]));
    }

    #[test]
    fn test_norm() {
        let q: Quaternion = (2f64.powf(-0.5), [0.0, 2f64.powf(-0.5), 0.0]);
        assert!((norm(q) - 1.0).abs() < EPSILON);
    }

    // 手計算した結果で動作確認
    #[test]
    fn test_vector_rotation() {
        let r: Vector3 = [2.0, 2.0, 0.0];
        let q: Quaternion = (2f64.powf(-0.5), [0.0, 2f64.powf(-0.5), 0.0]);
        let result = vector_rotation(q, r);
        assert!( (result[0] - 0.0).abs() < EPSILON);
        assert!( (result[1] - 2.0).abs() < EPSILON);
        assert!( (result[2] + 2.0).abs() < EPSILON);
    }

    // 二つの方法でSlerpを行う．
    #[test]
    fn test_slerp() {
        let p = axis_angle(PI, [2.0, 1.2, 3.5]);
        let q = axis_angle(PI/2.0, [3.0, 4.5, 1.0]);
        let mut t = 0.1;
        for _i in 0..10 {
            let p_t = slerp(p, q, t);
            let q_t = slerp_1(p, q, t);
            // check
            assert!( (p_t.0 - q_t.0).abs() < EPSILON );
            for i in 0..3 {
                assert!( (p_t.1[i] - q_t.1[i]).abs() < EPSILON );
            }

            t += 0.1;
        }
    }

    // 二つの方法で，回転を表すクォータニオンの単位ベクトルを求める．
    // どっちの方法でも大丈夫だけど，n_2の方が計算量が少ない．
    #[test]
    fn find_unit_vector() {
        let q: Quaternion = axis_angle(PI, [1.0, 4.0, 2.0]);
        let omega = q.0.acos();
        let n_1: Vector3 = mul_scalar_vec(1.0 / omega.sin(), q.1);
        let n_2: Vector3 = normalize_vec(q.1);
        for i in 0..3 {
            assert!( (n_2[i] - n_1[i]).abs() < EPSILON );
        }
    }

    #[test]
    fn integration_test() {
        let omega: Vector3 = [0.0, PI/2.0, 0.0];  // [rad/s]
        let dt = 1.0;  // [s]
        let dq = integration(omega, dt);
        let mut r: Vector3 = [2.0, 2.0, 0.0];
        r = vector_rotation(dq, r);
        assert!( (r[0] - 0.0).abs() < EPSILON );
        assert!( (r[1] - 2.0).abs() < EPSILON );
        assert!( (r[2] + 2.0).abs() < EPSILON );
    }

}