extern crate quaternion;

use quaternion::*;

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
    let p = axis_angle([2.0, 1.2, 3.5], PI);
    let q = axis_angle([3.0, 4.5, 1.0], PI/2.0);
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
    let q: Quaternion = axis_angle([1.0, 4.0, 2.0], PI);
    let omega = q.0.acos();
    let n_1: Vector3 = mul_scalar_vec(1.0 / omega.sin(), q.1);
    let n_2: Vector3 = normalize_vec(q.1);
    for i in 0..3 {
        assert!( (n_2[i] - n_1[i]).abs() < EPSILON );
    }
}

#[test]
fn integration_test() {
    let q0 = id();
    let omega: Vector3 = [0.0, PI/2.0, 0.0];  // [rad/s]
    let dt = 1.0;  // [s]
    let mut r: Vector3 = [2.0, 2.0, 0.0];

    let q = vector_integration(q0, omega, dt);
    r = vector_rotation(q, r);
    assert!( (r[0] - 0.0).abs() < EPSILON );
    assert!( (r[1] - 2.0).abs() < EPSILON );
    assert!( (r[2] + 2.0).abs() < EPSILON );
}

#[test]
fn integration_method_test() {
    let q0 = id();
    let omega: Vector3 = [0.0, PI/2.0, 0.0];  // [rad/s]
    let dt = 0.003;  // [s]

    // 理論的には正確な積分を行う．
    // dt間の角速度が一定であれば，dtを大きくしても正確に積分できる．
    let q_1 = vector_integration(q0, omega, dt);
    println!("{:?}", q_1);

    // この方法は積分結果が超球面上に存在しない．
    // 三角関数を使わないぶん計算量は少ないが，導出方法として正確ではない．
    // dtが大きすぎると誤差が大きくなる．
    let q_2 = vector_integration_euler(q0, omega, dt);
    println!("{:?}", q_2);

    assert!( (q_1.0 - q_2.0).abs() < EPSILON );
    assert!( (q_1.1[0] - q_2.1[0]).abs() < EPSILON );
    assert!( (q_1.1[1] - q_2.1[1]).abs() < EPSILON );
    assert!( (q_1.1[2] - q_2.1[2]).abs() < EPSILON );
}