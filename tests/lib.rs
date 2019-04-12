extern crate quaternion;

use quaternion::*;

const PI: f64 = std::f64::consts::PI;
const EPSILON: f64 = 0.00000001;

fn deg_to_rad(deg: f64) -> f64 {
    (deg / 180.0) * PI
}


#[test]
fn test_euler_quaternion() {
    let roll_ori = deg_to_rad(20.0);
    let pitch_ori = deg_to_rad(18.0);  // 180度を超えると値がおかしくなる
    let yaw_ori = deg_to_rad(0.0);  // 180度，360度で符号が反転する．

    let q = from_euler_angles(roll_ori, pitch_ori, yaw_ori);
    let eulers = to_euler_angles(q);
    
    // 値チェック
    println!("roll: {}, ori: {}", eulers[0], roll_ori);
    assert!( (eulers[0] - roll_ori).abs() < EPSILON );
    println!("pitch: {}, ori: {}", eulers[1], pitch_ori);
    assert!( (eulers[1] - pitch_ori).abs() < EPSILON );
    println!("yaw: {}, ori: {}", eulers[2], yaw_ori);
    assert!( (eulers[2] - yaw_ori).abs() < EPSILON );
}

#[test]
fn test_add() {
    let a = (0.5, [1.0, 1.0, 1.0]);
    let b = (0.5, [1.0, 1.0, 1.0]);
    assert_eq!(add(a, b), (1.0, [2.0, 2.0, 2.0]));
}

#[test]
fn test_norm() {
    let q: Quaternion<f64> = (2f64.powf(-0.5), [0.0, 2f64.powf(-0.5), 0.0]);
    assert!((norm(q) - 1.0).abs() < EPSILON);
}

// 手計算した結果で動作確認
#[test]
fn test_vector_rotation() {
    let r: Vector3<f64> = [2.0, 2.0, 0.0];
    let q: Quaternion<f64> = (2f64.powf(-0.5), [0.0, 2f64.powf(-0.5), 0.0]);
    let result = vector_rotation(q, r);
    assert!( (result[0] - 0.0).abs() < EPSILON);
    assert!( (result[1] - 2.0).abs() < EPSILON);
    assert!( (result[2] + 2.0).abs() < EPSILON);
}

// 二つの方法で，回転を表すクォータニオンの単位ベクトルを求める．
// どっちの方法でも大丈夫だけど，n_2の方が計算量が少ない．
#[test]
fn find_unit_vector() {
    let q: Quaternion<f64> = from_axis_angle([1.0, 4.0, 2.0], PI);
    let omega = q.0.acos();
    let n_1: Vector3<f64> = scale_vec(1.0 / omega.sin(), q.1);
    let n_2: Vector3<f64> = normalize_vec(q.1);
    for i in 0..3 {
        assert!( (n_2[i] - n_1[i]).abs() < EPSILON );
    }
}

/// 角速度を正しく積分出来ているかテスト
#[test]
fn test_integration() {
    let q0 = id();
    let omega: Vector3<f64> = [0.0, PI/2.0, 0.0];  // [rad/s]
    let dt = 1.0;  // [s]
    let mut r: Vector3<f64> = [2.0, 2.0, 0.0];

    let q = integration(q0, omega, dt);
    r = coordinate_rotation(q, r);
    assert!( (r[0] - 0.0).abs() < EPSILON );
    assert!( (r[1] - 2.0).abs() < EPSILON );
    assert!( (r[2] - 2.0).abs() < EPSILON );
}

// 二つの積分方法を試す
#[test]
fn test_integration_method() {
    let omega: Vector3<f64> = [PI/6.0, PI/2.0, PI/4.0];  // [rad/s]
    let dt = 0.001;  // [s]

    let mut q_1 = id();
    let mut q_2 = id();
    for _ in 0..1000 {
        // 理論的には正確な積分を行う．
        // dt間の角速度が一定であれば，dtを大きくしても正確に積分できる．
        q_1 = integration(q_1, omega, dt);

        // この方法は積分結果が超球面上に存在しない．
        // 三角関数を使わないぶん計算量は少ないが，導出方法として正確ではない．
        // dtが大きすぎると誤差が大きくなる．
        //
        // 空間角速度で計算しても同じ結果になる．何故...？
        // ↓
        // 同一軸での回転を表す四元数を合成した場合，積の順序に関わらず結果は等しくなるため．
        // 途中で角速度を変えればズレるはず．
        q_2 = integration_euler(q_2, omega, dt);
    }
    println!("q_1: {:?}", q_1);
    println!("q_2: {:?}", q_2);

    let epsilon = 0.0001;
    assert!( (q_1.0 - q_2.0).abs() < epsilon );
    assert!( (q_1.1[0] - q_2.1[0]).abs() < epsilon );
    assert!( (q_1.1[1] - q_2.1[1]).abs() < epsilon );
    assert!( (q_1.1[2] - q_2.1[2]).abs() < epsilon );
}

// 二つの方法でSlerpを行う．
#[test]
fn test_slerp_eq() {
    let p = from_axis_angle([2.0, 1.2, 3.5], PI);
    let q = from_axis_angle([3.0, 4.5, 1.0], PI/2.0);
    let mut t = 0.1;
    for _i in 0..10 {
        let p_t = slerp(p, q, t);
        let q_t = slerp_1(p, q, t);
        println!("p_slerp: {:?}", p_t);
        println!("q_slerp_1: {:?}", q_t);
        // check
        assert!( (p_t.0 - q_t.0).abs() < EPSILON );
        for i in 0..3 {
            assert!( (p_t.1[i] - q_t.1[i]).abs() < EPSILON );
        }

        t += 0.1;
    }
}