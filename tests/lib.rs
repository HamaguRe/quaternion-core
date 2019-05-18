extern crate quaternion;
use quaternion::*;

const PI: f64 = std::f64::consts::PI;
const EPSILON: f64 = 1e-8;

fn deg_to_rad(deg: f64) -> f64 {
    (deg / 180.0) * PI
}


#[test]
fn test_euler_quaternion() {
    let roll_ori = deg_to_rad(20.0);
    let pitch_ori = deg_to_rad(15.0);  // 180度を超えると値がおかしくなる
    let yaw_ori = deg_to_rad(0.0);  // 180度，360度で符号が反転する．

    let q = from_euler_angles(roll_ori, pitch_ori, yaw_ori);
    let eulers = to_euler_angles(q);
    
    // 値チェック
    assert!( (eulers[0] - roll_ori).abs() < EPSILON );
    assert!( (eulers[1] - pitch_ori).abs() < EPSILON );
    assert!( (eulers[2] - yaw_ori).abs() < EPSILON );
    println!("roll: {}, ori: {}", eulers[0], roll_ori);
    println!("pitch: {}, ori: {}", eulers[1], pitch_ori);
    println!("yaw: {}, ori: {}", eulers[2], yaw_ori);
}

#[test]
fn test_add() {
    let a = (0.5, [1.0, 1.0, 1.0]);
    let b = (0.5, [1.0, 1.0, 1.0]);
    assert_eq!(add(a, b), (1.0, [2.0, 2.0, 2.0]));
}

#[test]
fn test_norm() {
    let q = from_axis_angle([1.0, 2.0, 0.5], 1.5);
    assert!( (norm(q) - 1.0f64).abs() < EPSILON);
}

#[test]
fn test_normalize() {
    let q = (1.0, [2.0, 3.0, 4.0]);
    let q_n = normalize(q);
    assert!( (norm(q_n) - 1.0f64).abs() < EPSILON );
}

#[test]
fn test_sign_inversion() {
    let q = (-1.0, [1.0, 2.0, -1.0]);
    let p = sign_inversion(q);
    assert_eq!(p.0, 1.0);
    assert_eq!(p.1[0], -1.0);
    assert_eq!(p.1[1], -2.0);
    assert_eq!(p.1[2],  1.0);
}

// 手計算した結果で動作確認
#[test]
fn test_vector_rotation() {
    let r: Vector3<f64> = [2.0, 2.0, 0.0];
    let q = from_axis_angle([0.0, 1.0, 0.0], PI/2.0);
    let result = vector_rotation(q, r);
    assert!( (result[0] - 0.0).abs() < EPSILON);
    assert!( (result[1] - 2.0).abs() < EPSILON);
    assert!( (result[2] + 2.0).abs() < EPSILON);
}

#[test]
fn test_coordinate_rotation() {
    let r = [2.0, 2.0, 0.0];
    let q = from_axis_angle([0.0, 1.0, 0.0], PI/2.0);
    let result = coordinate_rotation(q, r);
    assert!( (result[0] - 0.0).abs() < EPSILON);
    assert!( (result[1] - 2.0).abs() < EPSILON);
    assert!( (result[2] - 2.0).abs() < EPSILON);
}

// 回転を表すクォータニオンの単位ベクトル（回転軸）を求める．
#[test]
fn test_get_unit_vector() {
    let axis = [1.0, 4.0, 2.0];
    let q: Quaternion<f64> = from_axis_angle(axis, PI);
    assert!( (norm(q) - 1.0).abs() < EPSILON );

    // 一番計算量が少ないが，引数が単位四元数であることを前提とする．
    let n_1 = get_unit_vector(q);
    // n_1より少し計算量が多いが，必ず単位ベクトルを返す．
    let n_2 = normalize_vec(q.1);

    let n = normalize_vec(axis);
    for i in 0..3 {
        assert!( (n_1[i] - n[i]).abs() < EPSILON );
        assert!( (n_2[i] - n[i]).abs() < EPSILON );
    }
}

// 回転角を取り出す
#[test]
fn test_get_angle() {
    let axis = [1.0, 4.0, 2.0];
    let angle = PI;
    let q: Quaternion<f64> = from_axis_angle(axis, angle);

    let get_angle = get_angle(q);
    assert!( (get_angle - angle).abs() < EPSILON );
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

        // 三角関数を使わないぶん計算量は少ないが，計算方法として正確ではない．
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