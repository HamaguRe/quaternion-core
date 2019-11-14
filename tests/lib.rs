use quaternion::*;

const PI: f64 = std::f64::consts::PI;
const EPSILON: f64 = 1e-10;


#[test]
fn test_add() {
    let a: Quaternion<f64> = (0.5, [1.0, 1.0, 1.0]);
    let b: Quaternion<f64> = (0.5, [1.0, 1.0, 1.0]);
    assert_eq!( add(a, b), (1.0, [2.0; 3]) );
}

#[test]
fn test_sub() {
    let a = (0.5, [1.0, 1.0, 1.0]);
    let b = (0.5, [1.0, 1.0, 1.0]);
    assert_eq!( sub(a, b), (0.0, [0.0; 3]) );
}

#[test]
fn test_scale_add() {
    let s = 2.0;
    let a = (0.5, [1.0, 1.0, 1.0]);
    let b = (0.5, [1.0, 1.0, 1.0]);

    let result_1 = scale_add(s, a, b);
    let result_2 = add( scale(s, a), b );
    let diff = sub(result_1, result_2);
    assert!( diff.0.abs() < EPSILON );
    assert!( diff.1[0].abs() < EPSILON );
    assert!( diff.1[1].abs() < EPSILON );
    assert!( diff.1[2].abs() < EPSILON );
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
fn test_negate() {
    let q = (-1.0, [1.0, 2.0, -1.0]);
    let p = negate(q);

    assert_eq!( p, (1.0, [-1.0, -2.0, 1.0]) )
}

// 回転行列に変換してからベクトルを回転させる．
#[test]
fn test_rotate_by_direction_cosines() {
    let r: Vector3<f64> = [2.0, 2.0, 0.0];
    let q = from_axis_angle([0.0, 1.0, 0.0], PI/2.0);

    // 位置ベクトルの回転
    let m = to_direction_cosines_vector(q);
    let result = matrix_product(m, r);
    assert!( (result[0] - 0.0).abs() < EPSILON);
    assert!( (result[1] - 2.0).abs() < EPSILON);
    assert!( (result[2] + 2.0).abs() < EPSILON);

    // 座標系の回転
    let m = to_direction_cosines_frame(q);
    let result = matrix_product(m, r);
    assert!( (result[0] - 0.0).abs() < EPSILON);
    assert!( (result[1] - 2.0).abs() < EPSILON);
    assert!( (result[2] - 2.0).abs() < EPSILON);
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
fn test_frame_rotation() {
    let r = [2.0, 2.0, 0.0];
    let q = from_axis_angle([0.0, 1.0, 0.0], PI/2.0);
    let result = frame_rotation(q, r);
    assert!( (result[0] - 0.0).abs() < EPSILON);
    assert!( (result[1] - 2.0).abs() < EPSILON);
    assert!( (result[2] - 2.0).abs() < EPSILON);
}

// 回転を表す四元数の単位ベクトル（回転軸）を求める．
#[test]
fn test_get_unit_vector() {
    let axis = [1.0, 4.0, 2.0];
    let q = from_axis_angle(axis, PI);
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

// 回転軸と回転角を取り出す
#[test]
fn test_to_axis_angle() {
    let axis = [1.0, 4.0, 2.0];
    let angle = PI;
    let q = from_axis_angle(axis, angle);

    // 軸の方向は分かるが，元の大きさはわからない．
    let n = normalize_vec(axis);
    let f = to_axis_angle(q);
    assert!( (f.1 - angle).abs() < EPSILON );
    for i in 0..3 {
        assert!( (f.0[i] - n[i]).abs() < EPSILON );
    }
}

#[test]
fn test_cross() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];

    let r = cross_vec(r1, r2);
    assert_eq!( [0.0, 0.0, 1.0] , r );
}

#[test]
fn test_exp_ln() {
    let q_0 = IDENTITY;
    let omega = [1.0, 2.0, 1.0];
    let dt = 1.0f64;
    let q_1 = integration(q_0, omega, dt);

    // dt間の姿勢変化から角速度を復元
    // 角速度は一定
    let tmp = mul( q_1, conj(q_0) );
    let omega_re = scale( 2.0/dt, ln(tmp) ).1;

    for i in 0..3 {
        assert!( (omega[i] - omega_re[i]).abs() < EPSILON );
    }
}

/// 角速度を正しく積分出来ているかテスト
#[test]
fn test_integration() {
    let q_0 = IDENTITY;
    let omega = [0.0, PI/2.0, 0.0];  // [rad/s]
    let dt = 1.0;  // [s]
    let mut r = [2.0, 2.0, 0.0];

    let q = integration(q_0, omega, dt);
    r = frame_rotation(q, r);
    assert!( (r[0] - 0.0).abs() < EPSILON );
    assert!( (r[1] - 2.0).abs() < EPSILON );
    assert!( (r[2] - 2.0).abs() < EPSILON );
}

// 二つの積分方法を試す
#[test]
fn test_integration_method() {
    let omega = [PI/6.0, PI/2.0, PI/4.0];  // [rad/s]
    let dt = 0.001;  // [s]

    let mut q_1 = IDENTITY;
    let mut q_2 = IDENTITY;
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

    let epsilon = 1e-6;
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
    let mut t = 0.0;
    for _ in 0..10 {
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