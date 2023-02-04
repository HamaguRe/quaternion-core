// $ cargo test --features fma

use quaternion_core::*;

const PI: f64 = std::f64::consts::PI;
const EPSILON: f64 = 1e-12;


// 3通りの方法でVersorの軸回りの回転角を求める．
//
// 方法1と2は事前に正規化をしておく必要があるが，方法3は実部とベクトル部の比較を
// 行っているため事前の正規化が不要という特徴がある．
// 
// どれで求めても良いが，一般的には[-π, π]の範囲で計算した方が扱い易いのと
// atanの計算時間はasinやacosに比べて2〜3倍長いので方法2がおすすめ．
#[test]
fn test_get_angle() {
    // 適当なVersorを作る
    let axis = [0.0, 1.0, 1.0];
    let angle = -0.5 * PI;
    let q = from_axis_angle(axis, angle);
    assert!( ( norm(q) - 1.0 ).abs() < EPSILON );

    // 方法1
    // これが一番シンプル．実部のみを使うので計算量が少ない．
    // 計算精度が一つの変数に依存してしまうのは良くない...？
    // range: [0, 2π]
    let angle1 = 2.0 * q.0.acos();

    // 方法2
    // 実部の符号を反映することで幾何学的には方法1と同じ結果が得られる．
    // 実部と虚部両方の値を使っているのでなんとなく気持ちが良い．
    // range: (-π, π]
    let angle2 = ( 2.0 * norm(q.1).asin() ).copysign(q.0);

    // 方法3
    // 普通にゼロ除算が発生するが，atanなので計算できる．
    // range: (-π, π)
    let angle3 = 2.0 * (norm(q.1) / q.0).atan();  // これで正しい．atan2だと値域がおかしくなる．

    println!("axis: {:?}", normalize(q.1));
    println!("angle1: {}PI, angle2: {}PI, angle3: {}PI", angle1/PI, angle2/PI, angle3/PI);

    assert!( (angle1 - 0.5*PI).abs() < EPSILON );
    assert!( (angle2 - 0.5*PI).abs() < EPSILON );
    assert!( (angle3 - 0.5*PI).abs() < EPSILON );
}

#[test]
fn test_axis_angle() {
    // to_axis_angle関数はどのような回転角を入れても正常に変換できるが，
    // テストコードの実装上，angleの範囲は(-2π, 2π)とする．
    let axis = normalize([0.0, 1.0, 1.0]);
    let angle = -1.5*PI;
    let q = from_axis_angle(axis, angle);
    println!("q: {:?}", q);

    assert!( (1.0 - norm(q)).abs() < EPSILON );

    let (re_axis, re_angle) = to_axis_angle(q);

    println!("re_axis: {:?}", re_axis);
    println!("re_angle: {}*PI", re_angle / PI);

    // 軸ベクトルのチェック
    if angle.is_sign_positive() {
        assert_eq_vec(re_axis, axis);
    } else {
        // 負の回転角の場合には回転軸が反転する
        assert_eq_vec(re_axis, negate(axis));
    }

    // 回転角のチェック
    if angle > PI {
        let tmp = angle - 2.0*PI;
        assert!( (re_angle - tmp).abs() < EPSILON );
    } else if angle < -PI {
        let tmp = -2.0*PI - angle;
        assert!( (re_angle - tmp).abs() < EPSILON );
    } else {
        assert!( (re_angle - angle.abs()).abs() < EPSILON );
    }
}

#[test]
fn test_dcm() {
    // 色々な値でテスト
    let v = [2.5, 1.0, -3.0];
    let diff = [0.2, -0.1, 2.5];
    let mut axis = [1.0, -0.2, 0.9];
    for i in 0..20 {
        axis = add(axis, diff);
        let q = from_axis_angle(axis, PI * (i as f64));  // Versor
        let dcm = to_dcm(q);
        let q_rest = from_dcm(dcm);

        assert!( ( norm(q_rest) - 1.0 ).abs() < EPSILON );
    
        let rotated_q = point_rotation(q, v);
        let rotated_q_rest = point_rotation(q_rest, v);
        assert_eq_vec(rotated_q, rotated_q_rest);
    }

    // a <--> b の相互変換が正しく行えるかテスト
    let a = [1.0f64, 0.0, 0.0];
    let b = normalize([0.0f64, 1.0, 1.0]);
    let q = rotate_a_to_b(a, b).unwrap();

    let m_a2b = to_dcm(q);
    let b_check = matrix_product(m_a2b, a);
    assert_eq_vec(b, b_check);

    let m_b2a = to_dcm( conj(q) );
    let a_check = matrix_product(m_b2a, b);
    assert_eq_vec(a, a_check);
}

// ベタ実装の計算結果と比較して実装が正しいことを確認．
// オイラー角 --> 四元数の変換では特異点は生じない．
#[test]
fn test_from_intrinsic_euler_angles() {
    use RotationType::*;
    use RotationSequence::*;

    //let angles = [-PI/1.0, PI/4.0, PI/2.5];
    let angles = [1.8*PI, 1.5*PI, 2.0*PI];

    // ------ Proper Euler angles ------ //
    // ZXZ
    let q = from_euler_angles(Intrinsic, ZXZ, angles);
    let z1 = from_axis_angle([0.0, 0.0, 1.0], angles[0]);
    let x  = from_axis_angle([1.0, 0.0, 0.0], angles[1]);
    let z2 = from_axis_angle([0.0, 0.0, 1.0], angles[2]);
    let q_check = mul( mul(z1, x), z2 );
    assert_eq_quat(q, q_check);

    // XYX
    let q = from_euler_angles(Intrinsic, XYX, angles);
    let x1 = from_axis_angle([1.0, 0.0, 0.0], angles[0]);
    let y  = from_axis_angle([0.0, 1.0, 0.0], angles[1]);
    let x2 = from_axis_angle([1.0, 0.0, 0.0], angles[2]);
    let q_check = mul( mul(x1, y), x2 );
    assert_eq_quat(q, q_check);

    // YZY
    let q = from_euler_angles(Intrinsic, YZY, angles);
    let y1 = from_axis_angle([0.0, 1.0, 0.0], angles[0]);
    let z  = from_axis_angle([0.0, 0.0, 1.0], angles[1]);
    let y2 = from_axis_angle([0.0, 1.0, 0.0], angles[2]);
    let q_check = mul( mul(y1, z), y2 );
    assert_eq_quat(q, q_check);

    // ZYZ
    let q = from_euler_angles(Intrinsic, ZYZ, angles);
    let z1 = from_axis_angle([0.0, 0.0, 1.0], angles[0]);
    let y  = from_axis_angle([0.0, 1.0, 0.0], angles[1]);
    let z2 = from_axis_angle([0.0, 0.0, 1.0], angles[2]);
    let q_check = mul( mul(z1, y), z2 );
    assert_eq_quat(q, q_check);

    // XZX
    let q = from_euler_angles(Intrinsic, XZX, angles);
    let x1 = from_axis_angle([1.0, 0.0, 0.0], angles[0]);
    let z  = from_axis_angle([0.0, 0.0, 1.0], angles[1]);
    let x2 = from_axis_angle([1.0, 0.0, 0.0], angles[2]);
    let q_check = mul( mul(x1, z), x2 );
    assert_eq_quat(q, q_check);

    // YXY
    let q = from_euler_angles(Intrinsic, YXY, angles);
    let y1 = from_axis_angle([0.0, 1.0, 0.0], angles[0]);
    let x  = from_axis_angle([1.0, 0.0, 0.0], angles[1]);
    let y2 = from_axis_angle([0.0, 1.0, 0.0], angles[2]);
    let q_check = mul( mul(y1, x), y2 );
    assert_eq_quat(q, q_check);

    // ------ Tait–Bryan angles ------ //
    // XYZ
    let q = from_euler_angles(Intrinsic, XYZ, angles);
    let x = from_axis_angle([1.0, 0.0, 0.0], angles[0]);
    let y = from_axis_angle([0.0, 1.0, 0.0], angles[1]);
    let z = from_axis_angle([0.0, 0.0, 1.0], angles[2]);
    let q_check = mul( mul(x, y), z );
    assert_eq_quat(q, q_check);

    // YZX
    let q = from_euler_angles(Intrinsic, YZX, angles);
    let y = from_axis_angle([0.0, 1.0, 0.0], angles[0]);
    let z = from_axis_angle([0.0, 0.0, 1.0], angles[1]);
    let x = from_axis_angle([1.0, 0.0, 0.0], angles[2]);
    let q_check = mul( mul(y, z), x );
    assert_eq_quat(q, q_check);

    // ZXY
    let q = from_euler_angles(Intrinsic, ZXY, angles);
    let z = from_axis_angle([0.0, 0.0, 1.0], angles[0]);
    let x = from_axis_angle([1.0, 0.0, 0.0], angles[1]);
    let y = from_axis_angle([0.0, 1.0, 0.0], angles[2]);
    let q_check = mul( mul(z, x), y );
    assert_eq_quat(q, q_check);

    // XZY
    let q = from_euler_angles(Intrinsic, XZY, angles);
    let x = from_axis_angle([1.0, 0.0, 0.0], angles[0]);
    let z = from_axis_angle([0.0, 0.0, 1.0], angles[1]);
    let y = from_axis_angle([0.0, 1.0, 0.0], angles[2]);
    let q_check = mul( mul(x, z), y );
    assert_eq_quat(q, q_check);

    // ZYX
    let q = from_euler_angles(Intrinsic, ZYX, angles);
    let z = from_axis_angle([0.0, 0.0, 1.0], angles[0]);
    let y = from_axis_angle([0.0, 1.0, 0.0], angles[1]);
    let x = from_axis_angle([1.0, 0.0, 0.0], angles[2]);
    let q_check = mul( mul(z, y), x );
    assert_eq_quat(q, q_check);

    // YXZ
    let q = from_euler_angles(Intrinsic, YXZ, angles);
    let y = from_axis_angle([0.0, 1.0, 0.0], angles[0]);
    let x = from_axis_angle([1.0, 0.0, 0.0], angles[1]);
    let z = from_axis_angle([0.0, 0.0, 1.0], angles[2]);
    let q_check = mul( mul(y, x), z );
    assert_eq_quat(q, q_check);
}

#[test]
fn test_from_extrinsic_euler_angles() {
    use RotationType::*;
    use RotationSequence::*;

    let angles = [-PI/3.6, PI/4.5, -PI/2.5];

    // ------ Proper Euler angles ------ //
    // ZXZ
    let q = from_euler_angles(Extrinsic, ZXZ, angles);
    let z1 = from_axis_angle([0.0, 0.0, 1.0], angles[0]);
    let x  = from_axis_angle([1.0, 0.0, 0.0], angles[1]);
    let z2 = from_axis_angle([0.0, 0.0, 1.0], angles[2]);
    let q_check = mul( mul(z2, x), z1 );
    assert_eq_quat(q, q_check);

    // XYX
    let q = from_euler_angles(Extrinsic, XYX, angles);
    let x1 = from_axis_angle([1.0, 0.0, 0.0], angles[0]);
    let y  = from_axis_angle([0.0, 1.0, 0.0], angles[1]);
    let x2 = from_axis_angle([1.0, 0.0, 0.0], angles[2]);
    let q_check = mul( mul(x2, y), x1 );
    assert_eq_quat(q, q_check);

    // YZY
    let q = from_euler_angles(Extrinsic, YZY, angles);
    let y1 = from_axis_angle([0.0, 1.0, 0.0], angles[0]);
    let z  = from_axis_angle([0.0, 0.0, 1.0], angles[1]);
    let y2 = from_axis_angle([0.0, 1.0, 0.0], angles[2]);
    let q_check = mul( mul(y2, z), y1 );
    assert_eq_quat(q, q_check);

    // ZYZ
    let q = from_euler_angles(Extrinsic, ZYZ, angles);
    let z1 = from_axis_angle([0.0, 0.0, 1.0], angles[0]);
    let y  = from_axis_angle([0.0, 1.0, 0.0], angles[1]);
    let z2 = from_axis_angle([0.0, 0.0, 1.0], angles[2]);
    let q_check = mul( mul(z2, y), z1 );
    assert_eq_quat(q, q_check);

    // XZX
    let q = from_euler_angles(Extrinsic, XZX, angles);
    let x1 = from_axis_angle([1.0, 0.0, 0.0], angles[0]);
    let z  = from_axis_angle([0.0, 0.0, 1.0], angles[1]);
    let x2 = from_axis_angle([1.0, 0.0, 0.0], angles[2]);
    let q_check = mul( mul(x2, z), x1 );
    assert_eq_quat(q, q_check);

    // YXY
    let q = from_euler_angles(Extrinsic, YXY, angles);
    let y1 = from_axis_angle([0.0, 1.0, 0.0], angles[0]);
    let x  = from_axis_angle([1.0, 0.0, 0.0], angles[1]);
    let y2 = from_axis_angle([0.0, 1.0, 0.0], angles[2]);
    let q_check = mul( mul(y2, x), y1 );
    assert_eq_quat(q, q_check);

    // ------ Tait–Bryan angles ------ //
    // XYZ
    let q = from_euler_angles(Extrinsic, XYZ, angles);
    let x = from_axis_angle([1.0, 0.0, 0.0], angles[0]);
    let y = from_axis_angle([0.0, 1.0, 0.0], angles[1]);
    let z = from_axis_angle([0.0, 0.0, 1.0], angles[2]);
    let q_check = mul( mul(z, y), x );
    assert_eq_quat(q, q_check);

    // YZX
    let q = from_euler_angles(Extrinsic, YZX, angles);
    let y = from_axis_angle([0.0, 1.0, 0.0], angles[0]);
    let z = from_axis_angle([0.0, 0.0, 1.0], angles[1]);
    let x = from_axis_angle([1.0, 0.0, 0.0], angles[2]);
    let q_check = mul( mul(x, z), y );
    assert_eq_quat(q, q_check);

    // ZXY
    let q = from_euler_angles(Extrinsic, ZXY, angles);
    let z = from_axis_angle([0.0, 0.0, 1.0], angles[0]);
    let x = from_axis_angle([1.0, 0.0, 0.0], angles[1]);
    let y = from_axis_angle([0.0, 1.0, 0.0], angles[2]);
    let q_check = mul( mul(y, x), z );
    assert_eq_quat(q, q_check);

    // XZY
    let q = from_euler_angles(Extrinsic, XZY, angles);
    let x = from_axis_angle([1.0, 0.0, 0.0], angles[0]);
    let z = from_axis_angle([0.0, 0.0, 1.0], angles[1]);
    let y = from_axis_angle([0.0, 1.0, 0.0], angles[2]);
    let q_check = mul( mul(y, z), x );
    assert_eq_quat(q, q_check);

    // ZYX
    let q = from_euler_angles(Extrinsic, ZYX, angles);
    let z = from_axis_angle([0.0, 0.0, 1.0], angles[0]);
    let y = from_axis_angle([0.0, 1.0, 0.0], angles[1]);
    let x = from_axis_angle([1.0, 0.0, 0.0], angles[2]);
    let q_check = mul( mul(x, y), z );
    assert_eq_quat(q, q_check);

    // YXZ
    let q = from_euler_angles(Extrinsic, YXZ, angles);
    let y = from_axis_angle([0.0, 1.0, 0.0], angles[0]);
    let x = from_axis_angle([1.0, 0.0, 0.0], angles[1]);
    let z = from_axis_angle([0.0, 0.0, 1.0], angles[2]);
    let q_check = mul( mul(z, x), y );
    assert_eq_quat(q, q_check);
}

/// 四元数 ---> オイラー角の変換が正しく出来ているか確認する．
fn check_quat_to_euler(rt: RotationType, rs: RotationSequence, angles: Vector3<f64>) {
    let q = from_euler_angles(rt, rs, angles);
    let q2e = to_euler_angles(rt, rs, q);  // <-- ここで変換している
    println!("rt: {:?}, rs: {:?}", rt, rs);
    println!("angles: {:?}", angles);
    println!("q2e: {:?}", q2e);

    // 元のオイラー角と全く同じものを復元することは出来ないけど，回転としては同じものを表している．
    // 角度でチェックするよりベクトルの回転で見たほうが楽．
    let v = [1.0, 0.5, -0.2];
    let r_angles = point_rotation(q, v);
    let r_result = point_rotation(from_euler_angles(rt, rs, q2e), v);
    println!("roteta angles: {:?}", r_angles);
    println!("roteta q2e: {:?}", r_result);
    assert_eq_vec(r_angles, r_result);
    // 向きの異なるベクトルで再チェック
    let v = [-0.2, 1.5, 0.8];
    let r_angles = point_rotation(q, v);
    let r_result = point_rotation(from_euler_angles(rt, rs, q2e), v);
    println!("roteta angles: {:?}", r_angles);
    println!("roteta q2e: {:?}", r_result);
    assert_eq_vec(r_angles, r_result);
}

#[test]
fn test_to_intrinsic_euler_angles() {
    use RotationType::Intrinsic;
    use RotationSequence::*;

    let angles = [PI/1.6, PI/1.5, PI/2.5];

    // ------ Proper Euler angles ------ //
    check_quat_to_euler(Intrinsic, ZXZ, angles);
    check_quat_to_euler(Intrinsic, XYX, angles);
    check_quat_to_euler(Intrinsic, YZY, angles);
    check_quat_to_euler(Intrinsic, ZYZ, angles);
    check_quat_to_euler(Intrinsic, XZX, angles);
    check_quat_to_euler(Intrinsic, YXY, angles);

    // ------ Tait–Bryan angles ------ //
    check_quat_to_euler(Intrinsic, XYZ, angles);
    check_quat_to_euler(Intrinsic, YZX, angles);
    check_quat_to_euler(Intrinsic, ZXY, angles);
    check_quat_to_euler(Intrinsic, XZY, angles);
    check_quat_to_euler(Intrinsic, ZYX, angles);
    check_quat_to_euler(Intrinsic, YXZ, angles);
}

#[test]
fn test_to_extrinsic_euler_angles() {
    use RotationType::Extrinsic;
    use RotationSequence::*;

    let angles = [PI/0.6, PI, PI/2.5];

    // ------ Proper Euler angles ------ //
    check_quat_to_euler(Extrinsic, ZXZ, angles);
    check_quat_to_euler(Extrinsic, XYX, angles);
    check_quat_to_euler(Extrinsic, YZY, angles);
    check_quat_to_euler(Extrinsic, ZYZ, angles);
    check_quat_to_euler(Extrinsic, XZX, angles);
    check_quat_to_euler(Extrinsic, YXY, angles);

    // ------ Tait–Bryan angles ------ //
    check_quat_to_euler(Extrinsic, XYZ, angles);
    check_quat_to_euler(Extrinsic, YZX, angles);
    check_quat_to_euler(Extrinsic, ZXY, angles);
    check_quat_to_euler(Extrinsic, XZY, angles);
    check_quat_to_euler(Extrinsic, ZYX, angles);
    check_quat_to_euler(Extrinsic, YXZ, angles);
}

#[test]
fn test_rotation_vector() {
    // 適当なVersorを作る
    let q = from_axis_angle([1.0, 2.0, 3.0], PI);

    // 回転ベクトルに変換
    let rot_v = to_rotation_vector(q);

    // 回転ベクトルから復元したVecsor
    let mut q_rest = from_rotation_vector(rot_v);

    // 符号が反転していても，３次元空間上で表す回転は同じ
    if q.0.is_sign_positive() && q_rest.0.is_sign_negative() {
        q_rest = negate(q_rest);
    }

    assert_eq_quat(q, q_rest);
}

#[test]
fn test_mul() {
    let a = (1.0f64, [0.5, -1.2, 3.0]);
    let b = (0.5f64, [-0.2, 2.5, -3.3]);

    let v0 = scale(a.0, b.1);
    let v1 = scale(b.0, a.1);
    let p = (
        a.0 * b.0 - dot(a.1, b.1),
        add( add(v0, v1), cross(a.1, b.1) )
    );

    let q = mul(a, b);
    assert_eq_quat(q, p);
}

// 手計算した結果で動作確認
#[test]
fn test_point_rotation() {
    let r: Vector3<f64> = [2.0, 2.0, 0.0];
    let q = from_axis_angle([0.0, 1.0, 0.0], PI/2.0);
    let result = point_rotation(q, r);

    assert_eq_vec(result, [0.0, 2.0, -2.0]);
}

// 手計算した結果で動作確認
#[test]
fn test_frame_rotation() {
    let r = [2.0, 2.0, 0.0];
    let q = from_axis_angle([0.0, 1.0, 0.0], PI/2.0);
    let result = frame_rotation(q, r);

    assert_eq_vec(result, [0.0, 2.0, 2.0]);
}

// 回転軸と回転角を取り出す
#[test]
fn test_to_axis_angle() {
    let axis = [1.0, 4.0, 2.0];
    let angle = PI/4.0;
    let q = from_axis_angle(axis, angle);

    // 軸の方向は分かるが，元の大きさはわからない．
    let n = normalize(axis);
    let f = to_axis_angle(q);
    assert!( (f.1 - angle).abs() < EPSILON );
    for i in 0..3 {
        assert!( (f.0[i] - n[i]).abs() < EPSILON );
    }
}

#[test]
fn test_rotate_a_to_b() {
    let a = normalize([1.0f64, -0.5, -2.0]);

    // 適当にbを作って確認
    {
        let q = from_axis_angle([-1.0, 5.0, -0.5], 0.7);
        let b = point_rotation(q, a);
    
        let a_to_b = rotate_a_to_b(a, b).unwrap();
        let a_to_b_t = rotate_a_to_b_shortest(a, b, 1.0).unwrap();
        let b_rest = point_rotation(a_to_b, a);
        let b_rest_t = point_rotation(a_to_b_t, a);
        assert_eq_vec(b, b_rest);
        assert_eq_vec(b, b_rest_t);
        assert!((1.0 - norm(a_to_b)).abs() < EPSILON);
        assert!((1.0 - norm(a_to_b_t)).abs() < EPSILON);
    }

    // aとbが平行かつ逆向きな場合（回転軸の求め方が特殊になる）
    {
        let b = negate(a);
        let a_to_b = rotate_a_to_b(a, b).unwrap();
        let a_to_b_shortest = rotate_a_to_b_shortest(a, b, 1.0).unwrap();
        let b_rest = point_rotation(a_to_b, a);
        let b_rest_param = point_rotation(a_to_b_shortest, a);
        assert_eq_vec(b, b_rest);
        assert_eq_vec(b, b_rest_param);
        assert!((1.0 - norm(a_to_b)).abs() < EPSILON);
        assert!((1.0 - norm(a_to_b_shortest)).abs() < EPSILON);
    }

    // ゼロベクトルを入れた場合
    {
        let b = [0.0; 3];
        let a_to_b = rotate_a_to_b(a, b);
        assert_eq!(None, a_to_b);
        let a_to_b_shortest = rotate_a_to_b_shortest(a, b, 1.0);
        assert_eq!(None, a_to_b_shortest);
    }
}

fn assert_eq_vec(a: Vector3<f64>, b: Vector3<f64>) {
    for i in 0..a.len() {
        assert!((a[i] - b[i]).abs() < EPSILON);
    }
}

/// 二つの四元数を比較する．
/// 
/// 表す回転が等しければ良いので符号が反転していても気にしない．
fn assert_eq_quat(a: Quaternion<f64>, b: Quaternion<f64>) {
    if dot(a, b).is_sign_positive() {
        assert!((a.0 - b.0).abs() < EPSILON);
        assert_eq_vec(a.1, b.1);
    } else {
        // bの符号を反転
        assert!((a.0 + b.0).abs() < EPSILON);
        assert_eq_vec(a.1, negate(b.1));
    }
}
