use quaternion::*;

const PI: f64 = std::f64::consts::PI;
const EPSILON: f64 = 1e-14;  // libmを使う場合は1e-12に落とさないと通らない


// 二つの異なる方法でVersorの軸回りの回転角を求める．
// 正の回転角を指定してVersorを生成した場合は以下の方法で回転角と回転軸を復元できるが，
// 負の回転角を指定した場合には，回転角は正となり回転軸が反転するので完全には復元できない．
// とはいえ，回転軸と回転角を合わせて考えれば三次元空間上での回転は変化していないので
// そこまで大きな問題はないと思われる．
#[test]
fn test_get_angle() {
    // 適当なVersorを作る
    let axis = [0.0, 0.0, 2.0];
    let angle = -1.5 * PI;
    let q = from_axis_angle(axis, angle);
    assert!( ( norm(q) - 1.0 ).abs() < EPSILON );

    // 方法1
    // 実部のみを使うので計算量が少ない．
    // 計算精度が一つの変数に依存してしまうのは良くない...？
    // 0 <= angle1 <= 2π
    let angle1 = 2.0 * q.0.acos();
    // 方法2
    // 実部の符号を反映することで幾何学的には方法1と同じ結果が得られる．
    // 実部と虚部両方の値を使っているのでなんとなく気持ちが良い．
    // -π <= angle2 <= π
    let angle2 = ( 2.0 * norm_vec(q.1).asin() ).copysign(q.0);

    println!("axis: {:?}", normalize_vec(q.1));
    println!("angle1: {}PI, angle2: {}PI", angle1/PI, angle2/PI);

    assert!( (angle1 - 1.5*PI).abs() < EPSILON );
    assert!( (angle2 + 0.5*PI).abs() < EPSILON );
}

#[test]
fn test_axis_angle() {
    // to_axis_angle関数はどのような回転角を入れても正常に変換できるが，
    // テストコードの実装上，angleの範囲は(-2π, 2π)とする．
    let axis = normalize_vec([0.0, 1.0, 1.0]);
    let angle = -1.5*PI;
    let q = from_axis_angle(axis, angle);

    assert!( (1.0 - norm(q)).abs() < EPSILON );

    let (re_axis, re_angle) = to_axis_angle(q);

    println!("re_axis: {:?}", re_axis);
    println!("re_angle: {}*PI", re_angle / PI);

    // 軸ベクトルのチェック
    if angle.is_sign_positive() {
        for i in 0..3 {
            assert!( (re_axis[i] - axis[i]).abs() < EPSILON );
        }
    } else {
        // 負の回転角の場合には回転軸が反転する
        for i in 0..3 {
            assert!( (re_axis[i] + axis[i]).abs() < EPSILON );
        }
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
    let v = [2.5, 1.0, -3.0];
    let diff = [0.2, -0.1, 2.5];
    let mut axis = [1.0, -0.2, 0.9];
    for i in 0..20 {
        axis = add_vec(axis, diff);
        let q = from_axis_angle(axis, PI * (i as f64));  // Versor
        let dcm = to_dcm(q);
        let q_rest = from_dcm(dcm);

        assert!( ( norm(q_rest) - 1.0 ).abs() < EPSILON );
    
        let rotated_q = vector_rotation(q, v);
        let rotated_q_rest = vector_rotation(q_rest, v);
    
        assert!( (rotated_q[0] - rotated_q_rest[0]).abs() < EPSILON );
        assert!( (rotated_q[1] - rotated_q_rest[1]).abs() < EPSILON );
        assert!( (rotated_q[2] - rotated_q_rest[2]).abs() < EPSILON );
    }
}

#[test]
fn test_euler() {
    let q_z = from_axis_angle([0.0, 0.0, 1.0f64], PI*(3.0/4.0));
    let q_y = from_axis_angle([0.0, 1.0, 0.0f64], PI*(1.0/4.0));
    let q_x = from_axis_angle([1.0, 0.0, 0.0f64], PI/2.0);

    let q = mul(q_x, mul(q_y, q_z));
    let [yaw, pitch, roll] = to_euler_angle(q);

    let tmp = 180.0 / PI;
    println!("yaw: {}, pitch: {}, roll: {}", yaw*tmp, pitch*tmp, roll*tmp);

    assert!( (yaw + (PI/4.0)).abs() < EPSILON );
    assert!( (pitch - (PI/4.0)).abs() < EPSILON );
    assert!( (roll + (PI/2.0)).abs() < EPSILON );
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

    let diff = sub(q, q_rest);
    assert!(diff.0.abs() < EPSILON);
    for i in 0..3 {
        assert!(diff.1[i].abs() < EPSILON);
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
    let s = 2.0_f64;
    let a = (0.5, [1.0, 2.5, 1.2]);
    let b = (2.2, [6.5, 1.0, 3.4]);

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
fn test_rotate_by_dcm() {
    let r: Vector3<f64> = [2.0, 2.0, 0.0];
    let q = from_axis_angle([0.0, 1.0, 0.0], PI/2.0);

    // 位置ベクトルの回転
    let m = to_dcm(q);
    let result = matrix_product(m, r);
    let diff = sub_vec(result, [0.0, 2.0, -2.0]);
    for i in 0..3 {
        assert!( diff[i].abs() < EPSILON );
    }

    // 座標系の回転
    let m = to_dcm( conj(q) );
    let result = matrix_product(m, r);
    let diff = sub_vec(result, [0.0, 2.0, 2.0]);
    for i in 0..3 {
        assert!( diff[i].abs() < EPSILON );
    }
}

// 手計算した結果で動作確認
#[test]
fn test_vector_rotation() {
    let r: Vector3<f64> = [2.0, 2.0, 0.0];
    let q = from_axis_angle([0.0, 1.0, 0.0], PI/2.0);
    let result = vector_rotation(q, r);

    let diff = sub_vec(result, [0.0, 2.0, -2.0]);
    for i in 0..3 {
        assert!( diff[i].abs() < EPSILON );
    }
}

#[test]
fn test_frame_rotation() {
    let r = [2.0, 2.0, 0.0];
    let q = from_axis_angle([0.0, 1.0, 0.0], PI/2.0);
    let result = frame_rotation(q, r);

    let diff = sub_vec(result, [0.0, 2.0, 2.0]);
    for i in 0..3 {
        assert!( diff[i].abs() < EPSILON );
    }
}

// 回転軸と回転角を取り出す
#[test]
fn test_to_axis_angle() {
    let axis = [1.0, 4.0, 2.0];
    let angle = PI/4.0;
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
fn test_rotate_a_to_b() {
    let a = normalize_vec([1.0f64, -0.5, -2.0]);
    
    let q = from_axis_angle([1.0, 0.5, -0.5], PI);
    let b = vector_rotation(q, a);

    let a_to_b = rotate_a_to_b(a, b, 1.0);
    let b_rest = vector_rotation(a_to_b, a);

    let diff = sub_vec(b, b_rest);
    for i in 0..3 {
        assert!( diff[i].abs() < EPSILON );
    }
}