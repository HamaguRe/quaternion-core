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

#[test]
fn test_dcm() {
    let v = [2.5, 1.0, -3.0];
    let diff = [0.2, -0.1, 2.5];
    let mut axis = [1.0, -0.2, 0.9];
    for i in 0..20 {
        axis = add_vec(axis, diff);
        let q = from_axis_angle(axis, PI * (i as f64));
        let dcm = to_dcm_vector(q);
        let q_rest = from_dcm_vector(dcm);
    
        let rotated_q = vector_rotation(q, v);
        let rotated_q_rest = vector_rotation(q_rest, v);
    
        assert!( (rotated_q[0] - rotated_q_rest[0]).abs() < EPSILON );
        assert!( (rotated_q[1] - rotated_q_rest[1]).abs() < EPSILON );
        assert!( (rotated_q[2] - rotated_q_rest[2]).abs() < EPSILON );
    }
}

// 回転行列に変換してからベクトルを回転させる．
#[test]
fn test_rotate_by_dcm() {
    let r: Vector3<f64> = [2.0, 2.0, 0.0];
    let q = from_axis_angle([0.0, 1.0, 0.0], PI/2.0);

    // 位置ベクトルの回転
    let m = to_dcm_vector(q);
    let result = matrix_product(m, r);
    assert!( (result[0] - 0.0).abs() < EPSILON);
    assert!( (result[1] - 2.0).abs() < EPSILON);
    assert!( (result[2] + 2.0).abs() < EPSILON);

    // 座標系の回転
    let m = to_dcm_frame(q);
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
    //let angle2 = 2.0 * norm_vec(q.1).asin();

    println!("axis: {:?}", normalize_vec(q.1));
    println!("angle1: {}PI, angle2: {}PI", angle1/PI, angle2/PI);

    assert!( (angle1 - 1.5*PI).abs() < EPSILON );
    assert!( (angle2 + 0.5*PI).abs() < EPSILON );
}

#[test]
fn test_cross() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];

    let r = cross_vec(r1, r2);
    assert_eq!( [0.0, 0.0, 1.0] , r );
}