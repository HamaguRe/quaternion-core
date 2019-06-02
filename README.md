# Quaternion library
#### 【日本語】
  Rustで作成したクォータニオン計算用のライブラリです。

  関数の動作については、ドキュメンテーションコメントを参照してください。

  また、クォータニオン自体について知りたい方は、以下のページを参照してください。

  * [四元数まとめ資料を書いた（宇宙電波実験室）](https://space-denpa.jp/2019/03/26/quaternion-doc/)

#### 【English】
  It's a quaternion library written in Rust.

  Refer to the documentation comments for function behavior.

# Example of use
#### Cargo.toml
```
[dependencies]
quaternion = {git="https://github.com/HamaguRe/quaternion.git", rev="d0a504af8068a60c6b7ad2d8ae875eab7ba8e9e8"}
```

#### quaternion_example.rs
```
extern crate quaternion;
use quaternion as quat;

const PI: f64 = std::f64::consts::PI;
const EPSILON: f64 = 1e-8;

fn main() {
    // Position vector
    let r = [2.0, 2.0, 0.0];
    // π/2[rad] rotation around the y axis.
    let q = quat::from_axis_angle([0.0, 1.0, 0.0], PI/2.0);

    let rotated = quat::vector_rotation(q, r);
    assert!( (rotated[0] - 0.0).abs() < EPSILON );
    assert!( (rotated[1] - 2.0).abs() < EPSILON );
    assert!( (rotated[2] + 2.0).abs() < EPSILON );
}
```