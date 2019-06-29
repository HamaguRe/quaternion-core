# Quaternion library
#### 【日本語】
  Rustで作成した四元数（クォータニオン）計算用のライブラリです。

  関数の動作については、ドキュメンテーションコメントを参照してください。

  また、四元数自体について詳しく知りたい方は、以下のページをご覧ください。

  * [四元数まとめ資料を書いた（宇宙電波実験室）](https://space-denpa.jp/2019/03/26/quaternion-doc/)

#### 【English】
  It's a quaternion library written in Rust.

  Refer to the documentation comments for function behavior.

# Example of use
#### Cargo.toml
```
[dependencies]
quaternion = {git="https://github.com/HamaguRe/quaternion.git", rev="f0100da10378f36b7079b0452485460577e0ae04"}
```

#### src/main.rs
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