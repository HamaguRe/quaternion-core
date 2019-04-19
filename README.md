# Quaternion library
#### 【日本語】
  Rustで作成したクォータニオン計算用のライブラリです。
  
  クォータニオンの勉強を兼ねて作成したものなので、計算結果の正当性は保証できません。

  関数の動作については、ドキュメンテーションコメントを参照してください。

#### 【English】
  It is a quaternion library written in Rust.
  
  Because it was created as a study of quaternion, I can not guarantee the calculation result.

  Refer to the documentation comments for function behavior.

# Examle of use
#### Cargo.toml
```
[dependencies]
quaternion = {git="https://github.com/HamaguRe/quaternion.git", rev="d79af72367b644357c6f3dd481123c221778aa2c"}
```

#### example.rs
```
extern crate quaternion;
use quaternion as quat;

const PI: f64 = std::f64::consts::PI;
const EPSILON: f64 = 0.00000001;

fn main() {
    // vector rotation
    let r = [2.0, 2.0, 0.0];
    // π/2[rad] rotation around the y axis.
    let q = quat::from_axis_angle([0.0, 1.0, 0.0], PI/2.0);

    let rotated = quat::vector_rotation(q, r);
    assert!( (rotated[0] - 0.0).abs() < EPSILON );
    assert!( (rotated[1] - 2.0).abs() < EPSILON );
    assert!( (rotated[2] + 2.0).abs() < EPSILON );
}
```