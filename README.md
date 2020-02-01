# Quaternion library
## 【日本語】
  Rustで作成した四元数（クォータニオン）計算用のライブラリです。

  実装の都合上、計算精度はf64のみとなっています。

  関数の動作については、ドキュメンテーションコメントを参照してください。

  また、四元数自体について詳しく知りたい方は、以下のページをご覧ください。

  * [四元数まとめ資料を書いた（宇宙電波実験室）](https://space-denpa.jp/2019/03/26/quaternion-doc/)

### 例外処理に関して
 零ベクトルを入力した場合には例外処理を行いますが、すべての要素が零となる四元数が入力される状況（ベクトルを純虚四元数として扱った場合などに起こる）は考慮していません。

## 【English】
  It's a Quaternion library written in Rust.

  The calculation accuracy is only f64 for implementation reasons.

  Refer to the documentation comments for function behavior.

### About exception handling
When a zero vector is input, exception processing is performed. However, a situation where a quaternion in which all elements are zero (which occurs when a vector is treated as a pure quaternion) is not considered.

# Example of use
#### Cargo.toml
```
[dependencies]
quaternion = {git="https://github.com/HamaguRe/quaternion.git"}
```

#### src/main.rs
```
use quaternion as quat;
use quat::Vector3;

const PI: f64 = std::f64::consts::PI;
const EPSILON: f64 = 1e-10;

fn main() {
    // Position vector
    let r: Vector3<f64> = [2.0, 2.0, 0.0];

    // π/2[rad] rotation around the y axis.
    let q = quat::from_axis_angle([0.0, 1.0, 0.0], PI/2.0);

    let rotated = quat::vector_rotation(q, r);
    assert!( (rotated[0] - 0.0).abs() < EPSILON );
    assert!( (rotated[1] - 2.0).abs() < EPSILON );
    assert!( (rotated[2] + 2.0).abs() < EPSILON );
}
```