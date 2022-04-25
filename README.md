# Quaternion

## About

### 【日本語】

Rustで作成した四元数（クォータニオン）計算用のライブラリです。

`version=2.3.0`以降では`no_std`環境で使用できます。

関数の動作については、`src/lib.rs`内のドキュメンテーションコメントを参照してください。
もしくは、`cargo doc`コマンドから確認することもできます。

また、四元数自体について詳しく知りたい方は以下の資料をご覧ください。この資料中に出てくる演算はすべて実装してあります。

- [四元数まとめ資料（宇宙電波実験室）](https://space-denpa.jp/2019/03/26/quaternion-doc/)

### 【English】

It's a Quaternion library written in Rust.

In `version=2.3.0` or later, it is available in the `no_std` environment.

For function behavior, see the documentation comments in `src/lib.rs`.

## FMA (Fused Multiply-Add)

In this library, we use the [mul_add](https://doc.rust-lang.org/std/primitive.f64.html#method.mul_add) method mainly to improve the performance, but by default it is replace with a unfused multiply-add (`s*a + b`) . If you wish to use mul_add method, enable the `fma` feature in the `Cargo.toml` file.

If your CPU does not support FMA instructions, or if you use libm (running in no_std environment), enabling the `fma` feature may cause slowdown of execution speed. Also, due to rounding error, `s.mul_add(a, b)` and `s*a + b` will not match perfectly.

## Dependencies

- [num-traits](https://crates.io/crates/num-traits)

- [libm](https://crates.io/crates/libm) (optional)

## License

- [MIT license](https://opensource.org/licenses/MIT)

## Example of use

Cargo.toml

```toml
[dependencies.quaternion]
git = "https://github.com/HamaguRe/quaternion.git"
version = "2.6"

# Uncomment if you wish to use "mul_add" method.
#features = ["fma"]

# Uncomment if you wish to use in "no_std" environment.
#default-features = false
#features = ["libm"]
```

src/main.rs

```rust
use quaternion as quat;
use quat::Vector3;

const PI: f64 = std::f64::consts::PI;
const EPSILON: f64 = 1e-14;

fn main() {
    // Position vector
    let r: Vector3<f64> = [2.0, 2.0, 0.0];

    // Generates a quaternion representing the
    // rotation of π/2[rad] around the y-axis.
    let q = quat::from_axis_angle([0.0, 1.0, 0.0], PI/2.0);

    let result = quat::vector_rotation(q, r);

    // Check if the calculation is correct.
    let diff = quat::sub_vec(result, [0.0, 2.0, -2.0]);
    for val in diff.iter() {
        assert!( val.abs() < EPSILON );
    }
}
```
