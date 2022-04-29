# Quaternion

## About

### 【日本語】

Rustで作成した四元数（クォータニオン）計算用のライブラリです。

`version="2.3.0"`以降では`no_std`環境で使用できます。

関数の動作については、`src/lib.rs`内のドキュメンテーションコメントを参照してください。
もしくは、`cargo doc`コマンドから確認することもできます。

また、四元数自体について詳しく知りたい方は以下の資料をご覧ください。この資料中に出てくる演算はすべて実装してあります。

- [四元数まとめ資料（宇宙電波実験室）](https://space-denpa.jp/2019/03/26/quaternion-doc/)

### 【English】

It's a Quaternion library written in Rust.

In `version="2.3.0"` or later, it is possible to use in the `no_std` environment.

For function behavior, see the documentation comments in `src/lib.rs`. Alternatively, you can check it from the `cargo doc` command.

## Features

### The `fma` feature

This library uses the [mul_add](https://doc.rust-lang.org/std/primitive.f64.html#method.mul_add) method mainly to improve the performance, but by default it is replace with a unfused multiply-add (`s*a + b`) . If you wish to use mul_add method, enable the `fma` feature.

If your CPU does not support FMA instructions, or if you use `libm` (running in no_std environment), enabling the `fma` feature may cause slowdown of computation speed. Also, due to rounding error, results of `s.mul_add(a, b)` and `s*a + b` will not match perfectly.

### The `simd` feature

By enabling this feature, the SIMD implementation using the [std::arch](https://docs.rs/rustc-std-workspace-std/1.0.1/std/arch/index.html) module can be used in some functions.

Currently (`version="2.7.0"`) only `x86` and `x86_64` architectures are supported.

To enable this feature, CPU must support these instruction sets:
```
SSE, SSE2, SSE3, SSE4.1, AVX, FMA
```

Also, specify the `-C target-cpu` flag to the compiler as follows:

```console
$ RUSTFLAGS='-C target-cpu=native' cargo build
```

### The `libm` feature and `default-feature = false`

These options allow for use in the `no_std` environment. In this case, mathematical functions (e.g. sin, cos, sqrt ...) are provided by `libm`.

## Dependencies

- [num-traits](https://crates.io/crates/num-traits)

- [libm](https://crates.io/crates/libm) (optional)

## License

Licensed under either of
[Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0)
or
[MIT License](https://opensource.org/licenses/MIT)
at your option.

Since `version="2.7.0"`, the Apache License was added and became dual license.

## Example of use

Cargo.toml

```toml
[dependencies.quaternion]
git = "https://github.com/HamaguRe/quaternion.git"
version = "2.7"

# Uncomment if you wish to use FMA and SIMD.
#features = ["fma", "simd"]

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
