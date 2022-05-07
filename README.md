# quaternion-core

`Quaternion` library written in Rust.

## About

### Japanese

四元数に関する基本的な演算や各種回転表現との変換を関数として提供します。

### English

It provides basic operations on quaternions and conversions to and from several attitude representations as functions.

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
quaternion-core = "0.1"
```

## Features

Cargo.toml

```toml
[dependencies.quaternion-core]
version = "0.1"

# Uncomment if you wish to use FMA and SIMD.
#features = ["fma", "simd"]

# Uncomment if you wish to use in "no_std" environment.
#default-features = false
#features = ["libm"]
```

### The `fma` feature

This library uses the [mul_add](https://doc.rust-lang.org/std/primitive.f64.html#method.mul_add) method mainly to improve the performance, but by default it is replace with a unfused multiply-add (`s*a + b`) . If you wish to use mul_add method, enable the `fma` feature.

If your CPU does not support FMA instructions, or if you use `libm` (running in no_std environment), enabling the `fma` feature may cause slowdown of computation speed. Also, due to rounding error, results of `s.mul_add(a, b)` and `s*a + b` will not match perfectly.

### The `simd` feature

___Attension!! : This feature may have bugs and should not be enabled at first.___

By enabling this feature, the SIMD implementation using the [std::arch](https://docs.rs/rustc-std-workspace-std/1.0.1/std/arch/index.html) module can be used in some functions.

Currently (`version="0.1.0"`) only `x86` and `x86_64` architectures are supported.

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

## Example

src/main.rs

```rust
use quaternion_core as quat;

const PI: f64 = std::f64::consts::PI;
const EPSILON: f64 = 1e-14;

fn main() {
    // Position vector
    let r = [2.0, 2.0, 0.0];

    // Generates a quaternion representing the
    // rotation of π/2[rad] around the y-axis.
    let q = quat::from_axis_angle([0.0, 1.0, 0.0], PI/2.0);

    let result = quat::vector_rotation(q, r);

    // Check if the calculation is correct.
    let true_val = [0.0, 2.0, -2.0];
    let diff = quat::sub_vec(true_val, result);
    for val in diff.iter() {
        assert!( val.abs() < EPSILON );
    }
}
```

## License

Licensed under either of
[Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0)
or
[MIT License](https://opensource.org/licenses/MIT)
at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in the work by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any additional terms or conditions.
