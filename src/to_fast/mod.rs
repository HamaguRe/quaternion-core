// SIMDで高速化を目指す
// 
// アトリビュートでの判定にバグがあるらしく，
// 使える命令をコンパイル時に教えてあげる必要がある．
// avxとfmaの両方に対応している場合は，
// $ RUSTFLAGS='-C target-cpu=native' cargo build
// もしくは
// $ RUSTFLAGS='-C target-feature=+avx,fma' cargo build
// とする．avxのみ対応している場合はこっち↓
// $ RUSTFLAGS='-C target-feature=+avx' cargo build

#[allow(unused_imports)]
use super::Quaternion;

#[allow(unused_imports)]
#[cfg(target_arch = "x86")]
use std::arch::x86::*;

#[allow(unused_imports)]
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;


#[inline(always)]
#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "avx"))]
pub fn dot(a: Quaternion<f64>, b: Quaternion<f64>) -> f64 {
    let c: [f64; 4] = unsafe {
        let a = _mm256_setr_pd(a.0, (a.1)[0], (a.1)[1], (a.1)[2]);
        let b = _mm256_setr_pd(b.0, (b.1)[0], (b.1)[1], (b.1)[2]);
        std::mem::transmute( _mm256_mul_pd(a, b) )
    };
    c[0] + c[1] + c[2] + c[3]
}

#[inline(always)]
#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "avx"))]
pub fn add(a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
    unsafe {
        let a = _mm256_setr_pd(a.0, (a.1)[0], (a.1)[1], (a.1)[2]);
        let b = _mm256_setr_pd(b.0, (b.1)[0], (b.1)[1], (b.1)[2]);
        std::mem::transmute( _mm256_add_pd(a, b) )
    }
}

#[inline(always)]
#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "avx"))]
pub fn sub(a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
    unsafe {
        let a = _mm256_setr_pd(a.0, (a.1)[0], (a.1)[1], (a.1)[2]);
        let b = _mm256_setr_pd(b.0, (b.1)[0], (b.1)[1], (b.1)[2]);
        std::mem::transmute( _mm256_sub_pd(a, b) )
    }
}

#[inline(always)]
#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "avx"))]
pub fn scale(s: f64, q: Quaternion<f64>) -> Quaternion<f64> {
    unsafe {
        let s = _mm256_set1_pd(s);
        let q = _mm256_setr_pd(q.0, (q.1)[0], (q.1)[1], (q.1)[2]);
        std::mem::transmute( _mm256_mul_pd(s, q) )
    }
}

#[inline(always)]
#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "fma"))]
pub fn scale_add(s: f64, a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
    unsafe {
        let s = _mm256_set1_pd(s);
        let a = _mm256_setr_pd(a.0, (a.1)[0], (a.1)[1], (a.1)[2]);
        let b = _mm256_setr_pd(b.0, (b.1)[0], (b.1)[1], (b.1)[2]);
        std::mem::transmute( _mm256_fmadd_pd(s, a, b) )
    }
}