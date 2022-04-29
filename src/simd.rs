//! SIMD module
//! 
//! Support architectures: x86, x86_64

use super::Quaternion;

#[cfg(not(all(feature = "std", feature = "simd")))]
use super::Float;

#[cfg(all(feature = "std", feature = "simd"))]
use std::mem;

#[cfg(all(feature = "std", feature = "simd", target_arch = "x86"))]
use std::arch::x86::*;
#[cfg(all(feature = "std", feature = "simd", target_arch = "x86_64"))]
use std::arch::x86_64::*;

pub trait FloatSimd<T> {
    fn sum(q: Quaternion<T>) -> T;
    fn add(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>;
    fn sub(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>;
    fn scale(s: T, q: Quaternion<T>) -> Quaternion<T>;
    fn scale_add(s: T, a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>;
    fn hadamard(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T>;
    fn hadamard_add(a: Quaternion<T>, b: Quaternion<T>, c: Quaternion<T>) -> Quaternion<T>;
    fn dot(a: Quaternion<T>, b: Quaternion<T>) -> T;
    fn negate(q: Quaternion<T>) -> Quaternion<T>;
}

#[cfg(all(feature = "std", feature = "simd"))]
impl<T> FloatSimd<f32> for T {
    #[inline]
    fn sum(q: Quaternion<f32>) -> f32 {
        unsafe {
            let a = _mm_set_ps(q.0, q.1[0], q.1[1], q.1[2]);  // SSE
            
            let mut c: [f32; 4] = mem::MaybeUninit::uninit().assume_init();
            _mm_storeu_ps(&mut c[0], _mm_hadd_ps(a, a));  // SSE3
            c[0] + c[1]
        }
    }

    #[inline]
    fn add(a: Quaternion<f32>, b: Quaternion<f32>) -> Quaternion<f32> {
        unsafe {
            let a = _mm_setr_ps(a.0, a.1[0], a.1[1], a.1[2]);  // SSE
            let b = _mm_setr_ps(b.0, b.1[0], b.1[1], b.1[2]);
            
            let mut c: [f32; 4] = mem::MaybeUninit::uninit().assume_init();
            _mm_storeu_ps(&mut c[0], _mm_add_ps(a, b));  // SSE
            mem::transmute::<[f32; 4], Quaternion<f32>>(c)
        }
    }

    #[inline]
    fn sub(a: Quaternion<f32>, b: Quaternion<f32>) -> Quaternion<f32> {
        unsafe {
            let a = _mm_setr_ps(a.0, a.1[0], a.1[1], a.1[2]);  // SSE
            let b = _mm_setr_ps(b.0, b.1[0], b.1[1], b.1[2]);
            
            let mut c: [f32; 4] = mem::MaybeUninit::uninit().assume_init();
            _mm_storeu_ps(&mut c[0], _mm_sub_ps(a, b));  // SSE
            mem::transmute::<[f32; 4], Quaternion<f32>>(c)
        }
    }

    #[inline]
    fn scale(s: f32, q: Quaternion<f32>) -> Quaternion<f32> {
        unsafe {
            let s = _mm_set1_ps(s);  // SSE
            let q = _mm_setr_ps(q.0, q.1[0], q.1[1], q.1[2]);  // SSE
            
            let mut c: [f32; 4] = mem::MaybeUninit::uninit().assume_init();
            _mm_storeu_ps(&mut c[0], _mm_mul_ps(s, q));  // SSE
            mem::transmute::<[f32; 4], Quaternion<f32>>(c)
        }
    }

    #[inline]
    fn scale_add(s: f32, a: Quaternion<f32>, b: Quaternion<f32>) -> Quaternion<f32> {
        unsafe {
            let s = _mm_set1_ps(s);  // SSE
            let a = _mm_setr_ps(a.0, a.1[0], a.1[1], a.1[2]);  // SSE
            let b = _mm_setr_ps(b.0, b.1[0], b.1[1], b.1[2]);
            
            let mut c: [f32; 4] = mem::MaybeUninit::uninit().assume_init();
            if cfg!(feature = "fma") {
                _mm_storeu_ps(&mut c[0], _mm_fmadd_ps(s, a, b));  // SSE, FMA
            } else {
                _mm_storeu_ps(&mut c[0], _mm_add_ps(_mm_mul_ps(s, a), b));  // SSE
            }
            mem::transmute::<[f32; 4], Quaternion<f32>>(c)
        }
    }

    #[inline]
    fn hadamard(a: Quaternion<f32>, b: Quaternion<f32>) -> Quaternion<f32> {
        unsafe {
            let a = _mm_setr_ps(a.0, a.1[0], a.1[1], a.1[2]);  // SSE
            let b = _mm_setr_ps(b.0, b.1[0], b.1[1], b.1[2]);
            
            let mut c: [f32; 4] = mem::MaybeUninit::uninit().assume_init();
            _mm_storeu_ps(&mut c[0], _mm_mul_ps(a, b));  // SSE
            mem::transmute::<[f32; 4], Quaternion<f32>>(c)
        }
    }

    #[inline]
    fn hadamard_add(a: Quaternion<f32>, b: Quaternion<f32>, c: Quaternion<f32>) -> Quaternion<f32> {
        unsafe {
            let a = _mm_setr_ps(a.0, a.1[0], a.1[1], a.1[2]);  // SSE
            let b = _mm_setr_ps(b.0, b.1[0], b.1[1], b.1[2]);
            let c = _mm_setr_ps(c.0, c.1[0], c.1[1], c.1[2]);
            
            let mut d: [f32; 4] = mem::MaybeUninit::uninit().assume_init();
            if cfg!(feature = "fma") {
                _mm_storeu_ps(&mut d[0], _mm_fmadd_ps(a, b, c));  // SSE, FMA
            } else {
                _mm_storeu_ps(&mut d[0], _mm_add_ps(_mm_mul_ps(a, b), c));  // SSE
            }
            mem::transmute::<[f32; 4], Quaternion<f32>>(d)
        }
    }

    #[inline]
    fn dot(a: Quaternion<f32>, b: Quaternion<f32>) -> f32 {
        unsafe {
            let a = _mm_setr_ps(a.0, a.1[0], a.1[1], a.1[2]);  // SSE
            let b = _mm_setr_ps(b.0, b.1[0], b.1[1], b.1[2]);

            let mut c: f32 = mem::MaybeUninit::uninit().assume_init();
            _mm_store_ss(&mut c, _mm_dp_ps(a, b, 0b11110001));  // SSE, SSE4.1
            c
        }
    }

    #[inline]
    fn negate(q: Quaternion<f32>) -> Quaternion<f32> {
        unsafe {
            let zero = _mm_setzero_ps();  // SSE
            let q = _mm_setr_ps(q.0, q.1[0], q.1[1], q.1[2]);
            
            let mut c: [f32; 4] = mem::MaybeUninit::uninit().assume_init();
            _mm_storeu_ps(&mut c[0], _mm_sub_ps(zero, q));  // SSE
            mem::transmute::<[f32; 4], Quaternion<f32>>(c)
        }
    }
}

#[cfg(all(feature = "std", feature = "simd"))]
impl<T> FloatSimd<f64> for T {
    #[inline]
    fn sum(q: Quaternion<f64>) -> f64 {
        unsafe {
            let q_h = _mm_set_pd(q.0, q.1[0]); // SSE2
            let q_l = _mm_loadu_pd(&q.1[1]);   // SSE2

            let mut c: [f64; 2] = mem::MaybeUninit::uninit().assume_init();
            _mm_storeu_pd(&mut c[0], _mm_hadd_pd(q_h, q_l));  // SSE2, SSE3
            c[0] + c[1]
        }
    }

    #[inline]
    fn add(a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
        unsafe {
            let a = _mm256_setr_pd(a.0, a.1[0], a.1[1], a.1[2]);  // AVX
            let b = _mm256_setr_pd(b.0, b.1[0], b.1[1], b.1[2]);
            
            let mut c: [f64; 4] = mem::MaybeUninit::uninit().assume_init();
            _mm256_storeu_pd(&mut c[0], _mm256_add_pd(a, b));  // AVX
            mem::transmute::<[f64; 4], Quaternion<f64>>(c)
        }
    }

    #[inline]
    fn sub(a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
        unsafe {
            let a = _mm256_setr_pd(a.0, a.1[0], a.1[1], a.1[2]);  // AVX
            let b = _mm256_setr_pd(b.0, b.1[0], b.1[1], b.1[2]);
            
            let mut c: [f64; 4] = mem::MaybeUninit::uninit().assume_init();
            _mm256_storeu_pd(&mut c[0], _mm256_sub_pd(a, b));  // AVX
            mem::transmute::<[f64; 4], Quaternion<f64>>(c)
        }
    }

    #[inline]
    fn scale(s: f64, q: Quaternion<f64>) -> Quaternion<f64> {
        unsafe {
            let s = _mm256_set1_pd(s);  // AVX
            let q = _mm256_setr_pd(q.0, q.1[0], q.1[1], q.1[2]);  // AVX
            
            let mut c: [f64; 4] = mem::MaybeUninit::uninit().assume_init();
            _mm256_storeu_pd(&mut c[0], _mm256_mul_pd(s, q));  // AVX
            mem::transmute::<[f64; 4], Quaternion<f64>>(c)
        }
    }

    #[inline]
    fn scale_add(s: f64, a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
        unsafe {
            let s = _mm256_set1_pd(s);  // AVX
            let a = _mm256_setr_pd(a.0, a.1[0], a.1[1], a.1[2]);  // AVX
            let b = _mm256_setr_pd(b.0, b.1[0], b.1[1], b.1[2]);
            
            let mut c: [f64; 4] = mem::MaybeUninit::uninit().assume_init();
            if cfg!(feature = "fma") {
                _mm256_storeu_pd(&mut c[0], _mm256_fmadd_pd(s, a, b));  // AVX, FMA
            } else {
                _mm256_storeu_pd(&mut c[0], _mm256_add_pd(_mm256_mul_pd(s, a), b));  // AVX
            }
            mem::transmute::<[f64; 4], Quaternion<f64>>(c)
        }
    }

    #[inline]
    fn hadamard(a: Quaternion<f64>, b: Quaternion<f64>) -> Quaternion<f64> {
        unsafe {
            let a = _mm256_setr_pd(a.0, a.1[0], a.1[1], a.1[2]);  // AVX
            let b = _mm256_setr_pd(b.0, b.1[0], b.1[1], b.1[2]);
            
            let mut c: [f64; 4] = mem::MaybeUninit::uninit().assume_init();
            _mm256_storeu_pd(&mut c[0], _mm256_mul_pd(a, b));  // AVX
            mem::transmute::<[f64; 4], Quaternion<f64>>(c)
        }
    }

    #[inline]
    fn hadamard_add(a: Quaternion<f64>, b: Quaternion<f64>, c: Quaternion<f64>) -> Quaternion<f64> {
        unsafe {
            let a = _mm256_setr_pd(a.0, a.1[0], a.1[1], a.1[2]);  // AVX
            let b = _mm256_setr_pd(b.0, b.1[0], b.1[1], b.1[2]);
            let c = _mm256_setr_pd(c.0, c.1[0], c.1[1], c.1[2]);
            
            let mut d: [f64; 4] = mem::MaybeUninit::uninit().assume_init();
            if cfg!(feature = "fma") {
                _mm256_storeu_pd(&mut d[0], _mm256_fmadd_pd(a, b, c));  // AVX, FMA
            } else {
                _mm256_storeu_pd(&mut d[0], _mm256_add_pd(_mm256_mul_pd(a, b), c));  // AVX
            }
            mem::transmute::<[f64; 4], Quaternion<f64>>(d)
        }
    }

    #[inline]
    fn dot(a: Quaternion<f64>, b: Quaternion<f64>) -> f64 {
        unsafe {
            let a_h = _mm_setr_pd(a.0, a.1[0]); // SSE2
            let a_l = _mm_loadu_pd(&a.1[1]);    // SSE2
            let b_h = _mm_setr_pd(b.0, b.1[0]);
            let b_l = _mm_loadu_pd(&b.1[1]);

            let mut c: [f64; 2] = mem::MaybeUninit::uninit().assume_init();
            let dp_h = _mm_dp_pd(a_h, b_h, 0b00110010);  // SSE4.1
            let dp_l = _mm_dp_pd(a_l, b_l, 0b00110001);
            _mm_storeu_pd(&mut c[0], _mm_move_sd(dp_h, dp_l));  // SSE2
            c[0] + c[1]
        }
    }

    #[inline]
    fn negate(q: Quaternion<f64>) -> Quaternion<f64> {
        unsafe {
            let zero = _mm256_setzero_pd();  // AVX
            let q = _mm256_setr_pd(q.0, q.1[0], q.1[1], q.1[2]);  // AVX
            
            let mut c: [f64; 4] = mem::MaybeUninit::uninit().assume_init();
            _mm256_storeu_pd(&mut c[0], _mm256_sub_pd(zero, q));  // AVX
            mem::transmute::<[f64; 4], Quaternion<f64>>(c)
        }
    }
}

#[cfg(not(all(feature = "std", feature = "simd")))]
impl<T: Float> FloatSimd<T> for T {
    #[inline]
    fn sum(q: Quaternion<T>) -> T {
        q.0 + super::sum_vec(q.1)
    }

    #[inline]
    fn add(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T> {
        ( a.0 + b.0, super::add_vec(a.1, b.1) )
    }

    #[inline]
    fn sub(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T> {
        ( a.0 - b.0, super::sub_vec(a.1, b.1) )
    }

    #[inline]
    fn scale(s: T, q: Quaternion<T>) -> Quaternion<T> {
        ( s * q.0, super::scale_vec(s, q.1) )
    }

    #[inline]
    fn scale_add(s: T, a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T> {
        ( super::mul_add(s, a.0, b.0), super::scale_add_vec(s, a.1, b.1) )
    }

    #[inline]
    fn hadamard(a: Quaternion<T>, b: Quaternion<T>) -> Quaternion<T> {
        ( a.0 * b.0, super::hadamard_vec(a.1, b.1) )
    }

    #[inline]
    fn hadamard_add(a: Quaternion<T>, b: Quaternion<T>, c: Quaternion<T>) -> Quaternion<T> {
        ( super::mul_add(a.0, b.0, c.0), super::hadamard_add_vec(a.1, b.1, c.1) )
    }

    #[inline]
    fn dot(a: Quaternion<T>, b: Quaternion<T>) -> T {
        super::sum( super::hadamard(a, b) )
    }

    #[inline]
    fn negate(q: Quaternion<T>) -> Quaternion<T> {
        ( -q.0, super::negate_vec(q.1) )
    }
}