//! 四元数と純虚四元数で共通する処理をまとめる

// Pure Quaternionのexpは使うけど，lnの引数にPure Quaternionを取ることは
// 普通無いと思うから実装しない．同じ理由でpowの実装も無し．
// どうしてもPure Quaternionを入れたければln((0.0, pure))みたいにすれば良い．

use super::{Vector3, Quaternion, dot, norm, cross, conj, mul_add, ZERO_VECTOR};
use super::Float;


/// This trait provides operations common to Quaternions and Pure Quaternions (Vector3).
pub trait QuaternionOps<T>: Copy {
    fn sum(self) -> T;
    /// `self + rhs`
    fn add(self, rhs: Self) -> Self;
    /// `self - rhs`
    fn sub(self, rhs: Self) -> Self;
    /// `s * self`
    fn scale(self, s: T) -> Self;
    /// `s * self + b`
    fn scale_add(self, s: T, b: Self) -> Self;
    /// `self ∘ rhs`
    fn hadamard(self, rhs: Self) -> Self;
    /// `self ∘ b + c`
    fn hadamard_add(self, b: Self, c: Self) -> Self;
    fn normalize(self) -> Self;
    /// `-self`
    fn negate(self) -> Self;
    /// `self * rhs`
    fn mul(self, rhs: Self) -> Quaternion<T>;
    fn inv(self) -> Self;
    fn exp(self) -> Quaternion<T>;
}

impl<T: Float> QuaternionOps<T> for Vector3<T> {
    #[inline]
    fn sum(self) -> T {
        self[0] + self[1] + self[2]
    }

    #[inline]
    fn add(self, rhs: Self) -> Self {
        [
            self[0] + rhs[0], 
            self[1] + rhs[1], 
            self[2] + rhs[2] 
        ]
    }

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        [
            self[0] - rhs[0], 
            self[1] - rhs[1], 
            self[2] - rhs[2]
        ]
    }

    #[inline]
    fn scale(self, s: T) -> Self {
        [
            s * self[0],
            s * self[1],
            s * self[2]
        ]
    }

    #[inline]
    fn scale_add(self, s: T, b: Self) -> Self {
        [
            mul_add(s, self[0], b[0]),
            mul_add(s, self[1], b[1]),
            mul_add(s, self[2], b[2]),
        ]
    }

    #[inline]
    fn hadamard(self, rhs: Self) -> Self {
        [
            self[0] * rhs[0],
            self[1] * rhs[1],
            self[2] * rhs[2]
        ]
    }

    #[inline]
    fn hadamard_add(self, b: Self, c: Self) -> Self {
        [
            mul_add(self[0], b[0], c[0]),
            mul_add(self[1], b[1], c[1]),
            mul_add(self[2], b[2], c[2]),
        ]
    }

    // 零ベクトルを入力した時は零ベクトルを返す
    #[inline]
    fn normalize(self) -> Self {
        let coef = norm(self).recip();
        if coef.is_infinite() {
            ZERO_VECTOR()
        } else {
            self.scale(coef)
        }
    }

    #[inline]
    fn negate(self) -> Self {
        [ -self[0], -self[1], -self[2] ]
    }

    // Product of Pure Quaternions
    #[inline]
    fn mul(self, rhs: Self) -> Quaternion<T> {
        ( -dot(self, rhs), cross(self, rhs) )
    }

    // 零ベクトルで特異点
    #[inline]
    fn inv(self) -> Self {
        self.negate().scale( dot(self, self).recip() )
    }

    #[inline]
    fn exp(self) -> Quaternion<T> {
        let norm_v = norm(self);
        let (sin, cos) = norm_v.sin_cos();
        ( cos, self.scale(sin / norm_v) )
    }
}

impl<T: Float> QuaternionOps<T> for Quaternion<T> {
    #[inline]
    fn sum(self) -> T {
        self.0 + self.1.sum()
    }

    #[inline]
    fn add(self, rhs: Self) -> Self {
        ( self.0 + rhs.0, self.1.add(rhs.1) )
    }

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        ( self.0 - rhs.0,  self.1.sub(rhs.1) )
    }

    #[inline]
    fn scale(self, s: T) -> Self {
        ( s * self.0, self.1.scale(s) )
    }

    #[inline]
    fn scale_add(self, s: T, b: Self) -> Self {
        ( mul_add(s, self.0, b.0), self.1.scale_add(s, b.1) )
    }

    #[inline]
    fn hadamard(self, rhs: Self) -> Self {
        ( self.0 * rhs.0, self.1.hadamard(rhs.1) )
    }

    #[inline]
    fn hadamard_add(self, b: Self, c: Self) -> Self {
        ( mul_add(self.0, b.0, c.0), self.1.hadamard_add(b.1, c.1) )
    }

    #[inline]
    fn normalize(self) -> Self {
        self.scale( norm(self).recip() )
    }

    #[inline]
    fn negate(self) -> Self {
        ( -self.0, self.1.negate() )
    }

    #[inline]
    fn mul(self, rhs: Self) -> Quaternion<T> {
        let self0_b = rhs.scale(self.0);
        (
            self0_b.0 - dot(self.1, rhs.1),
            self.1.scale_add(rhs.0, self0_b.1).add( cross(self.1, rhs.1) )
        )
    }

    #[inline]
    fn inv(self) -> Self {
        conj(self).scale( dot(self, self).recip() )
    }

    #[inline]
    fn exp(self) -> Quaternion<T> {
        let norm_v = norm(self.1);
        let (sin, cos) = norm_v.sin_cos();
        let coef = self.0.exp();
        ( coef * cos, self.1.scale((coef * sin) / norm_v) )
    }
}
