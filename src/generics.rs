//! 四元数と純虚四元数で共通する処理をまとめる

// Pure Quaternionのexpは使うけど，lnの引数にPure Quaternionを取ることは
// 普通無いと思うから実装しない．同じ理由でpowの実装も無し．
// どうしてもPure Quaternionを入れたければln((0.0, pure))みたいにすれば良い．

use super::{Vector3, Quaternion, Float};


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
    fn norm(self) -> T;
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
            super::mul_add(s, self[0], b[0]),
            super::mul_add(s, self[1], b[1]),
            super::mul_add(s, self[2], b[2]),
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
            super::mul_add(self[0], b[0], c[0]),
            super::mul_add(self[1], b[1], c[1]),
            super::mul_add(self[2], b[2], c[2]),
        ]
    }

    #[inline]
    fn norm(self) -> T {
        #[cfg(feature = "norm-sqrt")]
        {
            super::dot(self, self).sqrt()
        }

        #[cfg(not(feature = "norm-sqrt"))]
        {
            let mut s = self[0];
            for val in self.iter().skip(1) {
                s = super::pythag(s, *val);
            }
            s
        }
    }

    #[inline]
    fn negate(self) -> Self {
        [ -self[0], -self[1], -self[2] ]
    }

    // Product of Pure Quaternions
    #[inline]
    fn mul(self, rhs: Self) -> Quaternion<T> {
        ( -super::dot(self, rhs), super::cross(self, rhs) )
    }

    // 零ベクトルで特異点
    #[inline]
    fn inv(self) -> Self {
        self.negate().scale( super::dot(self, self).recip() )
    }

    #[inline]
    fn exp(self) -> Quaternion<T> {
        let norm_v = self.norm();
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
        ( super::mul_add(s, self.0, b.0), self.1.scale_add(s, b.1) )
    }

    #[inline]
    fn hadamard(self, rhs: Self) -> Self {
        ( self.0 * rhs.0, self.1.hadamard(rhs.1) )
    }

    #[inline]
    fn hadamard_add(self, b: Self, c: Self) -> Self {
        ( super::mul_add(self.0, b.0, c.0), self.1.hadamard_add(b.1, c.1) )
    }

    #[inline]
    fn norm(self) -> T {
        #[cfg(feature = "norm-sqrt")]
        {
            super::dot(self, self).sqrt()
        }

        #[cfg(not(feature = "norm-sqrt"))]
        {
            let mut s = self.0;
            for val in self.1 {
                s = super::pythag(s, val)
            }
            s
        }
    }

    #[inline]
    fn negate(self) -> Self {
        ( -self.0, self.1.negate() )
    }

    #[inline]
    fn mul(self, rhs: Self) -> Quaternion<T> {
        let self0_b = rhs.scale(self.0);
        (
            self0_b.0 - super::dot(self.1, rhs.1),
            self.1.scale_add(rhs.0, self0_b.1).add( super::cross(self.1, rhs.1) )
        )
    }

    #[inline]
    fn inv(self) -> Self {
        super::conj(self).scale( super::dot(self, self).recip() )
    }

    #[inline]
    fn exp(self) -> Quaternion<T> {
        let norm_v = self.1.norm();
        let (sin, cos) = norm_v.sin_cos();
        let coef = self.0.exp();
        ( coef * cos, self.1.scale((coef * sin) / norm_v) )
    }
}
