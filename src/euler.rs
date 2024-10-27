//! Implementation for Euler angles

// Extrinsic rotationとIntrinsic rotationは逆順の関係にあることを利用して実装量を減らした．
// Extrinsic rotationのほうが導出がやりやすいのでExtrinsic rotationベースで実装してある．

use super::{Float, FloatConst};
use super::{RotationSequence, Vector3, Quaternion};
use super::{scale, pfs::cast, to_dcm};

/// 特異点判定の閾値
/// 
/// THRESHOLD = sin(89.9 deg) = cos(0.1 deg)
const THRESHOLD: f64 = 0.9999984;

/// Extrinsic rotationのオイラー角を四元数に変換する．
#[inline]
pub fn from_extrinsic_euler_angles<T>(rs: RotationSequence, angles: Vector3<T>) -> Quaternion<T>
where T: Float {
    use RotationSequence::*;
    let [alpha, beta, gamma] = scale(cast(0.5), angles);
    let (sinb, cosb) = beta.sin_cos();

    let (q0, [q1, q2, q3]): Quaternion<T>;
    match rs {
        // ------ Proper Euler angles ------ //
        ZXZ | XYX | YZY | ZYZ | XZX | YXY => {
            let (sin_apg, cos_apg) = (alpha + gamma).sin_cos();
            let (sin_amg, cos_amg) = (alpha - gamma).sin_cos();
            match rs {
                ZXZ => {
                    q0 = cos_apg * cosb;
                    q1 = cos_amg * sinb;
                    q2 = sin_amg * -sinb;
                    q3 = sin_apg * cosb;
                },
                XYX => {
                    q0 = cos_apg * cosb;
                    q1 = sin_apg * cosb;
                    q2 = cos_amg * sinb;
                    q3 = sin_amg * -sinb;
                },
                YZY => {
                    q0 = cos_apg * cosb;
                    q1 = sin_amg * -sinb;
                    q2 = sin_apg * cosb;
                    q3 = cos_amg * sinb;
                },
                ZYZ => {
                    q0 = cos_apg * cosb;
                    q1 = sin_amg * sinb;
                    q2 = cos_amg * sinb;
                    q3 = sin_apg * cosb;
                },
                XZX => {
                    q0 = cos_apg * cosb;
                    q1 = sin_apg * cosb;
                    q2 = sin_amg * sinb;
                    q3 = cos_amg * sinb;
                },
                YXY => {
                    q0 = cos_apg * cosb;
                    q1 = cos_amg * sinb;
                    q2 = sin_apg * cosb;
                    q3 = sin_amg * sinb;
                },
                _ => unreachable!(),
            }
        },
        // ------- Tait–Bryan angles ------- //
        XYZ | YZX | ZXY | XZY | ZYX | YXZ => {
            let (sina, cosa) = alpha.sin_cos();
            let (sing, cosg) = gamma.sin_cos();
            let sina_sinb = sina * sinb;
            let cosa_cosb = cosa * cosb;
            let sina_cosb = sina * cosb;
            let cosa_sinb = cosa * sinb;
            match rs {
                XYZ => {
                    q0 = cosa_cosb*cosg + sina_sinb*sing;
                    q1 = sina_cosb*cosg - cosa_sinb*sing;
                    q2 = cosa_sinb*cosg + sina_cosb*sing;
                    q3 = cosa_cosb*sing - sina_sinb*cosg;
                },
                YZX => {
                    q0 = cosa_cosb*cosg + sina_sinb*sing;
                    q1 = cosa_cosb*sing - sina_sinb*cosg;
                    q2 = sina_cosb*cosg - cosa_sinb*sing;
                    q3 = cosa_sinb*cosg + sina_cosb*sing;
                },
                ZXY => {
                    q0 = cosa_cosb*cosg + sina_sinb*sing;
                    q1 = cosa_sinb*cosg + sina_cosb*sing;
                    q2 = cosa_cosb*sing - sina_sinb*cosg;
                    q3 = sina_cosb*cosg - cosa_sinb*sing;
                },
                XZY => {
                    q0 = cosa_cosb*cosg - sina_sinb*sing;
                    q1 = sina_cosb*cosg + cosa_sinb*sing;
                    q2 = cosa_cosb*sing + sina_sinb*cosg;
                    q3 = cosa_sinb*cosg - sina_cosb*sing;
                },
                ZYX => {
                    q0 = cosa_cosb*cosg - sina_sinb*sing;
                    q1 = cosa_cosb*sing + sina_sinb*cosg;
                    q2 = cosa_sinb*cosg - sina_cosb*sing;
                    q3 = sina_cosb*cosg + cosa_sinb*sing;
                },
                YXZ => {
                    q0 = cosa_cosb*cosg - sina_sinb*sing;
                    q1 = cosa_sinb*cosg - sina_cosb*sing;
                    q2 = sina_cosb*cosg + cosa_sinb*sing;
                    q3 = cosa_cosb*sing + sina_sinb*cosg;
                },
                _ => unreachable!(),
            }
        },
    }

    (q0, [q1, q2, q3])
}

/// Intrinsic rotationのオイラー角を四元数に変換する．
#[inline]
pub fn from_intrinsic_euler_angles<T>(rs: RotationSequence, angles: Vector3<T>) -> Quaternion<T>
where T: Float {
    use RotationSequence::*;

    let angles = [angles[2], angles[1], angles[0]];
    match rs {
        // ------ Proper Euler angles ------ //
        ZXZ | XYX | YZY | ZYZ | XZX | YXY => {
            from_extrinsic_euler_angles(rs, angles)
        },
        // ------- Tait–Bryan angles ------- //
        XYZ => from_extrinsic_euler_angles(ZYX, angles),
        YZX => from_extrinsic_euler_angles(XZY, angles),
        ZXY => from_extrinsic_euler_angles(YXZ, angles),
        XZY => from_extrinsic_euler_angles(YZX, angles),
        ZYX => from_extrinsic_euler_angles(XYZ, angles),
        YXZ => from_extrinsic_euler_angles(ZXY, angles),
    }
}

/// 四元数からExtrinsicオイラー角に変換する．
/// 
/// Properの場合は第二回転角の正弦が0（0, ±PI）のとき特異点となり，
/// Tait-Bryanの場合は第二回転角の余弦が0（±PI/2）のとき特異点となる．
/// 
/// 特異点では第一回転角を0[rad]にする．
#[inline]
pub fn to_extrinsic_euler_angles<T>(rs: RotationSequence, q: Quaternion<T>) -> Vector3<T>
where T: Float + FloatConst {
    use RotationSequence::*;
    let [
        [m11, m12, m13],
        [m21, m22, m23],
        [m31, m32, m33],
    ] = to_dcm(q);

    let [e1, e2, e3]: Vector3<T>;
    match rs {
        // ------ Proper Euler angles ------ //
        ZXZ => {
            if m33.abs() < cast(THRESHOLD) {
                e1 = m31.atan2(m32);
                e2 = m33.acos();
                e3 = m13.atan2(-m23);
            } else {  // Singularity
                e1 = T::zero();
                e2 = if m33.is_sign_positive() {T::zero()} else {T::PI()};
                e3 = m21.atan2(m11);
            }
        },
        XYX => {
            if m11.abs() < cast(THRESHOLD) {
                e1 = m12.atan2(m13);
                e2 = m11.acos();
                e3 = m21.atan2(-m31);
            } else {  // Singularity
                e1 = T::zero();
                e2 = if m11.is_sign_positive() {T::zero()} else {T::PI()};
                e3 = m32.atan2(m22);
            }
        },
        YZY => {
            if m22.abs() < cast(THRESHOLD) {
                e1 = m23.atan2(m21);
                e2 = m22.acos();
                e3 = m32.atan2(-m12);
            } else {  // Singularity
                e1 = T::zero();
                e2 = if m22.is_sign_positive() {T::zero()} else {T::PI()};
                e3 = m13.atan2(m33);
            }
        },
        ZYZ => {
            if m33.abs() < cast(THRESHOLD) {
                e1 = m32.atan2(-m31);
                e2 = m33.acos();
                e3 = m23.atan2(m13);
            } else {  // Singularity
                e1 = T::zero();
                e2 = if m33.is_sign_positive() {T::zero()} else {T::PI()};
                e3 = (-m12).atan2(m22);
            }
        },
        XZX => {
            if m11.abs() < cast(THRESHOLD) {
                e1 = m13.atan2(-m12);
                e2 = m11.acos();
                e3 = m31.atan2(m21);
            } else {  // Singularity
                e1 = T::zero();
                e2 = if m11.is_sign_positive() {T::zero()} else {T::PI()};
                e3 = (-m23).atan2(m33);
            }
        },
        YXY => {
            if m22.abs() < cast(THRESHOLD) {
                e1 = m21.atan2(-m23);
                e2 = m22.acos();
                e3 = m12.atan2(m32);
            } else {  // Singularity
                e1 = T::zero();
                e2 = if m22.is_sign_positive() {T::zero()} else {T::PI()};
                e3 = (-m31).atan2(m11);
            }
        },
        // ------- Tait–Bryan angles ------- //
        XYZ => {
            if m31.abs() < cast(THRESHOLD) {
                e1 = m32.atan2(m33);
                e2 = (-m31).asin();
                e3 = m21.atan2(m11);
            } else {  // Singularity
                e1 = T::zero();
                e2 = T::FRAC_PI_2().copysign(-m31);
                e3 = (-m12).atan2(m22);
            }
        },
        YZX => {
            if m12.abs() < cast(THRESHOLD) {
                e1 = m13.atan2(m11);
                e2 = (-m12).asin();
                e3 = m32.atan2(m22);
            } else {  // Singularity
                e1 = T::zero();
                e2 = T::FRAC_PI_2().copysign(-m12);
                e3 = (-m23).atan2(m33);
            }
        },
        ZXY => {
            if m23.abs() < cast(THRESHOLD) {
                e1 = m21.atan2(m22);
                e2 = (-m23).asin();
                e3 = m13.atan2(m33);
            } else {  // Singularity
                e1 = T::zero();
                e2 = T::FRAC_PI_2().copysign(-m23);
                e3 = (-m31).atan2(m11);
            }
        },
        XZY => {
            if m21.abs() < cast(THRESHOLD) {
                e1 = (-m23).atan2(m22);
                e2 = m21.asin();
                e3 = (-m31).atan2(m11);
            } else {  // Singularity
                e1 = T::zero();
                e2 = T::FRAC_PI_2().copysign(m21);
                e3 = m13.atan2(m33);
            }
        },
        ZYX => {
            if m13.abs() < cast(THRESHOLD) {
                e1 = (-m12).atan2(m11);
                e2 = m13.asin();
                e3 = (-m23).atan2(m33);
            } else {  // Singularity
                e1 = T::zero();
                e2 = T::FRAC_PI_2().copysign(m13);
                e3 = m32.atan2(m22);
            }
        },
        YXZ => {
            if m32.abs() < cast(THRESHOLD) {
                e1 = (-m31).atan2(m33);
                e2 = m32.asin();
                e3 = (-m12).atan2(m22);
            } else {  // Singularity
                e1 = T::zero();
                e2 = T::FRAC_PI_2().copysign(m32);
                e3 = m21.atan2(m11);
            }
        },
    }
    
    [e1, e2, e3]
}

/// 四元数をIntrinsicオイラー角に変換する．
/// 
/// Properの場合は第二回転角の正弦が0（0, ±PI）のとき特異点となり，
/// Tait-Bryanの場合は第二回転角の余弦が0（±PI/2）のとき特異点となる．
/// 
/// 特異点では第三回転角が0[rad]になる．
#[inline]
pub fn to_intrinsic_euler_angles<T>(rs: RotationSequence, q: Quaternion<T>) -> Vector3<T>
where T: Float + FloatConst {
    use RotationSequence::*;

    let angles = match rs {
        // ------ Proper Euler angles ------ //
        ZXZ | XYX | YZY | ZYZ | XZX | YXY => {
            to_extrinsic_euler_angles(rs, q)
        },
        // ------- Tait–Bryan angles ------- //
        XYZ => to_extrinsic_euler_angles(ZYX, q),
        YZX => to_extrinsic_euler_angles(XZY, q),
        ZXY => to_extrinsic_euler_angles(YXZ, q),
        XZY => to_extrinsic_euler_angles(YZX, q),
        ZYX => to_extrinsic_euler_angles(XYZ, q),
        YXZ => to_extrinsic_euler_angles(ZXY, q),
    };

    [angles[2], angles[1], angles[0]]
}
