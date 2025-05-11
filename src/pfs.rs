//! Private functions

use core::mem::MaybeUninit;
use super::{Float, Vector3};


/// 定数呼び出し以外に使わないのでエラー処理を省略．
#[inline(always)]
pub fn cast<T: Float>(x: f64) -> T {
    num_traits::cast::<f64, T>(x).unwrap()
}

/// 二次元ベクトル`[a, b]`のノルムを計算する．
/// 
/// `norm-sqrt`フィーチャーが有効の場合は`(a*a + b*b).sqrt()`とし，
/// そうでない場合は`a.hypot(b)`として計算する．
#[inline(always)]
pub fn norm2<T: Float>(a: T, b: T) -> T {
    #[cfg(feature = "norm-sqrt")]
    {
        (a*a + b*b).sqrt()
    }

    #[cfg(not(feature = "norm-sqrt"))]
    {
        a.hypot(b)
    }
}

/// `fma` featureを有効にした場合は`s.mul_add(a, b)`として展開され，
/// 有効にしなかった場合は単純な積和`s*a + b`に展開してコンパイルされる．
#[inline(always)]
pub fn mul_add<T: Float>(s: T, a: T, b: T) -> T {
    #[cfg(feature = "fma")]
    {
        s.mul_add(a, b)
    }

    #[cfg(not(feature = "fma"))]
    {
        (s * a) + b
    }
}

/// sinc(x) := sin(x) / x
/// 
/// sinc関数ではx=0において特異点とならずに計算できる（sin(x)/xにおけるx=0は可除特異点）．
/// 
/// sinc(Inf or NaN) is NaN
#[inline]
pub fn sinc<T: Float>(x: T) -> T {
    // テレスコーピング法で求めた多項式係数（ほぼ最良近似式）
    // f64で計算した場合，閉区間[-1, 1]において最大誤差 ±2.22e-16 程度．
    const COEF: [f64; 7] = [
        1.579349121783874e-10,   // s12
        -2.5048466629256012e-8,  // s10
        2.7557294426482064e-6,   // s8
        -0.00019841269754547076, // s6
        0.008333333333188874,    // s4
        -0.16666666666665766,    // s2
        0.9999999999999999,      // s0
    ];

    // 特異点回避できれば十分
    if x.abs() <= T::one() { // <-- false here when x is NaN or Inf.
        // y = s0 + s2*x^2 + s4*x^4 + ... + s12*x^12
        // 上式において，z=x^2とおけば以下のように変形できる．
        // y = (((((s12*z + s10)*z + s8)*z + s6)*z + s4)*z + s2)*z + s0
        // この変形はホーナー法（Horner's rule）と呼ばれるものである．
        let z = x * x;
        let mut y = cast(COEF[0]);
        for i in 1..(COEF.len() - 1) {
            y = mul_add(y, z, cast(COEF[i]));
        }
        y = mul_add(y, z, cast::<T>(COEF[COEF.len()-1]) - T::one());  // 精度改善のためのノウハウ
        y + T::one()
    } else {
        x.sin() / x
    }
}

/// 配列内の最大値とそのインデックスを返す．
#[inline]
pub fn max4<T: Float>(nums: [T; 4]) -> (usize, T) {
    let mut i_max = 0;
    for i in 1..nums.len() {
        if nums[i_max] < nums[i] {
            i_max = i;
        }
    }

    (i_max, nums[i_max])
}

/// `a`に直交し，ノルムが1であるベクトルを返す．ただし`norm(a) > 0`であること．
/// 
/// 【理論】
/// aとbが直交する場合には両者の内積がゼロになるので，そこから考えて内積が0に
/// なるようにbの要素を決めれば良い．
/// 具体的には，aの要素のうち２つを入れ替えて，残り一つの要素をゼロにする．また，
/// 入れ替える要素のうち片方の符号を反転させる．
/// 例えば，`a = [0.5, 1.0, -1.5]`であれば`b = [0.0, 1.5, 1.0]`とすることで
/// aに直交なベクトルbが求まる．
/// 
/// 注意点として，入れ替える要素は片方が必ずゼロ以外で無ければならない．
/// 例えば`a = [1.0, 0.0, 0.0]`のような場合に2つめと3つめの要素を入れ替えて
/// しまうと直交なベクトルとならない．
/// 
/// この関数では，入れ替える要素が一意に決まるように絶対値が最大のものと中間の
/// ものを入れ替え，最大値の符号を反転させる．例えば，`a = [0.5, -0.8, -1.5]`
/// に対しては`b = [0.0, 1.5, -0.8]`とする．
#[inline]
pub fn orthogonal_vector<T: Float>(a: Vector3<T>) -> Vector3<T> {
    let mut working_array: Vector3<T> = unsafe {MaybeUninit::uninit().assume_init()};

    // aの絶対値を入れておく
    for (i, val) in a.iter().enumerate() {
        working_array[i] = val.abs();
    }

    // aの絶対値が最大・最小となるインデックスを探す
    let mut i_min: usize = 1;
    let mut i_max: usize = 0;
    if working_array[0] < working_array[1] {
        i_min = 0;
        i_max = 1;
    }
    if working_array[i_max] < working_array[2] {
        i_max = 2;
    } else {
        if working_array[i_min] > working_array[2] {
            i_min = 2;
        }
    }

    let i_med = 3 - (i_max + i_min);
    let norm_inv = norm2(a[i_med], a[i_max]).recip();

    working_array[i_min] = T::zero();
    working_array[i_med] = -a[i_max] * norm_inv;
    working_array[i_max] =  a[i_med] * norm_inv;

    working_array
}
