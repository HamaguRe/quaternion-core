//! Private functions

use core::{mem, mem::MaybeUninit};
use super::{Float, Vector3, Quaternion};


#[inline]
#[allow(non_snake_case)]
pub fn IDENTITY<T: Float>() -> Quaternion<T> {
    (T::one(), [T::zero(); 3])
}

#[inline]
#[allow(non_snake_case)]
pub fn ZERO_VECTOR<T: Float>() -> Vector3<T> {
    [T::zero(); 3]
}

/// 定数呼び出し以外に使わないのでエラー処理を省略．
#[inline(always)]
pub fn cast<T: Float>(x: f64) -> T {
    num_traits::cast::<f64, T>(x).unwrap()
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

/// （主に呼び出し側の特異点近傍で）NaNにならないように定義域外の値をカットする．
/// 
/// 勝手に値をカットしても問題ないところでのみ使うこと．
#[inline]
pub fn acos_safe<T: Float>(x: T) -> T {
    // FloatConstを使いたくないからこの実装とする．
    (if x.abs() > T::one() { x.signum() } else { x }).acos()
}

/// sinc(x) := sin(x) / x
/// 
/// sinc関数ではx=0において特異点とならずに計算できる
/// （sin(x)/xにおけるx=0は可除去特異点）．
/// 
/// sinc(Inf or NaN) is NaN
#[inline]
pub fn sinc<T: Float>(x: T) -> T {
    // チェビシェフ近似多項式の係数
    // [-1, 1]の範囲で最大誤差1.55e-15程度
    const COEF: [f64; 7] = [
        1.5639670891687275e-10,  // s12
        -2.504349140508566e-8,   // s10
        2.7557233897823613e-6,   // s8
        -0.00019841269414627052, // s6
        0.00833333333232666,     // s4
        -0.16666666666657756,    // s2
        0.9999999999999986,      // s0
    ];

    // 特異点回避できれば十分
    if x.abs() <= T::one() { // <-- false here when x is NaN or Inf.
        // y = s0 + s2*x^2 + s4*x^4 + ... + s12*x^12
        // 上式において，z=x^2とおけば以下のように変形できる．
        // y = (((((s12*z + s10)*z + s8)*z + s6)*z + s4)*z + s2)*z + s0
        // この変形はホーナー法（Horner's rule）と呼ばれるものである．
        let z = x * x;
        let mut y = cast(COEF[0]);
        for i in 1..COEF.len() {
            y = mul_add(y, z, cast(COEF[i]));
        }
        y
    } else {
        x.sin() / x
    }
}

/// 配列内の最大値とそのインデックスを返す．
#[inline]
pub fn max4<T: Float>(nums: [T; 4]) -> (usize, T) {
    let mut index = 0;
    let mut max_num = nums[0];
    for (i, num) in nums.iter().enumerate().skip(1) {
        if *num > max_num {
            max_num = *num;
            index = i;
        }
    }
    (index, max_num)
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

    // aの絶対値が最大となるインデックスを探す（working_arrayにはaの絶対値を入れる）
    working_array[0] = a[0].abs();
    let mut i_max: usize = 0;
    let mut max_val = working_array[0];
    for (i, val) in a.iter().enumerate().skip(1) {
        working_array[i] = val.abs();
        if working_array[i] > max_val {
            max_val = working_array[i];
            i_max = i;
        }
    }
    // working_arrayの中央値を探す
    let idx1 = (i_max + 1) % 3;
    let idx2 = (i_max + 2) % 3;
    let i_med = if working_array[idx1] > working_array[idx2] {
        idx1
    } else {
        idx2
    };

    let norm_inv = {
        #[cfg(feature = "norm-sqrt")]
        {
            (a[i_med] * a[i_med] + a[i_max] * a[i_max]).sqrt().recip()
        }

        #[cfg(not(feature = "norm-sqrt"))]
        {
            pythag(a[i_med], a[i_max]).recip()
        }
    };
    working_array[i_med] = -a[i_max] * norm_inv;
    working_array[i_max] =  a[i_med] * norm_inv;
    working_array[3 - (i_med + i_max)] = T::zero();

    working_array
}

// f32,f64のメソッドには同様の機能を提供するhypot()が存在するが，
// このクレートはマイコン上での動作も想定しているため内部で平方根の
// 計算を行わないMoler-Morrison algorithmを採用した．収束が速いから
// ノルムの計算に使用しても丸め誤差が蓄積しにくいというのも理由の一つ．
// あと，単純に実装が美しい．
// 
/// Moler-Morrison algorithmにより，ピタゴラスの定理
/// `c^2 = a^2 + b^2`
/// を満たすcを求める．
/// 
/// `c = (a*a + b*b).sqrt()`のように実装した場合，有害なアンダーフロー
/// やオーバーフローが発生する可能性がある．
/// Moler-Morrison algorithmは他の方法に比べて速度面で若干劣るものの，
/// 堅牢かつ移植性の高い実装で正確な計算結果を得ることができる．
/// cがオーバーフローしない範囲の値である限り，本関数による演算で
/// オーバーフローが起こることはない．
/// 
/// 反復回数は有効数字6桁なら2回, 20桁なら3回, 60桁なら4回．
/// aとbの絶対値が大きく異なる場合（例えば|a| >> |b|）には，
/// より少ない反復回数で結果が求まる．
#[inline]
pub fn pythag<T: Float>(a: T, b: T) -> T {
    let mut a = a.abs();
    let mut b = b.abs();
    if a < b {
        mem::swap(&mut a, &mut b);
    }
    if b == T::zero() {
        return a;
    }

    let two : T = cast(2.0);
    let four: T = cast(4.0);
    for _ in 0..4 {  // loopにするとinfやNaNが入った際に抜けられなくなる
        let mut s = b / a;
        s = s * s;
        let tmp = four + s;
        if tmp == four {  // 収束判定（これでちゃんとbreakできる）
            break;
        }
        s = s / tmp;
        a = mul_add(two, a * s, a);
        b = s * b;
    }
    a
}
