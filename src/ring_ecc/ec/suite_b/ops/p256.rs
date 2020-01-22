// Copyright 2016 Brian Smith.
//
// Permission to use, copy, modify, and/or distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHORS DISCLAIM ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY
// SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
// OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
// CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

use super::{
    elem_sqr_mul, elem_sqr_mul_acc, Modulus, *,
};
use core::marker::PhantomData;

macro_rules! p256_limbs {
    [ $($limb:expr),+ ] => {
        limbs![$($limb),+, 0, 0, 0, 0]
    };
}

pub static COMMON_OPS: CommonOps = CommonOps {
    num_limbs: 256 / LIMB_BITS,

    q: Modulus {
        p: p256_limbs![
            0xffffffff, 0xffffffff, 0xffffffff, 0x00000000, 0x00000000, 0x00000000, 0x00000001,
            0xffffffff
        ],
        rr: p256_limbs![
            0x00000003, 0x00000000, 0xffffffff, 0xfffffffb, 0xfffffffe, 0xffffffff, 0xfffffffd,
            0x00000004
        ],
    },

    n: Elem {
        limbs: p256_limbs![
            0xfc632551, 0xf3b9cac2, 0xa7179e84, 0xbce6faad, 0xffffffff, 0xffffffff, 0x00000000,
            0xffffffff
        ],
        m: PhantomData,
        encoding: PhantomData, // Unencoded
    },

    a: Elem {
        limbs: p256_limbs![
            0xfffffffc, 0xffffffff, 0xffffffff, 0x00000003, 0x00000000, 0x00000000, 0x00000004,
            0xfffffffc
        ],
        m: PhantomData,
        encoding: PhantomData, // R
    },
    b: Elem {
        limbs: p256_limbs![
            0x29c4bddf, 0xd89cdf62, 0x78843090, 0xacf005cd, 0xf7212ed6, 0xe5a220ab, 0x04874834,
            0xdc30061d
        ],
        m: PhantomData,
        encoding: PhantomData, // R
    },

    elem_add_impl: GFp_nistz256_add,
    elem_mul_mont: GFp_nistz256_mul_mont,
    elem_sqr_mont: GFp_nistz256_sqr_mont,

    point_add_jacobian_impl: GFp_nistz256_point_add,
};

pub static PRIVATE_KEY_OPS: PrivateKeyOps = PrivateKeyOps {
    common: &COMMON_OPS,
    elem_inv_squared: p256_elem_inv_squared,
    point_mul_impl: GFp_nistz256_point_mul,
};

fn p256_elem_inv_squared(a: &Elem<R>) -> Elem<R> {
    // Calculate a**-2 (mod q) == a**(q - 3) (mod q)
    //
    // The exponent (q - 3) is:
    //
    //    0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc

    #[inline]
    fn sqr_mul(a: &Elem<R>, squarings: usize, b: &Elem<R>) -> Elem<R> {
        elem_sqr_mul(&COMMON_OPS, a, squarings, b)
    }

    #[inline]
    fn sqr_mul_acc(a: &mut Elem<R>, squarings: usize, b: &Elem<R>) {
        elem_sqr_mul_acc(&COMMON_OPS, a, squarings, b)
    }

    let b_1 = &a;
    let b_11 = sqr_mul(b_1, 1, b_1);
    let b_111 = sqr_mul(&b_11, 1, b_1);
    let f_11 = sqr_mul(&b_111, 3, &b_111);
    let fff = sqr_mul(&f_11, 6, &f_11);
    let fff_111 = sqr_mul(&fff, 3, &b_111);
    let fffffff_11 = sqr_mul(&fff_111, 15, &fff_111);
    let ffffffff = sqr_mul(&fffffff_11, 2, &b_11);

    // ffffffff00000001
    let mut acc = sqr_mul(&ffffffff, 31 + 1, b_1);

    // ffffffff00000001000000000000000000000000ffffffff
    sqr_mul_acc(&mut acc, 96 + 32, &ffffffff);

    // ffffffff00000001000000000000000000000000ffffffffffffffff
    sqr_mul_acc(&mut acc, 32, &ffffffff);

    // ffffffff00000001000000000000000000000000fffffffffffffffffffffff_11
    sqr_mul_acc(&mut acc, 30, &fffffff_11);

    // ffffffff00000001000000000000000000000000fffffffffffffffffffffffc
    COMMON_OPS.elem_square(&mut acc);
    COMMON_OPS.elem_square(&mut acc);

    acc
}

pub static PUBLIC_KEY_OPS: PublicKeyOps = PublicKeyOps {
    common: &COMMON_OPS,
};

pub static PUBLIC_SCALAR_OPS: PublicScalarOps = PublicScalarOps {
    public_key_ops: &PUBLIC_KEY_OPS,
    private_key_ops: &PRIVATE_KEY_OPS,

    q_minus_n: Elem {
        limbs: p256_limbs![0x039cdaae, 0x0c46353d, 0x58e8617b, 0x43190553, 0, 0, 0, 0],
        m: PhantomData,
        encoding: PhantomData, // Unencoded
    },
};

extern "C" {
    fn GFp_nistz256_add(
        r: *mut Limb,   // [COMMON_OPS.num_limbs]
        a: *const Limb, // [COMMON_OPS.num_limbs]
        b: *const Limb, // [COMMON_OPS.num_limbs]
    );
    fn GFp_nistz256_mul_mont(
        r: *mut Limb,   // [COMMON_OPS.num_limbs]
        a: *const Limb, // [COMMON_OPS.num_limbs]
        b: *const Limb, // [COMMON_OPS.num_limbs]
    );
    fn GFp_nistz256_sqr_mont(
        r: *mut Limb,   // [COMMON_OPS.num_limbs]
        a: *const Limb, // [COMMON_OPS.num_limbs]
    );

    fn GFp_nistz256_point_add(
        r: *mut Limb,   // [3][COMMON_OPS.num_limbs]
        a: *const Limb, // [3][COMMON_OPS.num_limbs]
        b: *const Limb, // [3][COMMON_OPS.num_limbs]
    );
    fn GFp_nistz256_point_mul(
        r: *mut Limb,          // [3][COMMON_OPS.num_limbs]
        p_scalar: *const Limb, // [COMMON_OPS.num_limbs]
        p_x: *const Limb,      // [COMMON_OPS.num_limbs]
        p_y: *const Limb,      // [COMMON_OPS.num_limbs]
    );
}

#[cfg(feature = "internal_benches")]
mod internal_benches {
    use super::{super::internal_benches::*, *};

    bench_curve!(&[
        Scalar {
            limbs: LIMBS_1,
            m: PhantomData,
            encoding: PhantomData,
        },
        Scalar {
            limbs: LIMBS_ALTERNATING_10,
            m: PhantomData,
            encoding: PhantomData,
        },
        Scalar {
            // n - 1
            limbs: p256_limbs![
                0xfc632551 - 1,
                0xf3b9cac2,
                0xa7179e84,
                0xbce6faad,
                0xffffffff,
                0xffffffff,
                0x00000000,
                0xffffffff
            ],
            m: PhantomData,
            encoding: PhantomData,
        },
    ]);
}
