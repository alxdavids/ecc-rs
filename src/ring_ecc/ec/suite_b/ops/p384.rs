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

macro_rules! p384_limbs {
    [$($limb:expr),+] => {
        limbs![$($limb),+]
    };
}

pub static P384_GENERATOR: (Elem<R>, Elem<R>) = (
    Elem {
        limbs: p384_limbs![
            0x49c0b528, 0x3dd07566, 0xa0d6ce38, 0x20e378e2, 0x541b4d6e, 0x879c3afc, 0x59a30eff,
            0x64548684, 0x614ede2b, 0x812ff723, 0x299e1513, 0x4d3aadc2
        ],
        m: PhantomData,
        encoding: PhantomData,
    },
    Elem {
        limbs: p384_limbs![
            0x4b03a4fe, 0x23043dad, 0x7bb4a9ac, 0xa1bfa8bf, 0x2e83b050, 0x8bade756, 0x68f4ffd9,
            0xc6c35219, 0x3969a840, 0xdd800226, 0x5a15c5e9, 0x2b78abc2
        ],
        m: PhantomData,
        encoding: PhantomData,
    },
);

pub static COMMON_OPS: CommonOps = CommonOps {
    num_limbs: 384 / LIMB_BITS,

    q: Modulus {
        p: p384_limbs![
            0xffffffff, 0x00000000, 0x00000000, 0xffffffff, 0xfffffffe, 0xffffffff, 0xffffffff,
            0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
        ],
        rr: p384_limbs![1, 0xfffffffe, 0, 2, 0, 0xfffffffe, 0, 2, 1, 0, 0, 0],
    },
    n: Elem {
        limbs: p384_limbs![
            0xccc52973, 0xecec196a, 0x48b0a77a, 0x581a0db2, 0xf4372ddf, 0xc7634d81, 0xffffffff,
            0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
        ],
        m: PhantomData,
        encoding: PhantomData, // Unencoded
    },

    a: Elem {
        limbs: p384_limbs![
            0xfffffffc, 0x00000003, 0x00000000, 0xfffffffc, 0xfffffffb, 0xffffffff, 0xffffffff,
            0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
        ],
        m: PhantomData,
        encoding: PhantomData, // Unreduced
    },
    b: Elem {
        limbs: p384_limbs![
            0x9d412dcc, 0x08118871, 0x7a4c32ec, 0xf729add8, 0x1920022e, 0x77f2209b, 0x94938ae2,
            0xe3374bee, 0x1f022094, 0xb62b21f4, 0x604fbff9, 0xcd08114b
        ],
        m: PhantomData,
        encoding: PhantomData, // Unreduced
    },

    elem_add_impl: GFp_p384_elem_add,
    elem_mul_mont: GFp_p384_elem_mul_mont,
    elem_sqr_mont: GFp_p384_elem_sqr_mont,

    point_add_jacobian_impl: GFp_nistz384_point_add,
};

pub static PRIVATE_KEY_OPS: PrivateKeyOps = PrivateKeyOps {
    common: &COMMON_OPS,
    elem_inv_squared: p384_elem_inv_squared,
    point_mul_impl: GFp_nistz384_point_mul,
};

fn p384_elem_inv_squared(a: &Elem<R>) -> Elem<R> {
    // Calculate a**-2 (mod q) == a**(q - 3) (mod q)
    //
    // The exponent (q - 3) is:
    //
    //    0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffe\
    //      ffffffff0000000000000000fffffffc

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

    let fffffffffffffff = sqr_mul(&fffffff_11, 30, &fffffff_11);

    let ffffffffffffffffffffffffffffff = sqr_mul(&fffffffffffffff, 60, &fffffffffffffff);

    // ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
    let mut acc = sqr_mul(
        &ffffffffffffffffffffffffffffff,
        120,
        &ffffffffffffffffffffffffffffff,
    );

    // fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff_111
    sqr_mul_acc(&mut acc, 15, &fff_111);

    // fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff
    sqr_mul_acc(&mut acc, 1 + 30, &fffffff_11);
    sqr_mul_acc(&mut acc, 2, &b_11);

    // fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff
    // 0000000000000000fffffff_11
    sqr_mul_acc(&mut acc, 64 + 30, &fffffff_11);

    // fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff
    // 0000000000000000fffffffc
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
        limbs: p384_limbs![
            0x333ad68c, 0x1313e696, 0xb74f5885, 0xa7e5f24c, 0x0bc8d21f, 0x389cb27e, 0, 0, 0, 0, 0,
            0
        ],

        m: PhantomData,
        encoding: PhantomData, // Unencoded
    },
};

unsafe extern "C" fn GFp_p384_elem_sqr_mont(
    r: *mut Limb,   // [COMMON_OPS.num_limbs]
    a: *const Limb, // [COMMON_OPS.num_limbs]
) {
    // XXX: Inefficient. TODO: Make a dedicated squaring routine.
    GFp_p384_elem_mul_mont(r, a, a);
}
extern "C" {
    fn GFp_p384_elem_add(
        r: *mut Limb,   // [COMMON_OPS.num_limbs]
        a: *const Limb, // [COMMON_OPS.num_limbs]
        b: *const Limb, // [COMMON_OPS.num_limbs]
    );
    fn GFp_p384_elem_mul_mont(
        r: *mut Limb,   // [COMMON_OPS.num_limbs]
        a: *const Limb, // [COMMON_OPS.num_limbs]
        b: *const Limb, // [COMMON_OPS.num_limbs]
    );

    fn GFp_nistz384_point_add(
        r: *mut Limb,   // [3][COMMON_OPS.num_limbs]
        a: *const Limb, // [3][COMMON_OPS.num_limbs]
        b: *const Limb, // [3][COMMON_OPS.num_limbs]
    );
    fn GFp_nistz384_point_mul(
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
            encoding: PhantomData,
            m: PhantomData
        },
        Scalar {
            limbs: LIMBS_ALTERNATING_10,
            encoding: PhantomData,
            m: PhantomData
        },
        Scalar {
            // n - 1
            limbs: p384_limbs![
                0xccc52973 - 1,
                0xecec196a,
                0x48b0a77a,
                0x581a0db2,
                0xf4372ddf,
                0xc7634d81,
                0xffffffff,
                0xffffffff,
                0xffffffff,
                0xffffffff,
                0xffffffff,
                0xffffffff
            ],
            encoding: PhantomData,
            m: PhantomData,
        },
    ]);
}
