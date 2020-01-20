//! utils module

use crate::ring_ecc;
use ring_ecc::ec::CurveID;
use ring_ecc::arithmetic::montgomery::Encoding as RingEncoding;
use ring_ecc::arithmetic::montgomery::R;
use ring_ecc::ec::suite_b::ops::{Elem,CommonOps,Scalar};
use ring_ecc::ec::suite_b::ops::PrivateKeyOps as CurveOps;
use ring_ecc::ec::suite_b::ops::p256::PUBLIC_KEY_OPS as P256_PK_OPS;
use ring_ecc::ec::suite_b::ops::p384::PUBLIC_KEY_OPS as P384_PK_OPS;
use ring_ecc::ec::suite_b::ops::elem::MAX_LIMBS;
use ring_ecc::ec::suite_b::ops::{elem_parse_big_endian_fixed_consttime,scalar_parse_big_endian_fixed_consttime};
use ring_ecc::limb::big_endian_from_limbs;
use ring_ecc::limb::LIMB_BYTES;
use ring_ecc::error::Unspecified;

use untrusted;
use subtle::ConditionallySelectable;

use num::{BigUint,Zero,One};
use num_traits::cast::ToPrimitive;
use core::marker::PhantomData;

use super::{P256,P384};

/// I2osp converts an integer to a big-endian octet-string of length `xLen` and
/// copies it into `out` (see: https://tools.ietf.org/html/rfc8017#section-4.1)
pub fn I2osp(x: BigUint, x_len: usize, out: &mut [u8]) {
    assert!(x >= Zero::zero());
    assert!(x < BigUint::from((1<<(8*x_len)) as u64));
    assert!(out.len() >= x_len);
    let mut val = x;
	for i in out.len()-x_len..out.len() {
        out[i] = (&val & BigUint::one()).to_u8().unwrap();
		val = &val >> 8;
    }
}

/// Os2ip converts an octet-string to a BigUint object
/// (https://tools.ietf.org/html/rfc8017#section-4.1)
pub fn Os2ip(x: &[u8]) -> BigUint {
	BigUint::from_bytes_be(x)
}

/// returns `1_i8` if the sign is positive, and `-1_i8` if it is negative. As
/// documented in draft-irtf-cfrg-hash-to-curve (Section
/// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-05#section-4.1.2).
/// Expects the input `x` to be in big-endian format.
pub fn sgn0_le(x: &[u8]) -> i8 {
    let pos = x[x.len()-1];
    let res = pos & 1;
    let sgn = i8::conditional_select(&1, &(-1), res.into());
    let zero_cmp = (pos == 0) as u8;
    let sgn = i8::conditional_select(&sgn, &0, zero_cmp.into());
    let zero_cmp_2 = (sgn == 0) as u8;
    i8::conditional_select(&sgn, &1, zero_cmp_2.into())
}

pub fn biguint_to_scalar(cops: &CommonOps, x: &BigUint) -> Scalar {
    let x_bytes = x.to_bytes_be();
    let MAX_LEN = cops.num_limbs*LIMB_BYTES;
    let mut inp = vec![0; MAX_LEN];
    inp[MAX_LEN-x_bytes.len()..].copy_from_slice(&x_bytes);
    scalar_parse_big_endian_fixed_consttime(cops, untrusted::Input::from(&inp)).unwrap()
}

/// Translates a value into a *ring* Elem<R>, where the value has not previously
/// been encoded
///
/// TODO: come up with better error type!!
pub fn biguint_to_elem_unenc(id: CurveID, x: &BigUint) -> Elem<R> {
    let bytes = x.to_bytes_be();
    let input = untrusted::Input::from(&bytes);
    input.read_all(Unspecified, |input| {
        Ok(
            match id {
                P256 => P256_PK_OPS.elem_parse(input).unwrap(),
                P384 => P384_PK_OPS.elem_parse(input).unwrap(),
                _ => panic!("blah")
            }
        )
    }).unwrap()
}

pub fn biguint_to_elem(cops: &CommonOps, x: &BigUint) -> Elem<R> {
    let x_bytes = x.to_bytes_be();
    let MAX_LEN = cops.num_limbs*LIMB_BYTES;
    let mut inp = vec![0; MAX_LEN];
    inp[MAX_LEN-x_bytes.len()..].copy_from_slice(&x_bytes);
    bytes_to_elem(cops, &inp)
}

pub fn elem_to_biguint<T: RingEncoding>(elem: Elem<T>) -> BigUint {
    BigUint::from_bytes_be(&elem_to_bytes(elem))
}

/// parses an elem from be_bytes
pub fn bytes_to_elem(cops: &CommonOps, bytes: &[u8]) -> Elem<R> {
    Elem {
        limbs: elem_parse_big_endian_fixed_consttime(cops, untrusted::Input::from(bytes)).unwrap().limbs,
        encoding: PhantomData,
        m: PhantomData
    }
}

/// transform elem to be_bytes
pub fn elem_to_bytes<T: RingEncoding>(elem: Elem<T>) -> Vec<u8> {
    let mut out: Vec<u8> = vec![0; MAX_LIMBS*LIMB_BYTES];
    big_endian_from_limbs(&elem.limbs, &mut out);
    out
}

/// gets the negative value of the element in the field
///
/// TODO: uses BigInt ops
pub fn minus_elem(cops: &CommonOps, p: &BigUint, elem: Elem<R>) -> Elem<R> {
    let elem_bu = elem_to_biguint(elem);
    biguint_to_elem(cops, &(p - elem_bu))
}

/// computes the element representing 1/x in the base field
pub fn invert_elem(ops: &CurveOps, x: Elem<R>) -> Elem<R> {
    let xx_inv = ops.elem_inverse_squared(&x);
    ops.common.elem_product(&xx_inv, &x)
}

/// Performs a square root operation in the underlying field using the input
/// exponent
pub fn sqrt(input: &BigUint, exp: &BigUint, modulus: &BigUint) -> BigUint {
    input.modpow(exp, modulus)
}

/// Performs an exponentiation in the underlying field that ascertains whether
/// the current element is a square, or not.
pub fn is_square(input: &BigUint, exp: &BigUint, modulus: &BigUint) -> BigUint {
    input.modpow(exp, modulus)
}