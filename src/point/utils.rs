//! utils module

use crate::ring_ecc;
use ring_ecc::arithmetic::montgomery::Encoding as RingEncoding;
use ring_ecc::arithmetic::montgomery::R;
use ring_ecc::ec::suite_b::ops::{Elem,CommonOps,Scalar};
use ring_ecc::ec::suite_b::ops::elem::MAX_LIMBS;
use ring_ecc::ec::suite_b::ops::{elem_parse_big_endian_fixed_consttime,scalar_parse_big_endian_fixed_consttime};
use ring_ecc::limb::big_endian_from_limbs;
use ring_ecc::limb::LIMB_BYTES;

use subtle::ConditionallySelectable;

use num::{BigUint,Zero,One};
use num_traits::cast::ToPrimitive;
use core::marker::PhantomData;

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

pub fn biguint_to_scalar(curve_params: &CommonOps, x: &BigUint) -> Scalar {
    let x_bytes = x.to_bytes_be();
    let MAX_LEN = curve_params.num_limbs*LIMB_BYTES;
    let mut inp = vec![0; MAX_LEN];
    inp[MAX_LEN-x_bytes.len()..].copy_from_slice(&x_bytes);
    scalar_parse_big_endian_fixed_consttime(curve_params, untrusted::Input::from(&inp)).unwrap()
}

pub fn biguint_to_elem(curve_params: &CommonOps, x: &BigUint) -> Elem<R> {
    let x_bytes = x.to_bytes_be();
    let MAX_LEN = curve_params.num_limbs*LIMB_BYTES;
    let mut inp = vec![0; MAX_LEN];
    inp[MAX_LEN-x_bytes.len()..].copy_from_slice(&x_bytes);
    bytes_to_elem(curve_params, &inp)
}

pub fn elem_to_biguint<T: RingEncoding>(elem: Elem<T>) -> BigUint {
    BigUint::from_bytes_be(&elem_to_bytes(elem))
}

/// parses an elem from be_bytes
pub fn bytes_to_elem(curve_params: &CommonOps, bytes: &[u8]) -> Elem<R> {
    Elem {
        limbs: elem_parse_big_endian_fixed_consttime(curve_params, untrusted::Input::from(bytes)).unwrap().limbs,
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