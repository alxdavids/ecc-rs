//! utils module

use crate::ring_ecc;
use ring_ecc::ec::CurveID;
use ring_ecc::arithmetic::montgomery::Encoding as RingEncoding;
use ring_ecc::arithmetic::montgomery::R;
use ring_ecc::ec::suite_b::ops::{Elem,CommonOps,Scalar,Modulus};
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
use core::marker::PhantomData;

use super::{P256,P384};

/// I2osp converts an integer to a big-endian octet-string of length `xLen` and
/// copies it into `out` (see: https://tools.ietf.org/html/rfc8017#section-4.1)
pub fn I2osp(x: u64, x_len: usize, out: &mut [u8]) {
    assert!(x >= Zero::zero());
    assert!(x < 1<<(8*x_len));
    assert!(out.len() >= x_len);
    let mut val = x;
	for i in out.len()-x_len..out.len() {
        out[i] = (&val & 1) as u8;
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

pub fn sgn0_le_elem(cops: &CommonOps, x: Elem<R>) -> i8 {
    sgn0_le(&elem_to_bytes(cops.elem_unencoded(&x)))
}

/// Returns the curve modulus (montgomery encoded) for performing modular
/// operations
pub fn get_modulus_as_biguint(q: &Modulus) -> BigUint {
    let p: Elem<R> = Elem {
        limbs: q.p,
        m: PhantomData,
        encoding: PhantomData
    };
    elem_to_biguint(p)
}

pub fn biguint_to_scalar(cops: &CommonOps, x: &BigUint) -> Scalar {
    let inp = get_filled_buffer(&x.to_bytes_be(), cops.num_limbs*LIMB_BYTES);
    scalar_parse_big_endian_fixed_consttime(cops, untrusted::Input::from(&inp)).unwrap()
}

/// Translates a value into a *ring* Elem<R>, where the value has not previously
/// been encoded
///
/// TODO: come up with better error type!!
pub fn biguint_to_elem_unenc(id: CurveID, cops: &CommonOps, x: &BigUint) -> Elem<R> {
    let mut bytes = x.to_bytes_be();
    let fill_len = cops.num_limbs*LIMB_BYTES;
    if bytes.len() < fill_len {
        bytes = get_filled_buffer(&x.to_bytes_be(), cops.num_limbs*LIMB_BYTES);
    }
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
    bytes_to_elem(cops, &get_filled_buffer(&x.to_bytes_be(), cops.num_limbs*LIMB_BYTES))
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

pub fn elem_one(id: CurveID, cops: &CommonOps) -> Elem<R> {
    biguint_to_elem_unenc(id, cops, &BigUint::one())
}

/// returns a new vector that is filled with 0's up to `fill_len` length
fn get_filled_buffer(unfilled: &[u8], fill_len: usize) -> Vec<u8> {
    let mut wrapper = vec![0; fill_len];
    wrapper[fill_len-unfilled.len()..].copy_from_slice(unfilled);
    wrapper
}

/// Performs a square root operation in the underlying field using the input
/// exponent
pub fn sqrt(input: &BigUint, exp: &BigUint, modulus: &BigUint) -> BigUint {
    input.modpow(exp, modulus)
}

/// Performs a square root operation in the underlying field using the input
/// exponent
pub fn elem_sqrt(id: CurveID, cops: &CommonOps, elem: Elem<R>, exp: &BigUint, modulus: &BigUint) -> Elem<R> {
    let input = elem_to_biguint(cops.elem_unencoded(&elem));
    let res = sqrt(&input, exp, modulus);
    biguint_to_elem_unenc(id, cops, &res)
}

/// Performs an exponentiation in the underlying field that ascertains whether
/// the current element is a square, or not.
pub fn is_square(input: &BigUint, exp: &BigUint, modulus: &BigUint) -> bool {
    let res = input.modpow(exp, modulus);
    let c = res == Zero::zero();
    let d = res == One::one();
    c || d
}

/// Returns `a` if `c` is `false`, and `b` if c is `true`
///
/// TODO: remove BigUint operations
pub fn elem_cmov(cops: &CommonOps, a: Elem<R>, b: Elem<R>, c: bool) -> Elem<R> {
    let c_bu = BigUint::from(c as u8);
    let a_bu = elem_to_biguint(a);
    let b_bu = elem_to_biguint(b);
    let res = (&c_bu * b_bu) + ((BigUint::one() - &c_bu) * a_bu);
    biguint_to_elem(cops, &res)
}

/// Moves the contents of `src` into the provided output buffer `dst`. Clears
/// the contents of `dst` first.
pub fn copy_into_cleared(src: &[u8], dst: &mut Vec<u8>) {
    dst.clear();
    dst.extend_from_slice(src)
}

pub mod hkdf {
    //! The hkdf module provides access to the functionality provoided HKDF as
    //! specified in [RFC5869](https://tools.ietf.org/html/rfc5869)
    //!
    //! A wrapper around the rust-crypto implementation of HKDF
    //! [RFC5869](https://tools.ietf.org/html/rfc5869) locked to using SHA512 as the
    //! underlying hash function. The abstraction is necessary so that we can
    //! provide multiple sized output buffers, where the HKDF instance.
    //!
    //! TODO: Rewrite to use the ring implementation. There were some difficulties
    //! around the way that ring does not give access to the raw bytes output by
    //! these algorithms

    use crypto::hkdf::{hkdf_extract,hkdf_expand};
    use crypto::sha2::Sha512;
    use super::copy_into_cleared;

    const SHA512_OUTPUT_BYTES_LENGTH: usize = 64;

    /// runs HKDF_Extract as specified in
    /// [RFC5869](https://tools.ietf.org/html/rfc5869).
    pub fn extract(seed: &[u8], secret: &[u8], out: &mut Vec<u8>) {
        if out.len() != SHA512_OUTPUT_BYTES_LENGTH {
            copy_into_cleared(&[0; SHA512_OUTPUT_BYTES_LENGTH], out);
        }
        hkdf_extract(Sha512::new(), &seed, &secret, out)
    }

    /// runs HKDF_Expand as specified in
    /// [RFC5869](https://tools.ietf.org/html/rfc5869). The value of `prk`
    /// should be uniformly sampled bytes
    pub fn expand(prk: &[u8], info: &[u8], out: &mut Vec<u8>) {
        hkdf_expand(Sha512::new(), &prk, &info, out)
    }
}