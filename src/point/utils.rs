//! utils module

use ring_ecc::ec::CurveID;
use ring_ecc::arithmetic::montgomery::Encoding as RingEncoding;
use ring_ecc::arithmetic::montgomery::R;
use ring_ecc::ec::suite_b::ops::{Elem,CommonOps,Scalar,ScalarOps};
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
use std::io::Error;
use core::marker::PhantomData;

use super::{P256,P384};
use super::errors;

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

/// runs sgn_le with initial *ring* input type of `Elem<R>`
pub fn sgn0_le_elem(cops: &CommonOps, x: Elem<R>) -> i8 {
    sgn0_le(&elem_to_bytes(cops.elem_unencoded(&x)))
}

/// Returns the curve modulus (montgomery encoded) for performing modular
/// operations
pub fn get_modulus_as_biguint(cops: &CommonOps) -> BigUint {
    let p: Elem<R> = Elem {
        limbs: cops.q.p,
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
pub fn biguint_to_elem_unenc(id: CurveID, cops: &CommonOps, x: &BigUint) -> Result<Elem<R>, Error> {
    let mut bytes = x.to_bytes_be();
    let fill_len = cops.num_limbs*LIMB_BYTES;
    if bytes.len() < fill_len {
        bytes = get_filled_buffer(&x.to_bytes_be(), cops.num_limbs*LIMB_BYTES);
    }
    let input = untrusted::Input::from(&bytes);
    let res = input.read_all(Unspecified, |input| {
        match id {
            P256 => Ok(P256_PK_OPS.elem_parse(input).unwrap()),
            P384 => Ok(P384_PK_OPS.elem_parse(input).unwrap()),
            _ => Err(Unspecified)
        }
    });

    if let Err(_) = res {
        return Err(errors::internal());
    }
    Ok(res.unwrap())
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
/// TODO: remove BigInt ops
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
    biguint_to_elem_unenc(id, cops, &BigUint::one()).unwrap()
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
    biguint_to_elem_unenc(id, cops, &res).unwrap()
}

/// Returns an unencoded Scalar object
pub fn scalar_unencoded(cops: &CommonOps, sc_ops: &ScalarOps, sc: &Scalar<R>) -> Scalar {
    sc_ops.scalar_product(sc, &biguint_to_scalar(cops, &BigUint::one()))
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

pub fn extend_bytes_to_field_len(id: CurveID, bytes: &[u8]) -> Vec<u8> {
    let length = match id {
        P256 => 32,
        P384 => 48,
        _ => panic!("bad curve specified"),
    };
    let mut v = vec![0; length];
    v[length-bytes.len()..].copy_from_slice(bytes);
    v
}