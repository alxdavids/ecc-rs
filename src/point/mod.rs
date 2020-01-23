//! The point module provides the main interface for performing elliptic curve
//! operations. Currently we only support the `secp256r1` and `secp384r1` curves,
//! both of which are implemented in the `suite_b` module of ring.
//!
//! # Examples
//!
//! Performing operations on elliptic curve points is as easy as running:
//!
//! ```
//! use ecc_rs::point::*;
//! use num::BigUint;
//!
//! // get new zero point
//! let _ = AffinePoint::new(P256);
//!
//! // get generator point
//! let g = AffinePoint::get_generator(P256);
//!
//! // perform scalar mult
//! let k = BigUint::parse_bytes(b"0a", 16).unwrap();
//! let k_base = AffinePoint::base_mul(P256, &k.to_bytes_be());
//! let k_g = g.scalar_mul(&k.to_bytes_be());
//!
//! // check validity and equality
//! assert_eq!(k_g.is_valid(), true);
//! assert_eq!(k_base.is_valid(), true);
//! assert_eq!(k_g.equals(&k_base), true);
//!
//! // add points
//! let k_1_g = g.to_jacobian().add(&k_g).to_affine();
//!
//! // check point is valid
//! // assert_eq!(k_1_g.is_valid(), true);
//!
//! // check points are equal
//! let k_1 = BigUint::parse_bytes(b"0b", 16).unwrap();
//! let k_1_chk_g = g.scalar_mul(&k_1.to_bytes_be()).to_affine();
//! // assert_eq!(k_1_g.equals(&k_1_chk_g), true);
//! ```
//!
//! Conversion between affine and jacobian points is also simple:
//!
//! ```
//! use ecc_rs::point::*;
//! let gen = AffinePoint::get_generator(P256);
//! let jac = gen.to_jacobian();
//! let gen_conv = jac.to_affine();
//! assert_eq!(gen.is_valid(), true);
//! assert_eq!(gen_conv.is_valid(), true);
//! assert_eq!(jac.is_valid(), true);
//! ```
//!
//! Point serialization/deserialization can be done for sending point data
//! across the wire. Note that serialization of a point converts into unencoded
//! format, and deserialization performs re-encoding internally:
//!
//! ```
//! use ecc_rs::point::*;
//! let gen = AffinePoint::get_generator(P256);
//!
//! // uncompressed
//! let ser_decmp = gen.serialize(false);
//! let deser_decmp = (AffinePoint::new(P256)).deserialize(&ser_decmp);
//!
//! // compressed
//! let ser_cmp = gen.serialize(true);
//! let deser_cmp = (AffinePoint::new(P256)).deserialize(&ser_cmp);
//! ```
//!
//! Get deterministic (randomly distributed) curve points using arbitrary input
//! bytes (see https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-05 for
//! full spec):
//!
//! ```
//! use ecc_rs::point::*;
//! let p = AffinePoint::new(P256);
//! let r = p.hash_to_curve("some_input".as_bytes());
//! assert!(r.is_valid());
//! ```

use num::{BigUint,BigInt,Zero};
use num_bigint::ToBigInt;
use core::marker::PhantomData;

use subtle::ConditionallySelectable;

use crate::ring_ecc;
use ring_ecc::ec::CurveID;
use ring_ecc::ec::suite_b::verify_affine_point_is_on_the_curve;
use ring_ecc::ec::suite_b::ops::Point as RingPoint;
use ring_ecc::ec::suite_b::ops::PrivateKeyOps as CurveOps;
use ring_ecc::ec::suite_b::ops::Elem;
use ring_ecc::ec::suite_b::ops::p256::PRIVATE_KEY_OPS as P256_OPS;
use ring_ecc::ec::suite_b::ops::p384::PRIVATE_KEY_OPS as P384_OPS;
use ring_ecc::ec::suite_b::ops::p256::PUBLIC_SCALAR_OPS as P256_SCALAR_OPS;
use ring_ecc::ec::suite_b::ops::p384::PUBLIC_SCALAR_OPS as P384_SCALAR_OPS;
use ring_ecc::ec::suite_b::ops::p384::P384_GENERATOR;
use ring_ecc::limb::{Limb,LIMB_BYTES};
use ring_ecc::ec::suite_b::ops::elem::MAX_LIMBS;
use ring_ecc::ec::suite_b::private_key::affine_from_jacobian;
use ring_ecc::arithmetic::montgomery::R;
use ring_ecc::rand;
use rand::sealed::SecureRandom;

mod h2c;
mod utils;

/// P256 constant for exposing *ring* CurveID::P256 enum
pub const P256: CurveID = CurveID::P256;
/// P384 constant for exposing *ring* CurveID::P384 enum
pub const P384: CurveID = CurveID::P384;

/// Used for indicating that a point is in montgomery encoding. This encoding is
/// used for perfoming all group operations since this is how *ring* internals
/// do it.
pub enum Encoded {}
/// Indiciates that the point is unencoded, aka transmission format.
pub enum Unencoded {}

/// The `AffinePoint` struct provides access to a base elliptic curve point
/// representation. The setting of `id` and `ops` determine which curve the
/// point is associated with. Currently we support `secp256r1` and `secp384r1`
/// as this is what is supported in ring.
///
/// The `AffinePoint` object can either be in montgomery encoding (`T =
/// Encoded`), or in unencoded (`T = Unencoded`) format. All operations are
/// performed on points in montgomery encoding as this is how *ring* does it.
/// Operations usually convert a point into Jacobian format, and `to_affine()`
/// should be called to recover the affine point respresentation. Unencoded
/// points are essentially private structs that are used when serializing and
/// deserializing the points.
///
/// # Examples
///
/// The best way to recover an `AffinePoint` object is to get the generator for
/// the requested curve (either `P256` or `P384`):
/// ```
/// use ecc_rs::point::*;
/// let g = AffinePoint::get_generator(P256);
/// ```
pub struct AffinePoint<T> {
    x: BigUint,
    y: BigUint,
    id: CurveID,
    ops: &'static CurveOps,
    encoding: PhantomData<T>,
}

impl AffinePoint<Encoded> {
    /// Returns a zero point, this isn't really a valid point representation but
    /// can be useful when constructing points via point addition
    pub fn new(id: CurveID) -> Self {
        match id {
            P256 => Self {
                x: Zero::zero(),
                y: Zero::zero(),
                id: P256,
                ops: &P256_OPS,
                encoding: PhantomData,
            },
            P384 => Self {
                x: Zero::zero(),
                y: Zero::zero(),
                id: P384,
                ops: &P384_OPS,
                encoding: PhantomData,
            },
            _ => panic!("unsupported"),
        }
    }

    /// Returns the fixed generator of the group instantiated by the curve (in
    /// montgomery encoding)
    pub fn get_generator(id: CurveID) -> Self {
        match id {
            P256 => Self {
                // *ring* doesn't explicitly make the P256 generator available,
                // but some applications need it. we use the montgomery encoded
                // values provided in (except with -y coordinate)
                // ring_ecc/ec/suite_b/ops/p256_point_mul_tests.txt
                x: BigUint::parse_bytes(b"18905f76a53755c679fb732b7762251075ba95fc5fedb60179e730d418a9143c", 16).unwrap(),
                y: BigUint::parse_bytes(b"8571ff1825885d85d2e88688dd21f3258b4ab8e4ba19e45cddf25357ce95560a", 16).unwrap(),
                id: P256,
                ops: &P256_OPS,
                encoding: PhantomData,
            },
            P384 => {
                // for some reason the *ring* generator does not agree with
                // https://www.secg.org/SEC2-Ver-1.0.pdf Section 2.8.1. using
                // the *ring* point for now
                Self {
                    x: utils::elem_to_biguint(P384_GENERATOR.0),
                    y: utils::elem_to_biguint(P384_GENERATOR.1),
                    id: P384,
                    ops: &P384_OPS,
                    encoding: PhantomData,
                }
            },
            _ => panic!("unsupported"),
        }
    }

    /// Performs `k*G`, where `G` is the fixed generator of the group, and `k`
    /// is the input scalar
    pub fn base_mul(id: CurveID, scalar: &[u8]) -> JacobianPoint<Encoded> {
        AffinePoint::get_generator(id).scalar_mul(scalar)
    }

    /// Performs `k*P` where `P` is some point on the curve and `k` is the input
    /// scalar
    pub fn scalar_mul(&self, scalar: &[u8]) -> JacobianPoint<Encoded> {
        let sc_bu = BigUint::from_bytes_be(scalar);
        let ring_sc = utils::biguint_to_scalar(self.ops.common, &sc_bu);
        let x_elem = utils::biguint_to_elem(self.ops.common, &self.x);
        let y_elem = utils::biguint_to_elem(self.ops.common, &self.y);
        let ring_pt = self.ops.point_mul(&ring_sc, &(x_elem, y_elem));
        JacobianPoint::from_ring_jac_point(ring_pt, self.id, self.ops)
    }

    /// Returns `true` if the two points are equal (share the same coordiantes
    /// and belong to the same curve), returns false otherwise
    pub fn equals(&self, other: &Self) -> bool {
        (self.x == other.x)
        && (self.y == other.y)
        && (self.id == other.id)
    }

    /// Returns `true` if the point satisfies the equation of the curve that it
    /// belongs to, returning `false` otherwise
    pub fn is_valid(&self) -> bool {
        let (x, y) = self.as_ring_affine();
        match verify_affine_point_is_on_the_curve(self.ops.common, (&x, &y)) {
            Ok(_) => true,
            Err(_) => false,
        }
    }

    /// Converts the `AffinePoint` object into a `JacobianPoint` object (i.e.
    /// representing it in jacobian coordinates)
    ///
    /// TODO: Do this in a less naive way
    pub fn to_jacobian(&self) -> JacobianPoint<Encoded> {
        self.scalar_mul(&vec![1_u8])
    }

    /// Converts the point into either uncompressed or compressed format
    /// depending on the value of the bool `compress`. The serialization of the
    /// point follows the format specified in
    /// [SEC1](https://www.secg.org/sec1-v2.pdf).
    pub fn serialize(&self, compress: bool) -> Vec<u8> {
        let p = self.to_unencoded();
        let mut out: Vec<u8> = Vec::new();
        let x_bytes = p.x.to_bytes_be();
        let y_bytes = p.y.to_bytes_be();
        if !compress {
            out.push(0x04);
            out.extend_from_slice(&x_bytes);
            out.extend_from_slice(&y_bytes);
        } else {
            let sign = utils::sgn0_le(&y_bytes);
            // add tag for compression
            match sign {
                1 => out.push(0x02),
                -1 => out.push(0x03),
                _ => panic!("not intended!")
            };
            out.extend_from_slice(&x_bytes);
        }
        out
    }

    /// Attempts to convert the provided byte array into a new encoded
    /// `AffinePoint` object. Fails if the recovered point is not valid.
    ///
    /// When a point is 'uncompressed' it should take the form `[0x04, <--
    /// coord_byte_length -->, <-- coord_byte_length -->]`. When it is
    /// compressed it should take the form `[b, <-- coord_byte_length -->]`,
    /// where `b == 0x02` or `b == 0x03`.
    ///
    /// # Warning
    ///
    /// Point deserialization for compressed points currently uses
    /// BigUint arithmetic and so it cannot be considered constant-time (even
    /// though it is implemented using straight-line arithmetic).
    pub fn deserialize(&mut self, input: &[u8]) -> AffinePoint<Encoded> {
        let coord_byte_length = self.ops.common.num_limbs*LIMB_BYTES;
        let cops = self.ops.common;
        let x = BigUint::from_bytes_be(&input[1..coord_byte_length+1]);
        let parity = input[0];
        // An array of bytes is in uncompressed format if the first byte of
        // `input` is `0x04`, and in compressed format if it is `0x02` or
        // `0x03`.
        let y = match parity {
            0x04 => {
                assert_eq!(input.len(), 2*coord_byte_length+1);
                BigUint::from_bytes_be(&input[1+coord_byte_length..])
            },
            0x02 | 0x03 => {
                assert_eq!(input.len(), coord_byte_length+1);
                let x_enc = utils::biguint_to_elem_unenc(self.id, cops, &x);
                let xx_enc = self.ops.common.elem_squared(&x_enc);
                let scalar_ops = match self.id {
                    P256 => &P256_SCALAR_OPS,
                    P384 => &P384_SCALAR_OPS,
                    _ => panic!("bad id")
                };
                let xx__a_unenc = scalar_ops.elem_sum(
                    &cops.elem_unencoded(&xx_enc),
                    &cops.elem_unencoded(&self.ops.common.a)
                );
                let xxx__ax_unenc = cops.elem_product(&x_enc, &xx__a_unenc);
                let yy_unenc = scalar_ops.elem_sum(&xxx__ax_unenc,
                                                &cops.elem_unencoded(&cops.b)
                                            );
                // perform square root
                // (TODO: don't do BigUint ops as they are not constant time)
                let q = utils::get_modulus_as_biguint(&self.ops.common);
                let sqrt_exp = (&q+BigUint::from(1_u64))/BigUint::from(4_u64);
                let y_sqrt = utils::sqrt(&utils::elem_to_biguint(yy_unenc), &sqrt_exp, &q);
                let y_sqrt_bytes = y_sqrt.to_bytes_be();

                // Check that the sign matches what is sent in the compressed
                // encoding. If not we need to get the alternative point with -y
                // coordinate, where y is the recovered coordinate
                let y_sqrt_sign = utils::sgn0_le(&y_sqrt_bytes);
                let y_cmp = (y_sqrt_sign == -1) as u8;
                let parity_bit = (0x03 == parity) as u8;
                let parity_cmp = (parity_bit == y_cmp) as u8;
                let out_sign = i8::conditional_select(&1, &(-1), parity_cmp.into());
                let mod_mask = u8::conditional_select(&1, &0, parity_cmp.into());
                let y_par = BigInt::from(out_sign) * y_sqrt.to_bigint().unwrap();
                let modulus = (BigUint::from(mod_mask) * &q).to_bigint().unwrap();
                (modulus - y_par).to_biguint().unwrap()
            },
            _ => panic!("invalid point serialization")
        };

        // returns an encoded point
        (AffinePoint {
            x: x,
            y: y,
            id: self.id,
            ops: self.ops,
            encoding: PhantomData
        }).to_encoded()
    }

    /// Hashes the input bytes `alpha` deterministically to a random point on
    /// the curve. Ensures that the output curve point does not reveal the
    /// discrete log with respect to the fixed group generator.
    pub fn hash_to_curve(&self, alpha: &[u8]) -> Self {
        let h2c = h2c::HashToCurve::new(self.id, self.ops, utils::get_modulus_as_biguint(self.ops.common));
        h2c.full(alpha)
    }

    /// Returns a uniform number of bytes equal to the length of elements in the
    /// base field
    pub fn uniform_bytes_from_field(&self) -> Vec<u8> {
        let fill_len = self.ops.common.num_limbs*LIMB_BYTES;
        let mut out = vec![0; fill_len];
        let rng = rand::SystemRandom::new();
        if let Err(_) = rng.fill_impl(&mut out) {
            panic!("error filling")
        }
        out
    }

    /// Returns the `AffinePoint` object in unencoded format, removing the
    /// montgomery encoding
    fn to_unencoded(&self) -> AffinePoint<Unencoded> {
        let ops = self.ops.common;
        let (x_r, y_r) = self.as_ring_affine();
        let (x_unenc, y_unenc) = (ops.elem_unencoded(&x_r), ops.elem_unencoded(&y_r));
        AffinePoint {
            x: utils::elem_to_biguint(x_unenc),
            y: utils::elem_to_biguint(y_unenc),
            id: self.id,
            ops: self.ops,
            encoding: PhantomData,
        }
    }

    /// Returns the affine coordinates of the point in the format used by
    /// *ring*. This is used when performing operations using the *ring*
    /// internals
    fn as_ring_affine(&self) -> (Elem<R>, Elem<R>) {
        (utils::biguint_to_elem(self.ops.common, &self.x), utils::biguint_to_elem(self.ops.common, &self.y))
    }
}

impl AffinePoint<Unencoded> {
    /// Re-encodes the `AffinePoint` object in montgomery encoding, for
    /// performing group operations.
    fn to_encoded(&self) -> AffinePoint<Encoded> {
        AffinePoint {
            x: utils::elem_to_biguint(utils::biguint_to_elem_unenc(self.id, self.ops.common, &self.x)),
            y: utils::elem_to_biguint(utils::biguint_to_elem_unenc(self.id, self.ops.common, &self.y)),
            id: self.id,
            ops: self.ops,
            encoding: PhantomData,
        }
    }
}

/// The `JacobianPoint` struct represents an elliptic curve point in jacobian
/// coordinates (in either montgomery encoding or not).
pub struct JacobianPoint<T> {
    x: BigUint,
    y: BigUint,
    z: BigUint,
    id: CurveID,
    ops: &'static CurveOps,
    encoding: PhantomData<T>
}

impl JacobianPoint<Encoded> {
    /// Returns the identity point
    pub fn new(id: CurveID) -> Self {
        AffinePoint::new(id).to_jacobian()
    }

    /// Performs ellptic curve point addition on the point `self` and `other`
    /// (both in jacobian coordinates) and returns the result as a new point
    pub fn add(&self, other: &JacobianPoint<Encoded>) -> Self {
        let p1 = self.as_ring_jac_point();
        let p2 = other.as_ring_jac_point();
        let p_add = self.ops.common.point_sum(&p1, &p2);
        Self::from_ring_jac_point(p_add, self.id, self.ops)
    }

    /// Returns `true` if the coordinates for `self` and `other` are the same,
    /// and the points are on the same curve. Returns `false` otherwise.
    pub fn equals(&self, other: &Self) -> bool {
        (self.x == other.x)
        && (self.y == other.y)
        && (self.z == other.z)
        && (self.id == other.id)
    }

    /// Converts the point to an `AffinePoint` object and runs `is_valid()`.
    pub fn is_valid(&self) -> bool {
        self.to_affine().is_valid()
    }

    /// Converts the `JacobianPoint` object to an `AffinePoint` object (i.e.
    /// affine coordinates) and returns the result.
    pub fn to_affine(&self) -> AffinePoint<Encoded> {
        let (x_ele, y_ele) = affine_from_jacobian(self.ops, &self.as_ring_jac_point()).unwrap();
        AffinePoint {
            x: utils::elem_to_biguint(x_ele),
            y: utils::elem_to_biguint(y_ele),
            id: self.id,
            ops: self.ops,
            encoding: PhantomData
        }
    }

    /// Converts the coordinates of the point to the representation used in
    /// *ring* internals
    fn as_ring_jac_point(&self) -> RingPoint {
        let (x, y, z) = self.use_for_ring_ops();
        let mut limbs: [Limb; MAX_LIMBS*3] = [0; MAX_LIMBS*3];
        limbs[..MAX_LIMBS].copy_from_slice(&x.limbs);
        limbs[self.ops.common.num_limbs..self.ops.common.num_limbs+MAX_LIMBS].copy_from_slice(&y.limbs);
        limbs[self.ops.common.num_limbs*2..self.ops.common.num_limbs*2+MAX_LIMBS].copy_from_slice(&z.limbs);
        RingPoint { xyz: limbs }
    }

    /// Returns a ring point representation back into a `JacobianPoint` struct
    fn from_ring_jac_point(ring_pt: RingPoint, id: CurveID, ops: &'static CurveOps) -> Self {
        let x = utils::elem_to_biguint(ops.common.point_x(&ring_pt));
        let y = utils::elem_to_biguint(ops.common.point_y(&ring_pt));
        let z = utils::elem_to_biguint(ops.common.point_z(&ring_pt));
        Self {
            x: x,
            y: y,
            z: z,
            id: id,
            ops: ops,
            encoding: PhantomData
        }
    }

    /// encodes the fields of the `JacobianPoint` into the internal *ring*
    /// format
    fn use_for_ring_ops(&self) -> (Elem<R>, Elem<R>, Elem<R>) {
        (
            utils::biguint_to_elem(self.ops.common, &self.x),
            utils::biguint_to_elem(self.ops.common, &self.y),
            utils::biguint_to_elem(self.ops.common, &self.z),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hex;

    #[test]
    fn scalar_mul_test() {
        for i in 0..2 {
            let gen = match i {
                0 => AffinePoint::get_generator(P256),
                1 => AffinePoint::get_generator(P384),
                _ => panic!("test failed")
            };

            // for some reason the *ring* test vectors use their bespoke generator
            // with the minus y coordinate
            fn get_minus_y_point(point: &AffinePoint<Encoded>) -> AffinePoint<Encoded> {
                let cops = point.ops.common;
                AffinePoint {
                    x: point.x.clone(),
                    y: utils::elem_to_biguint(utils::minus_elem(cops, &utils::get_modulus_as_biguint(cops), utils::biguint_to_elem(cops, &point.y))),
                    id: point.id,
                    ops: point.ops,
                    encoding: PhantomData,
                }
            }
            let test_gen = get_minus_y_point(&gen);
            for vector in MULT_TEST_VECTORS[i].iter() {
                let k = BigUint::parse_bytes(vector[0].as_bytes(), 16).unwrap();
                let aff = test_gen.scalar_mul(&k.to_bytes_be()).to_affine();

                // hex encode padding 0's (to_str_radix doesn't do this)
                let x_hex = hex::encode(aff.x.to_bytes_be());
                let y_hex = hex::encode(aff.y.to_bytes_be());

                assert_eq!(x_hex, vector[1], "x hex coordinate check for scalar_mult on curve: {:?} and scalar: {}", gen.id, vector[0]);
                assert_eq!(y_hex, vector[2], "y hex coordinate check for scalar_mult on curve: {:?} and scalar: {}", gen.id, vector[0]);
            }
        }
    }

    #[test]
    fn add_test() {
        for i in 0..2 {
            // TODO use generators
            let (ops, id) = match i {
                0 => (&P256_OPS, P256),
                1 => (&P384_OPS, P384),
                _ => panic!("test failed")
            };
            let vectors = ADD_TEST_VECTORS[i];
            for j in 0..vectors.len() {
                let pt_left: JacobianPoint<Encoded> = JacobianPoint {
                    x: BigUint::parse_bytes(vectors[j][0][0].as_bytes(), 16).unwrap(),
                    y: BigUint::parse_bytes(vectors[j][0][1].as_bytes(), 16).unwrap(),
                    z: BigUint::parse_bytes(vectors[j][0][2].as_bytes(), 16).unwrap(),
                    id: id,
                    ops: ops,
                    encoding: PhantomData,
                };
                let pt_right: JacobianPoint<Encoded> = JacobianPoint {
                    x: BigUint::parse_bytes(vectors[j][1][0].as_bytes(), 16).unwrap(),
                    y: BigUint::parse_bytes(vectors[j][1][1].as_bytes(), 16).unwrap(),
                    z: BigUint::parse_bytes(vectors[j][1][2].as_bytes(), 16).unwrap(),
                    id: id,
                    ops: ops,
                    encoding: PhantomData,
                };

                let aff = pt_left.add(&pt_right).to_affine();
                let x_hex = hex::encode(aff.x.to_bytes_be());
                let y_hex = hex::encode(aff.y.to_bytes_be());

                assert_eq!(aff.is_valid(), true);
                assert_eq!(x_hex, vectors[j][2][0], "x hex coordinate check for addition curve: {:?} and j: {}", id, j);
                assert_eq!(y_hex, vectors[j][2][1], "y hex coordinate check for addition curve: {:?} and j: {}", id, j);
            }
        }
    }

    #[test]
    fn validity_test() {
        let k = "0a".as_bytes();
        for &id in [P256, P384].iter() {
            let k_base = AffinePoint::base_mul(id, &k);
            let k_g = AffinePoint::get_generator(id).scalar_mul(&k);
            assert_eq!(k_g.is_valid(), true, "for id: {:?}", id);
            assert_eq!(k_base.is_valid(), true, "for id: {:?}", id);
            assert_eq!(k_g.equals(&k_base), true, "for id: {:?}", id);
        }
    }

    #[test]
    fn point_conversion() {
        for &id in [P256,P384].iter() {
            let gen = AffinePoint::get_generator(id);
            let jac = gen.to_jacobian();
            let gen_conv = jac.to_affine();
            assert_eq!(gen.is_valid(), true, "for id: {:?}", id);
            assert_eq!(gen_conv.is_valid(), true, "for id: {:?}", id);
            assert_eq!(jac.is_valid(), true, "for id: {:?}", id);
        }
    }

    #[test]
    fn point_encoding() {
        for &id in [P256,P384].iter() {
            let gen = AffinePoint::get_generator(id);
            let enc = gen.to_unencoded();
            // check unencoded point representation
            match id {
                P256 => {
                    assert_eq!(hex::encode(enc.x.to_bytes_be()), "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296");
                    assert_eq!(hex::encode(enc.y.to_bytes_be()), "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5");
                },
                P384 => {
                    assert_eq!(hex::encode(enc.x.to_bytes_be()), "aa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7");
                    assert_eq!(hex::encode(enc.y.to_bytes_be()), "3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f");
                },
                _ => panic!("bad curve chosen")
            };
            let unenc = enc.to_encoded();
            assert_eq!(unenc.is_valid(), true, "for id: {:?}", id);
            assert_eq!(unenc.equals(&gen), true, "for id: {:?}", id);
        }
    }

    #[test]
    fn point_serialization() {
        for &id in [P256,P384].iter() {
            let gen = AffinePoint::get_generator(id);
            let ser_de = gen.serialize(false);
            let deser_de = (AffinePoint::new(id)).deserialize(&ser_de);
            assert_eq!(deser_de.is_valid(), true, "decompressed point validity check for {:?}", id);
            assert_eq!(deser_de.equals(&gen), true, "decompressed point equality check for {:?}", id);
        }
    }

    #[test]
    fn point_serialization_compressed() {
        for &id in [P256,P384].iter() {
            let gen = AffinePoint::get_generator(id);
            let ser_cmp = gen.serialize(true);
            let deser_cmp = (AffinePoint::new(id)).deserialize(&ser_cmp);
            assert_eq!(deser_cmp.is_valid(), true, "compressed point validity check for {:?}", id);
            assert_eq!(deser_cmp.equals(&gen), true, "compressed point equality check for {:?}", id);
        }
    }

    #[test]
    fn hash_to_curve() {
        for &id in [P256,P384].iter() {
            let p = AffinePoint::new(id);
            let r = p.hash_to_curve("some_input".as_bytes());
            assert!(r.is_valid(), "curve: {:?}", id);
        }
    }

    #[test]
    fn uniform_bytes() {
        for &id in [P256,P384].iter() {
            let p = AffinePoint::new(id);
            let out = p.uniform_bytes_from_field();
            assert_eq!(out.len(), p.ops.common.num_limbs*LIMB_BYTES, "curve: {:?}", id);
        }
    }

    // taken from test vectors at ring_ecc/ec/suite_b/ops/
    const MULT_TEST_VECTORS: [[[&str; 3]; 10]; 2] = [
        // P-256 test vectors
        [
            ["1", "18905f76a53755c679fb732b7762251075ba95fc5fedb60179e730d418a9143c", "7a8e00e6da77a27b2d17797722de0cda74b5471c45e61ba3220daca8316aa9f5"],
            ["2", "f6bb32e43dcf3a3b732205038d1490d9aa6ae3c1a433827d850046d410ddd64d", "873a88adf5a475c5e65704f16dfbd241ead3283514dc9007d0c9b72c9e411e5a"],
            ["5", "c9079605890523c8941cb5aad076c20c90ec649a94b9537dbe1b8aaec45c61f5", "8c5f8943d22e16eacabf56788185e0978c3a97111a1477d414cf64b51845b0ef"],
            ["0a", "cc61c94724b3428f737d2cd648250b4998f9868a0fcf81392c18dbd19dc21ec8", "878e65aa14c4b54ba082d29a88641a2d3bc5766fc7c041f7f3d4bf877f226189"],
            ["14", "e2d1eeb6fe5bfb4e048099d95dd283ba5916868f0f862bdc8a979129d2337a31", "a11efb1d3a011f306e9e3ac07556a1995d2b13e9000034517d10e3be877bffa2"],
            ["ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc6324d1", "ad311f2c46d5a6173749bba4b3ad9db57ef2b6b9ac62ff5463c5cb817a2ad62a", "fa38d320ec008188f8aa266d75d6b138b46feaf3367834ffb77a8087c2ff3b6d"],
            ["ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc6324d2", "15fe6a86904a36cf6072a061ae619f2870e9016cdddfd92836e84bb6dee35b41", "76759223abe3c14bd0a8879244f403f2fd1c4a970ad602d09ab6968bf6005965"],
            ["ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc6324d3", "6941195b752838c39a7d703660ab52e9519a47b4807a9d289c9635be52bf127a", "91a4ea6d215215f65f153f56aa36de2d8459f5705276171860ffbe2e70da613f"],
            ["ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc6324d4", "a9fe2396bb85b9cb04b76d2d1ed32559f72dab6d225733faaab54cfc93740130", "b16d6af8c3febbc151dc5fac145ff0d52292393b579f3ce2128b0d24bf2219f0"],
            ["ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc6324d5", "d8c0472156db126649b1dfca7f1412750a10ed1576d7996d10f264bc85fc00a2", "8448ecae901289fe9b94c4ad4c99e43958d85051e4dd905a43dfabf3cd768437"],
        ],
        // P-384 test vectors
        [
            ["1", "4d3aadc2299e1513812ff723614ede2b6454868459a30eff879c3afc541b4d6e20e378e2a0d6ce383dd0756649c0b528", "d487543da5ea3a16227ffdd9c69657bf393cade6970b0026745218a9d17c4fae5e40573f844b5653dcfbc253b4fc5b01"],
            ["2", "db93b776427460c39c90a4fd2de4b506da821495f0687f503504e6f0ff9d48a18e6c8f2e022b53f0c8229e55783dde91", "1cb6b808edc20f3df8f2bcf6ff4f197bf60e01beae8d4526ea1b0e7423a77da617171b563d553327bd157b9dcebf4025"],
            ["5", "e1bb26b11d178ffbc676b987e2e8e4c660db2cf5aacc05e67f1eb1b5293c703bc185c0cbcc873466bb595eba68f1f0df", "b24e36291df61b5c193a6bbfdb21dc92b4dc805067d3fbee7f49e4cefe3ca465dd3e91d08d0cba99d496b4608f8c05de"],
            ["0a", "a91a5280db1c212ec2f646025c4845ec54a612d7eb32cabee9f6368eab225fa32a15a9a6e58d6787845539d3c8d99c02", "91466a4137fe8b17c1fcc8f5b85647557e69e473302802e2a6d53f4586b15401fc5acf13ec4617a34468e0881984a031"],
            ["14", "1d1cdf7b5f22346598ca8dd42d96c936f78cff0cd467f03a713466708cbaffc7cd96f20591e71d16ad610a2d94a70ec1", "546227efad1404fea1fe48d4b8e74438bdd259bc0019ecd466cc04da2332b20bda0ab953cde7d7001d85e76c9f2b249f"],
            ["ffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc528f3", "a6661cc4c241720a0336aab8777a16d2f313389118eb5195c0dd449e7c1c39840f4fa5eff21af80ae484fd9f8258030f", "239dcab2a277ac5b2657eedda5791ce3b01f1e0434a02a6b5d85a0e4a022c7a4228968656461e382678db970a7efacd0"],
            ["ffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc528f4", "bbe34e740f395db4a92bb11bd6e0f09d96fde63874231e0feef28f34522c62792aacfa4c569604aff7753246eca101ec", "c7e92a61e3948069f3d6832e367b2f949e29e339ad180b0e88d22a44a467c5195698bdbe018ad5aaaa6e59ddbf943cf9"],
            ["ffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc528f5", "ee82cb5d4c686f3145940572b53625a29c14d45c4d73654245b97c7577b60ca7135cfd8fd9f0f5f0ee101c5ade346cdf", "07b33602c5e1ecacbdf3a4d1afba7f8cf99d6754e04af855f9732cbb4f41b9c31525eb38bfca96d3aa1136b6b0c47624"],
            ["ffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc528f6", "03b2dab245cf48ae7ad6ebad036703ac9436a8cf3450356d50ea65aa46f614e270f25677719d6ada9612b8bc3cfb44bb", "8933bff572f35ccb7cd5a74bb0095802551f41bc219187aaecf53641e119ab8edf3fb730ed4f55a421953d2a4117061e"],
            ["ffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc528f7", "b4bb1133ad9267aa51cb20d25f083167238df4864cb95a00bdb4385fa6f8277692f3dbfda1435507f9e5ed57d99f5989", "49a89c9cf7a12be9b9e99d533cfb098b9544dbb2d4eb6e31045ce25c5e93210c55a4951bc8a66567bf4d13ab4f115bfc"],
        ],
    ];

    // taken from test vectors at ring_ecc/ec/suite_b/ops/
    const ADD_TEST_VECTORS: [[[[&str; 3]; 3]; 5]; 2] = [
        // P-256 test vectors
        [
            [
                ["0000000000000000000000000000000000000000000000000000000000000000", "0000000000000000000000000000000000000000000000000000000000000000", "0000000000000000000000000000000000000000000000000000000000000000"],
                ["18905f76a53755c679fb732b7762251075ba95fc5fedb60179e730d418a9143c", "8571ff1825885d85d2e88688dd21f3258b4ab8e4ba19e45cddf25357ce95560a", "00000000fffffffeffffffffffffffffffffffff000000000000000000000001"],
                ["18905f76a53755c679fb732b7762251075ba95fc5fedb60179e730d418a9143c", "8571ff1825885d85d2e88688dd21f3258b4ab8e4ba19e45cddf25357ce95560a", ""],
            ],
            [
                ["2b11cb945c8cf152ffa4c9c2b1c965b019b35d0b7626919ef0ae6cb9d232f8af", "6d333da42e30f7011245b6281015ded14e0f100968e758a1b6c3c083afc14ea0", "0000000000000000000000000000000000000000000000000000000000000000"],
                ["18905f76a53755c679fb732b7762251075ba95fc5fedb60179e730d418a9143c", "8571ff1825885d85d2e88688dd21f3258b4ab8e4ba19e45cddf25357ce95560a", "00000000fffffffeffffffffffffffffffffffff000000000000000000000001"],
                ["18905f76a53755c679fb732b7762251075ba95fc5fedb60179e730d418a9143c", "8571ff1825885d85d2e88688dd21f3258b4ab8e4ba19e45cddf25357ce95560a", ""],
            ],
            [
                ["18905f76a53755c679fb732b7762251075ba95fc5fedb60179e730d418a9143c", "8571ff1825885d85d2e88688dd21f3258b4ab8e4ba19e45cddf25357ce95560a", "00000000fffffffeffffffffffffffffffffffff000000000000000000000001"],
                ["18905f76a53755c679fb732b7762251075ba95fc5fedb60179e730d418a9143c", "8571ff1825885d85d2e88688dd21f3258b4ab8e4ba19e45cddf25357ce95560a", "00000000fffffffeffffffffffffffffffffffff000000000000000000000001"],
                ["f6bb32e43dcf3a3b732205038d1490d9aa6ae3c1a433827d850046d410ddd64d", "78c577510a5b8a3b19a8fb0e92042dbe152cd7cbeb236ff82f3648d361bee1a5", ""],
            ],
            [
                ["2b11cb945c8cf152ffa4c9c2b1c965b019b35d0b7626919ef0ae6cb9d232f8af", "6d333da42e30f7011245b6281015ded14e0f100968e758a1b6c3c083afc14ea0", "3c396f06c1dc69e4f4b2dce51cd660f761064a4ab098ef61ba3868961f0ef178"],
                ["2b11cb945c8cf152ffa4c9c2b1c965b019b35d0b7626919ef0ae6cb9d232f8af", "6d333da42e30f7011245b6281015ded14e0f100968e758a1b6c3c083afc14ea0", "3c396f06c1dc69e4f4b2dce51cd660f761064a4ab098ef61ba3868961f0ef178"],
                ["f6bb32e43dcf3a3b732205038d1490d9aa6ae3c1a433827d850046d410ddd64d", "873a88adf5a475c5e65704f16dfbd241ead3283514dc9007d0c9b72c9e411e5a", ""],
            ],
            [
                ["18905f76a53755c679fb732b7762251075ba95fc5fedb60179e730d418a9143c", "7a8e00e6da77a27b2d17797722de0cda74b5471c45e61ba3220daca8316aa9f5", "00000000fffffffeffffffffffffffffffffffff000000000000000000000001"],
                ["2b11cb945c8cf152ffa4c9c2b1c965b019b35d0b7626919ef0ae6cb9d232f8af", "6d333da42e30f7011245b6281015ded14e0f100968e758a1b6c3c083afc14ea0", "3c396f06c1dc69e4f4b2dce51cd660f761064a4ab098ef61ba3868961f0ef178"],
                ["f6bb32e43dcf3a3b732205038d1490d9aa6ae3c1a433827d850046d410ddd64d", "873a88adf5a475c5e65704f16dfbd241ead3283514dc9007d0c9b72c9e411e5a", ""],
            ]
        ],
        [
            [
                ["000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"],
                ["4d3aadc2299e1513812ff723614ede2b6454868459a30eff879c3afc541b4d6e20e378e2a0d6ce383dd0756649c0b528", "2b78abc25a15c5e9dd8002263969a840c6c3521968f4ffd98bade7562e83b050a1bfa8bf7bb4a9ac23043dad4b03a4fe", "000000000000000000000000000000000000000000000000000000000000000100000000ffffffffffffffff00000001"],
                ["4d3aadc2299e1513812ff723614ede2b6454868459a30eff879c3afc541b4d6e20e378e2a0d6ce383dd0756649c0b528", "2b78abc25a15c5e9dd8002263969a840c6c3521968f4ffd98bade7562e83b050a1bfa8bf7bb4a9ac23043dad4b03a4fe", ""],
            ],
            [
                ["4a0fd63f894499928e4b2b72aced45cfc589976f4ff86f78c904d59da9379a62b702d968c1184834c11db28c7356ceb6", "be113b04484cd4bc215a9f2a33a674c3764c38ca4de135dd50ce8dcf3c85d55a5aad0e171860bdb6c58201e6212d9ac5", "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"],
                ["4d3aadc2299e1513812ff723614ede2b6454868459a30eff879c3afc541b4d6e20e378e2a0d6ce383dd0756649c0b528", "2b78abc25a15c5e9dd8002263969a840c6c3521968f4ffd98bade7562e83b050a1bfa8bf7bb4a9ac23043dad4b03a4fe", "000000000000000000000000000000000000000000000000000000000000000100000000ffffffffffffffff00000001"],
                ["4d3aadc2299e1513812ff723614ede2b6454868459a30eff879c3afc541b4d6e20e378e2a0d6ce383dd0756649c0b528", "2b78abc25a15c5e9dd8002263969a840c6c3521968f4ffd98bade7562e83b050a1bfa8bf7bb4a9ac23043dad4b03a4fe", ""],
            ],
            [
                ["4d3aadc2299e1513812ff723614ede2b6454868459a30eff879c3afc541b4d6e20e378e2a0d6ce383dd0756649c0b528", "2b78abc25a15c5e9dd8002263969a840c6c3521968f4ffd98bade7562e83b050a1bfa8bf7bb4a9ac23043dad4b03a4fe", "000000000000000000000000000000000000000000000000000000000000000100000000ffffffffffffffff00000001"],
                ["4d3aadc2299e1513812ff723614ede2b6454868459a30eff879c3afc541b4d6e20e378e2a0d6ce383dd0756649c0b528", "2b78abc25a15c5e9dd8002263969a840c6c3521968f4ffd98bade7562e83b050a1bfa8bf7bb4a9ac23043dad4b03a4fe", "000000000000000000000000000000000000000000000000000000000000000100000000ffffffffffffffff00000001"],
                ["db93b776427460c39c90a4fd2de4b506da821495f0687f503504e6f0ff9d48a18e6c8f2e022b53f0c8229e55783dde91", "e34947f7123df0c2070d430900b0e68409f1fe415172bad915e4f18bdc588258e8e8e4a8c2aaccd842ea84633140bfda", ""],
            ],
            [
                ["f3ee335326d22614d01b5d7cd0be73f1bfdd75982c9c273f72d0abfeecbca0431601a1bcafcdeb07e21ecf4d91c7b520", "57b82ca1527c5a01b78bc8ccb9febe74178b04b7c6fde1c1c4ef9a220c4320bb560cb078542256a3900df61c107de6c5", "53b3adc887551c0e17c07ecb42d1a5ec105aeec6b0f040a936ed4f756e83939226232b4e11191b3eb1d841c650682ca0"],
                ["4d3aadc2299e1513812ff723614ede2b6454868459a30eff879c3afc541b4d6e20e378e2a0d6ce383dd0756649c0b528", "d487543da5ea3a16227ffdd9c69657bf393cade6970b0026745218a9d17c4fae5e40573f844b5653dcfbc253b4fc5b01", "000000000000000000000000000000000000000000000000000000000000000100000000ffffffffffffffff00000001"],
                ["db93b776427460c39c90a4fd2de4b506da821495f0687f503504e6f0ff9d48a18e6c8f2e022b53f0c8229e55783dde91", "1cb6b808edc20f3df8f2bcf6ff4f197bf60e01beae8d4526ea1b0e7423a77da617171b563d553327bd157b9dcebf4025", ""],
            ],
            [
                ["f3ee335326d22614d01b5d7cd0be73f1bfdd75982c9c273f72d0abfeecbca0431601a1bcafcdeb07e21ecf4d91c7b520", "57b82ca1527c5a01b78bc8ccb9febe74178b04b7c6fde1c1c4ef9a220c4320bb560cb078542256a3900df61c107de6c5", "53b3adc887551c0e17c07ecb42d1a5ec105aeec6b0f040a936ed4f756e83939226232b4e11191b3eb1d841c650682ca0"],
                ["f3ee335326d22614d01b5d7cd0be73f1bfdd75982c9c273f72d0abfeecbca0431601a1bcafcdeb07e21ecf4d91c7b520", "57b82ca1527c5a01b78bc8ccb9febe74178b04b7c6fde1c1c4ef9a220c4320bb560cb078542256a3900df61c107de6c5", "53b3adc887551c0e17c07ecb42d1a5ec105aeec6b0f040a936ed4f756e83939226232b4e11191b3eb1d841c650682ca0"],
                ["db93b776427460c39c90a4fd2de4b506da821495f0687f503504e6f0ff9d48a18e6c8f2e022b53f0c8229e55783dde91", "1cb6b808edc20f3df8f2bcf6ff4f197bf60e01beae8d4526ea1b0e7423a77da617171b563d553327bd157b9dcebf4025", ""],
            ]
        ]
    ];
}