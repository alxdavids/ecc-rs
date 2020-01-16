/// The point module provides the main interface for performing elliptic curve
/// operations. Currently we only support the `secp256r1` and `secp384r1` curves,
/// both of which are implemented in the `suite_b` module of ring.
///
/// # Examples
///
/// Performing operations on elliptic curve points is as easy as running:
///
/// ```
/// use ecc_rs::point::*;
/// use num::{BigUint,One};
///
/// // get new zero point
/// let _ = AffinePoint::new(P256);
///
/// // get generator point
/// let g = AffinePoint::get_generator(P256);
///
/// // perform scalar mult
/// let k = BigUint::parse_bytes(b"0a", 16).unwrap();
/// let k_g = g.scalar_mul(&k);
///
/// // add points
/// let k_1_g = g.to_jacobian().add(&k_g).to_affine();
///
/// // check point is valid
/// assert_eq!(k_1_g.is_valid(), true);
///
/// // check points are equal
/// let k_1 = BigUint::parse_bytes(b"0b", 16).unwrap();
/// let k_1_chk_g = g.scalar_mul(&k_1).to_affine();
/// assert_eq!(k_1_g.equals(k_1_chk_g), true);
/// ```
///
/// We can also create operate over Jacobian points:
///
/// ```
/// use ecc_rs::point::*;
/// use num::{BigUint,One};
///
/// // get generator point
/// let g = AffinePoint::get_generator(P256);
/// let g_jac = AffinePoint::get_generator(P256).to_jacobian();
///
/// // add points
/// let two_g = g_jac.add(&g_jac);
///
/// // check point is valid
/// assert_eq!(two_g.is_valid(), true);
///
/// // check points are equal
/// let two_sc = BigUint::parse_bytes(b"2", 16).unwrap();
/// let two_sc_g = g.scalar_mul(&two_sc);
/// assert_eq!(two_g.equals(two_sc_g), true);
/// ```

use num::{BigUint,Zero,One};
use untrusted;
use core::marker::PhantomData;

use crate::ring_ecc;
use ring_ecc::ec::CurveID;
use ring_ecc::ec::suite_b::verify_affine_point_is_on_the_curve;
use ring_ecc::ec::suite_b::ops::Point as RingPoint;
use ring_ecc::ec::suite_b::ops::PrivateKeyOps as CurveAuxOps;
use ring_ecc::ec::suite_b::ops::{Elem,CommonOps,Scalar};
use ring_ecc::ec::suite_b::ops::p256::COMMON_OPS as P256_COMMON_OPS;
use ring_ecc::ec::suite_b::ops::p256::PRIVATE_KEY_OPS as P256_AUX_OPS;
use ring_ecc::ec::suite_b::ops::p384::COMMON_OPS as P384_COMMON_OPS;
use ring_ecc::ec::suite_b::ops::p384::PRIVATE_KEY_OPS as P384_AUX_OPS;
use ring_ecc::ec::suite_b::ops::p384::P384_GENERATOR;
use ring_ecc::limb::{Limb,LIMB_BYTES};
use ring_ecc::limb::big_endian_from_limbs;
use ring_ecc::ec::suite_b::ops::elem::MAX_LIMBS;
use ring_ecc::ec::suite_b::ops::{elem_parse_big_endian_fixed_consttime,scalar_parse_big_endian_fixed_consttime};
use ring_ecc::arithmetic::montgomery::R;

/// P256 constant for exposing ring CurveID::P256 enum
pub const P256: CurveID = CurveID::P256;
/// P384 constant for exposing ring CurveID::P384 enum
pub const P384: CurveID = CurveID::P384;

/// The `AffinePoint` struct provides access to a base elliptic curve point
/// representation. The setting of `id`, `curve_params` and `aux` determine
/// which curve the point is associated with. Currently we support `secp256r1`
/// and `secp384r1` as this is what is supported in ring.
pub struct AffinePoint {
    x: BigUint,
    y: BigUint,
    id: CurveID,
    curve_params: &'static CommonOps,
    aux: &'static CurveAuxOps,
}

impl AffinePoint {
    /// Returns the identity point
    pub fn new(id: CurveID) -> Self {
        match id {
            P256 => Self {
                x: Zero::zero(),
                y: Zero::zero(),
                id: P256,
                curve_params: &P256_COMMON_OPS,
                aux: &P256_AUX_OPS,
            },
            P384 => Self {
                x: Zero::zero(),
                y: Zero::zero(),
                id: P384,
                curve_params: &P384_COMMON_OPS,
                aux: &P384_AUX_OPS,
            },
            _ => panic!("unsupported"),
        }
    }

    pub fn get_generator(id: CurveID) -> Self {
        match id {
            P256 => Self {
                // ring doesn't explicitly make the P256 generator available,
                // but some applications need it. we use the one in
                // ring_ecc/ec/suite_b/ops/p256_point_mul_tests.txt because it
                // differs from the standard generator
                x: BigUint::parse_bytes(b"18905f76a53755c679fb732b7762251075ba95fc5fedb60179e730d418a9143c", 16).unwrap(),
                y: BigUint::parse_bytes(b"8571ff1825885d85d2e88688dd21f3258b4ab8e4ba19e45cddf25357ce95560a", 16).unwrap(),
                id: P256,
                curve_params: &P256_COMMON_OPS,
                aux: &P256_AUX_OPS,
            },
            P384 => {
                // for some reason the ring generator does not agree with
                // https://www.secg.org/SEC2-Ver-1.0.pdf Section 2.8.1. using
                // the ring point for now
                Self {
                    x: elem_to_biguint(P384_GENERATOR.0),
                    y: elem_to_biguint(P384_GENERATOR.1),
                    id: P384,
                    curve_params: &P384_COMMON_OPS,
                    aux: &P384_AUX_OPS,
                }
            },
            _ => panic!("unsupported"),
        }
    }

    pub fn base_mul(id: CurveID, scalar: &BigUint) -> JacobianPoint {
        AffinePoint::get_generator(id).scalar_mul(scalar)
    }

    pub fn scalar_mul(&self, scalar: &BigUint) -> JacobianPoint {
        let ring_sc = biguint_to_scalar(self.curve_params, scalar);
        let x_elem = biguint_to_elem(self.curve_params, &self.x);
        let y_elem = biguint_to_elem(self.curve_params, &self.y);
        let ring_pt = self.aux.point_mul(&ring_sc, &(x_elem, y_elem));
        JacobianPoint::from_ring_jac_point(ring_pt, self.id, self.curve_params, self.aux)
    }

    pub fn equals(&self, other: Self) -> bool {
        (self.x == other.x)
        && (self.y == other.y)
        && (self.id == other.id)
    }

    pub fn is_valid(&self) -> bool {
        let (x, y) = self.as_ring_affine_point();
        match verify_affine_point_is_on_the_curve(self.curve_params, (&x, &y)) {
            Ok(_) => true,
            Err(_) => false,
        }
    }

    /// returns the point with y = -y
    pub fn get_minus_y_point(&self) -> AffinePoint {
        AffinePoint {
            x: self.x.clone(),
            y: self.get_curve_modulus_as_biguint() - &self.y,
            id: self.id,
            curve_params: self.curve_params,
            aux: self.aux
        }
    }

    /// to jacobian
    pub fn to_jacobian(&self) -> JacobianPoint {
        JacobianPoint {
            x: self.x.clone(),
            y: self.y.clone(),
            z: One::one(),
            id: self.id,
            curve_params: self.curve_params,
            aux: self.aux
        }
    }

    fn as_ring_affine_point(&self) -> (Elem<R>, Elem<R>) {
        (biguint_to_elem(self.curve_params, &self.x), biguint_to_elem(self.curve_params, &self.y))
    }

    fn get_curve_modulus_as_biguint(&self) -> BigUint {
        let p = Elem {
            limbs: self.curve_params.q.p,
            m: PhantomData,
            encoding: PhantomData
        };
        elem_to_biguint( p)
    }
}

pub struct JacobianPoint {
    x: BigUint,
    y: BigUint,
    z: BigUint,
    id: CurveID,
    curve_params: &'static CommonOps,
    aux: &'static CurveAuxOps,
}

impl JacobianPoint {
    /// Returns the identity point
    pub fn new(id: CurveID) -> Self {
        match id {
            P256 => AffinePoint {
                x: Zero::zero(),
                y: Zero::zero(),
                id: P256,
                curve_params: &P256_COMMON_OPS,
                aux: &P256_AUX_OPS,
            },
            P384 => AffinePoint {
                x: Zero::zero(),
                y: Zero::zero(),
                id: P384,
                curve_params: &P384_COMMON_OPS,
                aux: &P384_AUX_OPS,
            },
            _ => panic!("unsupported"),
        }.to_jacobian()
    }

    pub fn add(&self, other: &JacobianPoint) -> JacobianPoint {
        let p1 = self.as_ring_jac_point();
        let p2 = other.as_ring_jac_point();
        let p_add = self.curve_params.point_sum(&p1, &p2);
        Self::from_ring_jac_point(p_add, self.id, self.curve_params, self.aux)
    }

    pub fn equals(&self, other: Self) -> bool {
        (self.x == other.x)
        && (self.y == other.y)
        && (self.z == other.z)
        && (self.id == other.id)
    }

    pub fn is_valid(&self) -> bool {
        self.to_affine().is_valid()
    }

    pub fn to_affine(&self) -> AffinePoint {
        let ops = self.curve_params;
        let z = biguint_to_elem(ops, &self.z);

        // compute z^{-2}
        let z_inv_sq = self.aux.elem_inverse_squared(&z);

        // compute z^{-3}
        let mut z_inv_cub = z_inv_sq.clone();
        ops.elem_mul(&mut z_inv_cub, &z_inv_sq);
        ops.elem_mul(&mut z_inv_cub, &z);

        // compute x*z^{-2}
        let mut x_aff = biguint_to_elem(ops, &self.x);
        ops.elem_mul(&mut x_aff, &z_inv_sq);

        // compute y*z^{-3}
        let mut y_aff = biguint_to_elem(ops, &self.y);
        ops.elem_mul(&mut y_aff, &z_inv_cub);

        // output affine
        AffinePoint {
            x: elem_to_biguint(x_aff),
            y: elem_to_biguint(y_aff),
            id: self.id,
            curve_params: ops,
            aux: self.aux
        }
    }

    fn as_ring_jac_point(&self) -> RingPoint {
        let (x, y, z) = self.encode_for_ring_ops();
        let mut limbs: [Limb; MAX_LIMBS*3] = [0; MAX_LIMBS*3];
        limbs[..MAX_LIMBS].copy_from_slice(&x.limbs);
        limbs[self.curve_params.num_limbs..self.curve_params.num_limbs+MAX_LIMBS].copy_from_slice(&y.limbs);
        limbs[self.curve_params.num_limbs*2..self.curve_params.num_limbs*2+MAX_LIMBS].copy_from_slice(&z.limbs);
        RingPoint { xyz: limbs }
    }

    fn from_ring_jac_point(ring_pt: RingPoint, id: CurveID, ops: &'static CommonOps, aux: &'static CurveAuxOps) -> Self {
        let x = elem_to_biguint(ops.point_x(&ring_pt));
        let y = elem_to_biguint(ops.point_y(&ring_pt));
        let z = elem_to_biguint(ops.point_z(&ring_pt));
        Self {
            x: x,
            y: y,
            z: z,
            id: id,
            curve_params: ops,
            aux: aux
        }
    }

    /// encodes the fields of the point
    fn encode_for_ring_ops(&self) -> (Elem<R>, Elem<R>, Elem<R>) {
        (
            biguint_to_elem(self.curve_params, &self.x),
            biguint_to_elem(self.curve_params, &self.y),
            biguint_to_elem(self.curve_params, &self.z),
        )
    }
}

#[allow(dead_code)]
fn biguint_to_scalar(curve_params: &CommonOps, x: &BigUint) -> Scalar {
    let x_bytes = x.to_bytes_be();
    let MAX_LEN = curve_params.num_limbs*LIMB_BYTES;
    let mut inp = vec![0; MAX_LEN];
    inp[MAX_LEN-x_bytes.len()..].copy_from_slice(&x_bytes);
    scalar_parse_big_endian_fixed_consttime(curve_params, untrusted::Input::from(&inp)).unwrap()
}

fn biguint_to_elem(curve_params: &CommonOps, x: &BigUint) -> Elem<R> {
    let x_bytes = x.to_bytes_be();
    let MAX_LEN = curve_params.num_limbs*LIMB_BYTES;
    let mut inp = vec![0; MAX_LEN];
    inp[MAX_LEN-x_bytes.len()..].copy_from_slice(&x_bytes);
    bytes_to_elem(curve_params, &inp)
}

fn elem_to_biguint(elem: Elem<R>) -> BigUint {
    BigUint::from_bytes_be(&elem_to_bytes(elem))
}

/// parses an elem from be_bytes
fn bytes_to_elem(curve_params: &CommonOps, bytes: &[u8]) -> Elem<R> {
    Elem {
        limbs: elem_parse_big_endian_fixed_consttime(curve_params, untrusted::Input::from(bytes)).unwrap().limbs,
        encoding: PhantomData,
        m: PhantomData
    }
}

/// transform elem to be_bytes
fn elem_to_bytes(elem: Elem<R>) -> Vec<u8> {
    let mut out: Vec<u8> = vec![0; MAX_LIMBS*LIMB_BYTES];
    big_endian_from_limbs(&elem.limbs, &mut out);
    out
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

            // for some reason the ring test vectors use their bespoke generator
            // with the minus y coordinate
            let test_gen = gen.get_minus_y_point();
            for vector in MULT_TEST_VECTORS[i].iter() {
                let k = BigUint::parse_bytes(vector[0].as_bytes(), 16).unwrap();
                let aff = test_gen.scalar_mul(&k).to_affine();

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
            let (ops, aux, id) = match i {
                0 => (&P256_COMMON_OPS, &P256_AUX_OPS, P256),
                1 => (&P384_COMMON_OPS, &P384_AUX_OPS, P384),
                _ => panic!("test failed")
            };
            let vectors = ADD_TEST_VECTORS[i];
            for j in 0..vectors.len() {
                let pt_left = JacobianPoint {
                    x: BigUint::parse_bytes(vectors[j][0][0].as_bytes(), 16).unwrap(),
                    y: BigUint::parse_bytes(vectors[j][0][1].as_bytes(), 16).unwrap(),
                    z: BigUint::parse_bytes(vectors[j][0][2].as_bytes(), 16).unwrap(),
                    id: id,
                    curve_params: ops,
                    aux: aux
                };
                let pt_right = JacobianPoint {
                    x: BigUint::parse_bytes(vectors[j][1][0].as_bytes(), 16).unwrap(),
                    y: BigUint::parse_bytes(vectors[j][1][1].as_bytes(), 16).unwrap(),
                    z: BigUint::parse_bytes(vectors[j][1][2].as_bytes(), 16).unwrap(),
                    id: id,
                    curve_params: ops,
                    aux: aux
                };

                let aff = pt_left.add(&pt_right).to_affine();
                let x_hex = hex::encode(aff.x.to_bytes_be());
                let y_hex = hex::encode(aff.y.to_bytes_be());

                assert_eq!(aff.is_valid(), true);
                assert_eq!(x_hex, vectors[j][2][0], "x hex coordinate check for addition curve: {:?} and j: {}", id, j);
                assert_eq!(y_hex, vectors[j][2][1], "x hex coordinate check for addition curve: {:?} and j: {}", id, j);
            }
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