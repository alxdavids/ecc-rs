use num::{BigUint,Zero,One};
use untrusted;
use core::marker::PhantomData;

use crate::ring_ecc;
use ring_ecc::ec::CurveID;
use ring_ecc::ec::suite_b::ops::Point as RingPoint;
use ring_ecc::ec::suite_b::ops::PrivateKeyOps as CurveAuxOps;
use ring_ecc::ec::suite_b::ops::{Elem,CommonOps,Scalar};
use ring_ecc::ec::suite_b::ops::p256::COMMON_OPS as P256_COMMON_OPS;
use ring_ecc::ec::suite_b::ops::p256::PRIVATE_KEY_OPS as P256_AUX_OPS;
use ring_ecc::ec::suite_b::ops::p384::COMMON_OPS as P384_COMMON_OPS;
use ring_ecc::ec::suite_b::ops::p384::PRIVATE_KEY_OPS as P384_AUX_OPS;
use ring_ecc::limb::{Limb,LIMB_BYTES};
use ring_ecc::limb::big_endian_from_limbs;
use ring_ecc::ec::suite_b::ops::elem::MAX_LIMBS;
use ring_ecc::ec::suite_b::ops::{elem_parse_big_endian_fixed_consttime,scalar_parse_big_endian_fixed_consttime};
use ring_ecc::arithmetic::montgomery::R;

/// AffinePoint
struct AffinePoint {
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
            CurveID::P256 => Self {
                x: Zero::zero(),
                y: Zero::zero(),
                id: CurveID::P256,
                curve_params: &P256_COMMON_OPS,
                aux: &P256_AUX_OPS,
            },
            CurveID::P384 => Self {
                x: Zero::zero(),
                y: Zero::zero(),
                id: CurveID::P384,
                curve_params: &P384_COMMON_OPS,
                aux: &P384_AUX_OPS,
            },
            _ => panic!("unsupported"),
        }
    }

    pub fn get_generator(id: CurveID) -> Self {
        match id {
            CurveID::P256 => Self {
                // ring doesn't explicitly make the P256 generator available,
                // but some applications need it. we use the one in
                // ring_ecc/ec/suite_b/ops/p256_point_mul_tests.txt because it
                // differs from the standard generator
                x: BigUint::parse_bytes(b"18905f76a53755c679fb732b7762251075ba95fc5fedb60179e730d418a9143c", 16).unwrap(),
                y: BigUint::parse_bytes(b"7a8e00e6da77a27b2d17797722de0cda74b5471c45e61ba3220daca8316aa9f5", 16).unwrap(),
                id: CurveID::P256,
                curve_params: &P256_COMMON_OPS,
                aux: &P256_AUX_OPS,
            },
            CurveID::P384 => {
                // for some reason the ring generator does not agree with
                // https://www.secg.org/SEC2-Ver-1.0.pdf Section 2.8.1. using
                // the point detailed in
                // ring_ecc/ec/suite_b/ops/p384_point_mul_tests.txt for now
                Self {
                    x: BigUint::parse_bytes(b"4d3aadc2299e1513812ff723614ede2b6454868459a30eff879c3afc541b4d6e20e378e2a0d6ce383dd0756649c0b528", 16).unwrap(),
                    y: BigUint::parse_bytes(b"d487543da5ea3a16227ffdd9c69657bf393cade6970b0026745218a9d17c4fae5e40573f844b5653dcfbc253b4fc5b01", 16).unwrap(),
                    id: CurveID::P384,
                    curve_params: &P384_COMMON_OPS,
                    aux: &P384_AUX_OPS,
                }
            },
            _ => panic!("unsupported"),
        }
    }

    pub fn base_mul(id: CurveID, scalar: &Scalar) -> JacobianPoint {
        AffinePoint::get_generator(id).scalar_mul(scalar)
    }

    pub fn scalar_mul(&self, scalar: &Scalar) -> JacobianPoint {
        let x_elem = big_uint_to_elem(self.curve_params, &self.x);
        let y_elem = big_uint_to_elem(self.curve_params, &self.y);
        let ring_pt = self.aux.point_mul(scalar, &(x_elem, y_elem));
        JacobianPoint::from_ring_jac_point(ring_pt, self.id, self.curve_params, self.aux)
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
}

struct JacobianPoint {
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
            CurveID::P256 => AffinePoint {
                x: Zero::zero(),
                y: Zero::zero(),
                id: CurveID::P256,
                curve_params: &P256_COMMON_OPS,
                aux: &P256_AUX_OPS,
            },
            CurveID::P384 => AffinePoint {
                x: Zero::zero(),
                y: Zero::zero(),
                id: CurveID::P384,
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

    pub fn to_affine(&self) -> AffinePoint {
        let ops = self.curve_params;
        let z = big_uint_to_elem(ops, &self.z);

        // compute z^{-2}
        let z_inv_sq = self.aux.elem_inverse_squared(&z);

        // compute z^{-3}
        let mut z_inv_cub = z_inv_sq.clone();
        ops.elem_mul(&mut z_inv_cub, &z_inv_sq);
        ops.elem_mul(&mut z_inv_cub, &z);

        // compute x*z^{-2}
        let mut x_aff = big_uint_to_elem(ops, &self.x);
        ops.elem_mul(&mut x_aff, &z_inv_sq);

        // compute y*z^{-3}
        let mut y_aff = big_uint_to_elem(ops, &self.y);
        ops.elem_mul(&mut y_aff, &z_inv_cub);

        // output affine
        AffinePoint {
            x: elem_to_big_uint(ops, x_aff),
            y: elem_to_big_uint(ops, y_aff),
            id: self.id,
            curve_params: ops,
            aux: self.aux
        }
    }

    fn as_ring_jac_point(&self) -> RingPoint {
        let (x, y, z) = self.encode_for_ring_ops();
        let mut limbs: [Limb; MAX_LIMBS*3] = [0; MAX_LIMBS*3];
        limbs[..self.curve_params.num_limbs].copy_from_slice(&x.limbs);
        limbs[MAX_LIMBS..MAX_LIMBS+self.curve_params.num_limbs].copy_from_slice(&y.limbs);
        limbs[MAX_LIMBS*2..MAX_LIMBS*2+self.curve_params.num_limbs].copy_from_slice(&z.limbs);
        RingPoint { xyz: limbs }
    }

    fn from_ring_jac_point(ring_pt: RingPoint, id: CurveID, ops: &'static CommonOps, aux: &'static CurveAuxOps) -> Self {
        let x = elem_to_big_uint(ops, ops.point_x(&ring_pt));
        let y = elem_to_big_uint(ops, ops.point_y(&ring_pt));
        let z = elem_to_big_uint(ops, ops.point_z(&ring_pt));
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
            big_uint_to_elem(self.curve_params, &self.x),
            big_uint_to_elem(self.curve_params, &self.y),
            big_uint_to_elem(self.curve_params, &self.z),
        )
    }
}

fn biguint_to_scalar(curve_params: &CommonOps, x: &BigUint) -> Scalar {
    let x_bytes = x.to_bytes_be();
    let MAX_LEN = curve_params.num_limbs*LIMB_BYTES;
    let mut inp = vec![0; MAX_LEN];
    inp[MAX_LEN-x_bytes.len()..].copy_from_slice(&x_bytes);
    scalar_parse_big_endian_fixed_consttime(curve_params, untrusted::Input::from(&inp)).unwrap()
}

fn big_uint_to_elem(curve_params: &CommonOps, x: &BigUint) -> Elem<R> {
    let x_bytes = x.to_bytes_be();
    let MAX_LEN = curve_params.num_limbs*LIMB_BYTES;
    let mut inp = vec![0; MAX_LEN];
    inp[MAX_LEN-x_bytes.len()..].copy_from_slice(&x_bytes);
    bytes_to_elem(curve_params, &inp)
}

fn elem_to_big_uint(ops: &'static CommonOps, elem: Elem<R>) -> BigUint {
    BigUint::from_bytes_be(&elem_to_bytes(ops, elem))
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
fn elem_to_bytes(ops: &CommonOps, elem: Elem<R>) -> Vec<u8> {
    println!("{}", ops.num_limbs);
    let mut out: Vec<u8> = vec![0; MAX_LIMBS*LIMB_BYTES];
    big_endian_from_limbs(&elem.limbs, &mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use hex;

    const MULT_TEST_VECTORS: [[[&str; 3]; 10]; 2] = [
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

    #[test]
    fn scalar_mul_test() {
        for i in 0..2 {
            let gen = match i {
                0 => AffinePoint::get_generator(CurveID::P256),
                1 => AffinePoint::get_generator(CurveID::P384),
                _ => panic!("test failed")
            };
            for vector in MULT_TEST_VECTORS[i].iter() {
                let k = biguint_to_scalar(gen.curve_params, &BigUint::parse_bytes(vector[0].as_bytes(), 16).unwrap());
                let aff = gen.scalar_mul(&k).to_affine();

                // hex encode padding 0's (to_str_radix doesn't do this)
                let x_hex = hex::encode(aff.x.to_bytes_be());
                let y_hex = hex::encode(aff.y.to_bytes_be());

                assert_eq!(x_hex, vector[1], "x hex coordinate check for curve: {} and scalar: {}", i, vector[0]);
                assert_eq!(y_hex, vector[2], "y hex coordinate check for curve: {} and scalar: {}", i, vector[0]);
            }
        }
    }

    fn add_test() {

    }
}