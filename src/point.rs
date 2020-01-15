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
use ring_ecc::ec::suite_b::ops::p384::P384_GENERATOR;
use ring_ecc::limb::{Limb,LIMB_BYTES};
use ring_ecc::limb::big_endian_from_limbs;
use ring_ecc::ec::suite_b::ops::elem::MAX_LIMBS;
use ring_ecc::ec::suite_b::ops::{elem_parse_big_endian_fixed_consttime,scalar_parse_big_endian_fixed_consttime};
use ring_ecc::arithmetic::montgomery::R;

/// AffinePoint
struct AffinePoint {
    x: BigUint,
    y: BigUint,
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
                curve_params: &P256_COMMON_OPS,
                aux: &P256_AUX_OPS,
            },
            CurveID::P384 => Self {
                x: Zero::zero(),
                y: Zero::zero(),
                curve_params: &P384_COMMON_OPS,
                aux: &P384_AUX_OPS,
            },
            _ => panic!("unsupported"),
        }
    }

    pub fn scalar_mul(&self, scalar: &Scalar) -> JacobianPoint {
        let x_elem = big_uint_to_elem(self.curve_params, &self.x);
        let y_elem = big_uint_to_elem(self.curve_params, &self.y);
        let ring_pt = self.aux.point_mul(scalar, &(x_elem, y_elem));
        println!("{:?}", ring_pt.xyz);
        JacobianPoint::from_ring_jac_point(ring_pt, self.curve_params, self.aux)
    }

    /// to jacobian
    pub fn to_jacobian(&self) -> JacobianPoint {
        JacobianPoint {
            x: self.x.clone(),
            y: self.y.clone(),
            z: One::one(),
            curve_params: self.curve_params,
            aux: self.aux
        }
    }
}

struct JacobianPoint {
    x: BigUint,
    y: BigUint,
    z: BigUint,
    curve_params: &'static CommonOps,
    aux: &'static CurveAuxOps,
}

impl JacobianPoint {
    /// Returns the identity point
    pub fn new(id: CurveID) -> Self {
        match id {
            CurveID::P256 => {
                AffinePoint {
                    x: Zero::zero(),
                    y: Zero::zero(),
                    curve_params: &P256_COMMON_OPS,
                    aux: &P256_AUX_OPS,
                }
            },
            CurveID::P384 => AffinePoint {
                x: Zero::zero(),
                y: Zero::zero(),
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
        Self::from_ring_jac_point(p_add, self.curve_params, self.aux)
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

    fn from_ring_jac_point(ring_pt: RingPoint, ops: &'static CommonOps, aux: &'static CurveAuxOps) -> Self {
        let x = elem_to_big_uint(ops, ops.point_x(&ring_pt));
        let y = elem_to_big_uint(ops, ops.point_y(&ring_pt));
        let z = elem_to_big_uint(ops, ops.point_z(&ring_pt));
        Self {
            x: x,
            y: y,
            z: z,
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
    println!("be: {:?}", x_bytes);
    let MAX_LEN = curve_params.num_limbs*LIMB_BYTES;
    let mut inp = vec![0; MAX_LEN];
    inp[MAX_LEN-x_bytes.len()..].copy_from_slice(&x_bytes);
    scalar_parse_big_endian_fixed_consttime(curve_params, untrusted::Input::from(&inp)).unwrap()
}

fn big_uint_to_elem(curve_params: &CommonOps, x: &BigUint) -> Elem<R> {
    let x_bytes = x.to_bytes_be();
    let MAX_LEN = curve_params.num_limbs*LIMB_BYTES;
    let mut inp = vec![0; MAX_LEN];
    println!("x: {:?}", x_bytes.len());
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
    let mut out: Vec<u8> = vec![0; ops.num_limbs*LIMB_BYTES];
    big_endian_from_limbs(&elem.limbs, &mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_VECTORS: [[[&str; 3]; 4]; 2] = [
        [
            ["1", "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296", "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5"],
            ["2", "7cf27b188d034f7e8a52380304b51ac3c08969e277f21b35a60b48fc47669978", "07775510db8ed040293d9ac69f7430dbba7dade63ce982299e04b79d227873d1"],
            ["5", "51590b7a515140d2d784c85608668fdfef8c82fd1f5be52421554a0dc3d033ed", "e0c17da8904a727d8ae1bf36bf8a79260d012f00d4d80888d1d0bb44fda16da4"],
            ["10", "cef66d6b2a3a993e591214d1ea223fb545ca6c471c48306e4c36069404c5723f", "878662a229aaae906e123cdd9d3b4c10590ded29fe751eeeca34bbaa44af0773"],
        ],
        [
            ["1", "4d3aadc2299e1513812ff723614ede2b6454868459a30eff879c3afc541b4d6e20e378e2a0d6ce383dd0756649c0b528", "d487543da5ea3a16227ffdd9c69657bf393cade6970b0026745218a9d17c4fae5e40573f844b5653dcfbc253b4fc5b01"],
            ["2", "db93b776427460c39c90a4fd2de4b506da821495f0687f503504e6f0ff9d48a18e6c8f2e022b53f0c8229e55783dde91", "1cb6b808edc20f3df8f2bcf6ff4f197bf60e01beae8d4526ea1b0e7423a77da617171b563d553327bd157b9dcebf4025"],
            ["5", "e1bb26b11d178ffbc676b987e2e8e4c660db2cf5aacc05e67f1eb1b5293c703bc185c0cbcc873466bb595eba68f1f0df", "b24e36291df61b5c193a6bbfdb21dc92b4dc805067d3fbee7f49e4cefe3ca465dd3e91d08d0cba99d496b4608f8c05de"],
            ["10", "a91a5280db1c212ec2f646025c4845ec54a612d7eb32cabee9f6368eab225fa32a15a9a6e58d6787845539d3c8d99c02", "91466a4137fe8b17c1fcc8f5b85647557e69e473302802e2a6d53f4586b15401fc5acf13ec4617a34468e0881984a031"],
        ],
    ];

    #[test]
    fn scalar_mul_test() {
        for i in 1..2 {
            let ops = match i {
                0 => &P256_COMMON_OPS,
                1 => &P384_COMMON_OPS,
                _ => panic!("blah")
            };
            let gen = AffinePoint {
                x: BigUint::parse_bytes(b"4d3aadc2299e1513812ff723614ede2b6454868459a30eff879c3afc541b4d6e20e378e2a0d6ce383dd0756649c0b528", 16).unwrap(),
                y: BigUint::parse_bytes(b"d487543da5ea3a16227ffdd9c69657bf393cade6970b0026745218a9d17c4fae5e40573f844b5653dcfbc253b4fc5b01", 16).unwrap(),
                curve_params: ops,
                aux: &P384_AUX_OPS,
            };
            for vector in TEST_VECTORS[i].iter() {
                let k = biguint_to_scalar(gen.curve_params, &BigUint::parse_bytes(vector[0].as_bytes(), 10).unwrap());
                let aff = gen.scalar_mul(&k).to_affine();
                let x_hex = aff.x.to_str_radix(16);
                let y_hex = aff.y.to_str_radix(16);

                assert_eq!(x_hex, vector[1], "x hex coordinate check for curve: {} and scalar: {}", i, vector[0]);
                assert_eq!(y_hex, vector[2], "y hex coordinate check for curve: {} and scalar: {}", i, vector[0]);
            }
        }
    }
}