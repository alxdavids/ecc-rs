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

    const TEST_VECTORS: [[[&str; 3]; 10]; 2] = [
        [
            ["1", "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296", "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5"],
            ["2", "7cf27b188d034f7e8a52380304b51ac3c08969e277f21b35a60b48fc47669978", "07775510db8ed040293d9ac69f7430dbba7dade63ce982299e04b79d227873d1"],
            ["5", "51590b7a515140d2d784c85608668fdfef8c82fd1f5be52421554a0dc3d033ed", "e0c17da8904a727d8ae1bf36bf8a79260d012f00d4d80888d1d0bb44fda16da4"],
            ["10", "cef66d6b2a3a993e591214d1ea223fb545ca6c471c48306e4c36069404c5723f", "878662a229aaae906e123cdd9d3b4c10590ded29fe751eeeca34bbaa44af0773"],
            ["112233445566778899", "339150844ec15234807fe862a86be77977dbfb3ae3d96f4c22795513aeaab82f", "b1c14ddfdc8ec1b2583f51e85a5eb3a155840f2034730e9b5ada38b674336a21"],
            ["112233445566778899112233445566778899", "1b7e046a076cc25e6d7fa5003f6729f665cc3241b5adab12b498cd32f2803264", "bfea79be2b666b073db69a2a241adab0738fe9d2dd28b5604eb8c8cf097c457b"],
            ["29852220098221261079183923314599206100666902414330245206392788703677545185283", "9eace8f4b071e677c5350b02f2bb2b384aae89d58aa72ca97a170572e0fb222f", "1bbdaec2430b09b93f7cb08678636ce12eaafd58390699b5fd2f6e1188fc2a78"],
            ["57896042899961394862005778464643882389978449576758748073725983489954366354431", "878f22cc6db6048d2b767268f22ffad8e56ab8e2dc615f7bd89f1e350500dd8d", "714a5d7bb901c9c5853400d12341a892ef45d87fc553786756c4f0c9391d763e"],
            ["1766845392945710151501889105729049882997660004824848915955419660366636031", "659a379625ab122f2512b8dada02c6348d53b54452dff67ac7ace4e8856295ca", "49d81ab97b648464d0b4a288bd7818fab41a16426e943527c4fed8736c53d0f6"],
            ["28948025760307534517734791687894775804466072615242963443097661355606862201087", "cbceaaa8a4dd44bbce58e8db7740a5510ec2cb7ea8da8d8f036b3fb04cda4de4", "4bd7aa301a80d7f59fd983fedbe59bb7b2863fe46494935e3745b360e32332fa"]
        ],
        [
            ["1", "aa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7", "3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f"],
            ["2", "08d999057ba3d2d969260045c55b97f089025959a6f434d651d207d19fb96e9e4fe0e86ebe0e64f85b96a9c75295df61", "8e80f1fa5b1b3cedb7bfe8dffd6dba74b275d875bc6cc43e904e505f256ab4255ffd43e94d39e22d61501e700a940e80"],
            ["5", "11de24a2c251c777573cac5ea025e467f208e51dbff98fc54f6661cbe56583b037882f4a1ca297e60abcdbc3836d84bc", "8fa696c77440f92d0f5837e90a00e7c5284b447754d5dee88c986533b6901aeb3177686d0ae8fb33184414abe6c1713a"],
            ["10", "a669c5563bd67eec678d29d6ef4fde864f372d90b79b9e88931d5c29291238cced8e85ab507bf91aa9cb2d13186658fb", "a988b72ae7c1279f22d9083db5f0ecddf70119550c183c31c502df78c3b705a8296d8195248288d997784f6ab73a21dd"],
            ["112233445566778899", "a499efe48839bc3abcd1c5cedbdd51904f9514db44f4686db918983b0c9dc3aee05a88b72433e9515f91a329f5f4fa60", "3b7ca28ef31f809c2f1ba24aaed847d0f8b406a4b8968542de139db5828ca410e615d1182e25b91b1131e230b727d36a"],
            ["112233445566778899112233445566778899", "90a0b1cac601676b083f21e07bc7090a3390fe1b9c7f61d842d27fa315fb38d83667a11a71438773e483f2a114836b24", "3197d3c6123f0d6cd65d5f0de106fef36656cb16dc7cd1a6817eb1d51510135a8f492f72665cfd1053f75ed03a7d04c9"],
            ["10158184112867540819754776755819761756724522948540419979637868435924061464745859402573149498125806098880003248619520", "f2a066bd332dc59bbc3d01da1b124c687d8bb44611186422de94c1da4ecf150e664d353ccdb5cb2652685f8eb4d2cd49", "d6ed0bf75fdd8e53d87765fa746835b673881d6d1907163a2c43990d75b454294f942ec571ad5aae1806caf2bb8e9a4a"],
            ["9850501551105991028245052605056992139810094908912799254115847683881357749738726091734403950439157209401153690566655", "5c7f9845d1c4aa44747f9137b6f9c39b36b26b8a62e8af97290434d5f3b214f5a0131550adb19058dc4c8780c4165c4a", "712f7fccc86f647e70db8798228cb16344af3d00b139b6f8502939c2a965af0eb4e39e2e16ab8f597b8d5630a50c9d85"],
            ["9850502723405747097317271194763310482462751455185699630571661657946308788426092983270628740691202018691293898608608", "dd5838f7ec3b8acf1becfd746f8b668c577107e93548ed93ed0d254c112e76b10f053109ef8428bfcd50d38c4c030c57", "33244f479cdac34f160d9e4ce2d19d2ff0e3305b5bf0eef29e91e9de6e28f678c61b773aa7e3c03740e1a49d1aa2493c"],
            ["1146189371817832990947611400450889406070215735255370280811736587845016396640969656447803207438173695115264", "cb8ed893530bfba04b4ca655923aaad109a62bc8411d5925316c32d33602459c33057a1fbcb5f70aeb295d90f9165fbc", "426aee3e91b08420f9b357b66d5afcbcf3956590bf5564dbf9086042eb880493d19da39aaa6436c6b5fc66ce5596b43f"],
        ],
    ];

    #[test]
    fn scalar_mul_test() {
        for i in 1..2 {
            for vector in TEST_VECTORS[i].iter() {
                let ops = match i {
                    0 => &P256_COMMON_OPS,
                    1 => &P384_COMMON_OPS,
                    _ => panic!("blah")
                };

                let gen = AffinePoint {
                    x: elem_to_big_uint(ops, P384_GENERATOR.0),
                    y: elem_to_big_uint(ops, P384_GENERATOR.1),
                    curve_params: ops,
                    aux: &P384_AUX_OPS,
                };
                let k = biguint_to_scalar(gen.curve_params, &BigUint::parse_bytes(vector[0].as_bytes(), 16).unwrap());
                let test_p = AffinePoint {
                    x: BigUint::parse_bytes(vector[1].as_bytes(), 16).unwrap(),
                    y: BigUint::parse_bytes(vector[2].as_bytes(), 16).unwrap(),
                    curve_params: ops,
                    aux: &P384_AUX_OPS,
                };

                let aff = test_p.scalar_mul(&k).to_affine();
                let x_hex = aff.x.to_str_radix(16);
                let y_hex = aff.y.to_str_radix(16);

                assert_eq!(x_hex, vector[1], "x hex coordinate check for curve: {} and scalar: {}", i, vector[0]);
                assert_eq!(y_hex, vector[2], "y hex coordinate check for curve: {} and scalar: {}", i, vector[0]);
            }
        }
    }
}