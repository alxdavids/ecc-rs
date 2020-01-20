//! h2c module
use crate::ring_ecc;
use ring_ecc::ec::CurveID;
use ring_ecc::ec::suite_b::ops::Elem;
use ring_ecc::ec::suite_b::ops::PrivateKeyOps as CurveOps;
use ring_ecc::arithmetic::montgomery::R;

use num::{BigUint,BigInt,One};

use super::{AffinePoint,P256,P384,Encoded};
use super::utils;

pub struct HashToCurve {
	dst: String,
	z: u32,
	a: Elem<R>,
	b: Elem<R>,
	p: BigUint,
	m: i32,
    l: i32,
    h_eff: BigUint,
    id: CurveID,
	is_sq_exp: BigUint,
    sqrt_exp: BigUint,
    ops: CurveOps,
}

impl HashToCurve {
    fn new(id: CurveID, ops: CurveOps, modulus: BigUint) -> Self {
        let sqrt_exp = (&modulus+BigUint::from(1_u64))/BigUint::from(4_u64); // (p+1)/4
        let is_sq_exp = (&modulus-BigUint::from(1_u64))/BigUint::from(2_u64); // (p-1)/2
        match id {
            P256 => Self {
                dst: String::from("VOPRF-P384-SHA512-SSWU-RO-"),
                z: 10, // should actually be negative but made positive in SSWU anyway
                m: 1,
                l: 48,
                a: ops.common.a,
                b: ops.common.b,
                p: modulus,
                h_eff: BigUint::one(),
                id: id,
                sqrt_exp: sqrt_exp,
                is_sq_exp: is_sq_exp,
                ops: ops,
            },
            P384 => Self {
                dst: String::from("VOPRF-P384-SHA512-SSWU-RO-"),
                z: 12, // should actually be negative but made positive in SSWU anyway
                m: 1,
                l: 72,
                a: ops.common.a,
                b: ops.common.b,
                p: modulus,
                h_eff: BigUint::one(),
                id: id,
                sqrt_exp: sqrt_exp,
                is_sq_exp: is_sq_exp,
                ops: ops,
            },
            _ => panic!("unsupported curve")
        }
    }

    /// simplified swu
    fn sswu(&self, u_arr: Vec<BigUint>) -> AffinePoint<Encoded> {
        /// if length is greater than 1 then something bad is happening
        assert!(u_arr.len() == 1);
        let u = u_arr[0];
        let z = BigUint::from(self.z);

        // c1 = -b/a
        let a_inv = utils::invert_elem(&self.ops, self.a);
        let minus_a_inv = utils::minus_elem(&self.ops.common, &self.p, a_inv);
        let c1 = self.ops.common.elem_product(&self.b, &minus_a_inv);

        // c2 = -1/z
        // NOTE: we use -z in HashToCurve so we don't need to use minus again here
        let minus_z_inv = utils::invert_elem(&self.ops, utils::biguint_to_elem(self.ops.common, &z));
        AffinePoint::new()
    }
}

