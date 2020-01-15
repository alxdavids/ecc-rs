use num::{BigUint,Zero,One};
use crate::ring_ecc;
use ring_ecc::ec::CurveID;
use ring_ecc::ec::suite_b::ops::Point as RingPoint;
use ring_ecc::ec::suite_b::ops::{Elem,CommonOps};
use ring_ecc::ec::suite_b::ops::p256::COMMON_OPS as P256_COMMON_OPS;
use ring_ecc::ec::suite_b::ops::p384::COMMON_OPS as P384_COMMON_OPS;
use ring_ecc::limb::Limb;
use ring_ecc::ec::suite_b::ops::elem::MAX_LIMBS;
use ring_ecc::ec::suite_b::ops::elem_parse_big_endian_fixed_consttime;
use ring_ecc::arithmetic::montgomery::R;

use untrusted;
use core::marker::PhantomData;

/// AffinePoint
struct AffinePoint {
    x: BigUint,
    y: BigUint,
    ops: &'static CommonOps,
}

impl AffinePoint {
    /// Returns the identity point
    pub fn new(id: CurveID) -> Self {
        match id {
            CurveID::P256 => Self {
                x: Zero::zero(),
                y: Zero::zero(),
                ops: &P256_COMMON_OPS,
            },
            CurveID::P384 => Self {
                x: Zero::zero(),
                y: Zero::zero(),
                ops: &P384_COMMON_OPS,
            },
            _ => panic!("unsupported"),
        }
    }

    /// to jacobian
    pub fn to_jacobian(&self) -> JacobianPoint {
        JacobianPoint {
            x: self.x.clone(),
            y: self.y.clone(),
            z: One::one(),
            ops: &P256_COMMON_OPS,
        }
    }
}

struct JacobianPoint {
    x: BigUint,
    y: BigUint,
    z: BigUint,
    ops: &'static CommonOps,
}

impl JacobianPoint {
    /// Returns the identity point
    pub fn new(id: CurveID) -> Self {
        match id {
            CurveID::P256 => {
                AffinePoint {
                    x: Zero::zero(),
                    y: Zero::zero(),
                    ops: &P256_COMMON_OPS,
                }
            },
            CurveID::P384 => AffinePoint {
                x: Zero::zero(),
                y: Zero::zero(),
                ops: &P384_COMMON_OPS,
            },
            _ => panic!("unsupported"),
        }.to_jacobian()
    }

    pub fn add(&self, other: &JacobianPoint) -> JacobianPoint {
        let p1 = self.as_ring_jac_point();
        let p2 = other.as_ring_jac_point();
        let p_add = self.ops.point_sum(&p1, &p2);
        Self {
            x: Zero::zero(),
            y: Zero::zero(),
            z: Zero::zero(),
            ops: &P256_COMMON_OPS,
        }
    }

    fn as_ring_jac_point(&self) -> RingPoint {
        let (x, y, z) = self.encode_for_ops();
        let mut limbs: [Limb; MAX_LIMBS*3] = [0; MAX_LIMBS*3];
        limbs[..self.ops.num_limbs].copy_from_slice(&x.limbs);
        limbs[MAX_LIMBS..MAX_LIMBS+self.ops.num_limbs].copy_from_slice(&y.limbs);
        limbs[MAX_LIMBS*2..MAX_LIMBS*2+self.ops.num_limbs].copy_from_slice(&z.limbs);
        RingPoint { xyz: limbs }
    }

    fn from_ring_point(ring_pt: RingPoint, id: CurveID) -> Self {
        let pt = Self::new(id);
        let x = pt.ops.point_x(&ring_pt);
        let y = pt.ops.point_y(&ring_pt);
        let z = pt.ops.point_z(&ring_pt);
        // TODO
        Self {
            x: Zero::zero(),
            y: Zero::zero(),
            z: Zero::zero(),
            ops: &P256_COMMON_OPS,
        }
    }

    /// encodes the fields of the point
    fn encode_for_ops(&self) -> (Elem<R>, Elem<R>, Elem<R>) {
        let x_elem = bytes_to_elem(self.ops, &self.x.to_bytes_be());
        let y_elem = bytes_to_elem(self.ops, &self.y.to_bytes_be());
        let z_elem = bytes_to_elem(self.ops, &self.z.to_bytes_be());
        (x_elem, y_elem, z_elem)
    }
}

/// parses an elem
fn bytes_to_elem(ops: &CommonOps, bytes: &[u8]) -> Elem<R> {
    Elem {
        limbs: elem_parse_big_endian_fixed_consttime(ops, untrusted::Input::from(bytes)).unwrap().limbs,
        encoding: PhantomData,
        m: PhantomData
    }
}

/// back to bytes
fn elem_to_bytes() {}