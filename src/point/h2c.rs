//! h2c module
use crate::ring_ecc;
use ring_ecc::ec::CurveID;
use ring_ecc::ec::suite_b::ops::Elem;
use ring_ecc::ec::suite_b::ops::PrivateKeyOps as CurveOps;
use ring_ecc::arithmetic::montgomery::R;

use num::{BigUint,One,Zero};

use super::{AffinePoint,P256,P384,Encoded};
use super::utils;
use super::utils::hkdf;

pub struct HashToCurve {
	dst: &'static str,
	z: Elem<R>,
	a: Elem<R>,
	b: Elem<R>,
	p: BigUint,
	m: usize,
    l: usize,
    h_eff: BigUint,
    id: CurveID,
	is_sq_exp: BigUint,
    sqrt_exp: BigUint,
    ops: &'static CurveOps,
}

impl HashToCurve {
    /// creates a new HashToCurve object
    pub fn new(id: CurveID, ops: &'static CurveOps, modulus: BigUint) -> Self {
        let sqrt_exp = (&modulus+BigUint::from(1_u64))/BigUint::from(4_u64); // (p+1)/4
        let is_sq_exp = (&modulus-BigUint::from(1_u64))/BigUint::from(2_u64); // (p-1)/2
        match id {
            P256 => Self {
                dst: "VOPRF-P384-SHA512-SSWU-RO-",
                z: Self::get_z_value(id, ops, &modulus, 10),
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
                dst: "VOPRF-P384-SHA512-SSWU-RO-",
                z: Self::get_z_value(id, ops, &modulus, 12),
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

    /// performs the hash_to_base algorithm, as specified in
    /// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-05#section-5
    fn hash_to_base(&self, msg: &[u8], ctr: u8) -> Vec<BigUint> {
        // 1. msg_prime = HKDF-Extract(DST, msg || I2OSP(0, 1))
        let mut exp_inp = msg.to_vec();
        let mut oct = vec![0; 1];
        utils::I2osp(0, 1, &mut oct);
        exp_inp.extend_from_slice(&oct);
        let mut msg_prime = Vec::new();
        hkdf::extract(self.dst.as_bytes(), &exp_inp, &mut msg_prime);

        // 2. info_pfx = "H2C" || I2OSP(ctr, 1)   // "H2C" is a 3-byte ASCII string
        let mut info_pfx = Vec::new();
        let mut ctr_buf = vec![0; 1];
        utils::I2osp(ctr as u64, 1, &mut ctr_buf);
        info_pfx.extend_from_slice("H2C".as_bytes());
        info_pfx.extend_from_slice(&ctr_buf);

        // 3.
        let mut u = Vec::new();
        for i in 1..self.m+1 {
            // 4. info = info_pfx || I2OSP(i, 1)
            let mut i_buf = vec![0; 1];
            utils::I2osp(i as u64, 1, &mut i_buf);
            let mut info = Vec::new();
            info.extend_from_slice(&info_pfx);
            info.extend_from_slice(&i_buf);

            // 5. t = HKDF-Expand(msg_prime, info, L)
            let mut t = vec![0; self.l];
            hkdf::expand(&msg_prime, &info, &mut t);

            // 6. e_i = OS2IP(t) mod p
            let e_i = utils::Os2ip(&t) % &self.p;
            u.push(e_i);
        }

        // 7. u = (e_1, ..., e_m)
        u
    }

    /// simplified swu as described in
    /// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-05#section-6.6.2
    fn sswu(&self, u_vec: &[BigUint]) -> AffinePoint<Encoded> {
        let ops = self.ops;
        let cops = self.ops.common;
        // if length is greater than 1 then something bad is happening
        assert!(u_vec.len() == 1);
        let u = utils::biguint_to_elem_unenc(self.id, cops, &u_vec[0]);

        // c1 = -b/a
        let a_inv = utils::invert_elem(&ops, self.a);
        let minus_a_inv = utils::minus_elem(&cops, &self.p, a_inv);
        let c1 = cops.elem_product(&self.b, &minus_a_inv);

        // c2 = -1/z
        let z_inv = utils::invert_elem(&ops, self.z);
        let c2 = utils::minus_elem(&cops, &self.p, z_inv);

        // 1. t1 = z * u^2
        let uu = cops.elem_squared(&u);
        let t1 = cops.elem_product(&self.z, &uu);

        let mut t2 = cops.elem_squared(&t1); // 2. t2 = t1^2

        // 3. x1 = t1 + t2
        let mut x1 = t1;
        cops.elem_add(&mut x1, &t2);

        x1 = utils::invert_elem(&ops, x1); // 4. x1 = x1^{-1}
        let e1 = utils::elem_to_biguint(x1) == Zero::zero(); // 5. e1 = (x1 == 0)
        cops.elem_add(&mut x1, &utils::elem_one(self.id, cops)); // 6. x1 = x1 + 1
        x1 = utils::elem_cmov(cops, x1, c2, e1); // 7. x1 = cmov(x1, c2, e1)
        x1 = cops.elem_product(&x1, &c1); // 8. x1 = x1 * c1
        let mut gx1 = cops.elem_squared(&x1); // 9. gx1 = x1^2
        cops.elem_add(&mut gx1, &self.a); // 10. gx1 = gx1 + a
        gx1 = cops.elem_product(&gx1, &x1); // 11. gx1 = gx1 * x1
        cops.elem_add(&mut gx1, &self.b); // 12. gx1 = gx1 + b
        let x2 = cops.elem_product(&t1, &x1); // 13. x2 = t1 * x1
        t2 = cops.elem_product(&t1, &t2); // 14. t2 = t1 * t2
        let gx2 = cops.elem_product(&gx1, &t2); // 15. gx2 = gx1 * t2
        let e2 = utils::is_square(&utils::elem_to_biguint(gx1),
                            &self.is_sq_exp, &self.p); // 16. e2 = is_square(gx1)
        let x = utils::elem_cmov(cops, x2, x1, e2); // 17. x = cmov(x2, x1, e2)
        let yy = utils::elem_cmov(cops, gx2, gx1, e2); // 18. x = cmov(gx2, gx1, e2)
        let mut y = utils::elem_sqrt(self.id, cops, yy, &self.sqrt_exp, &self.p); // 19. y = sqrt(yy)

        // 20. sgn0(u) == sgn0(y)
        let u_sgn = utils::sgn0_le_elem(cops, u);
        let y_sgn = utils::sgn0_le_elem(cops, y);
        let e3 = u_sgn == y_sgn;

        y = utils::elem_cmov(cops, utils::minus_elem(cops, &self.p, y), y, e3); // 21. y = cmov(-y, y, e3x)

        // construct point output object
        let mut point = AffinePoint::new(self.id);
        point.x = utils::elem_to_biguint(x);
        point.y = utils::elem_to_biguint(y);
        assert!(point.is_valid());
        point
    }

    /// returns the appropriate z value by creating an element and using the
    /// negative version. The reason that I do it like this is because I don't
    /// want to get involved with BigInt.
    ///
    /// TODO: uses BigUint ops
    fn get_z_value(id: CurveID, ops: &CurveOps, p: &BigUint, minus_z: u32) -> Elem<R> {
        utils::minus_elem(ops.common, p, utils::biguint_to_elem_unenc(id, ops.common, &BigUint::from(minus_z)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ring_ecc::ec::suite_b::ops::p384::PRIVATE_KEY_OPS as P384_OPS;

    #[test]
    fn h2b_test() {
        let h2c = HashToCurve::new(P384, &P384_OPS, utils::get_modulus_as_biguint(&P384_OPS.common.q));
        let mut count = 0;
        for vector in HASH_TO_BASE_VECTORS.iter() {
            let input = vector[0].as_bytes();
            let expected = BigUint::parse_bytes(vector[1].as_bytes(), 10).unwrap();
            let u_vec = h2c.hash_to_base(input, 0);
            assert_eq!(u_vec[0].to_bytes_be(), expected.to_bytes_be(), "h2b test for count: {:?}", count);
            count = count+1;
        }
    }

    #[test]
    fn sswu_test() {
        let h2c = HashToCurve::new(P384, &P384_OPS, utils::get_modulus_as_biguint(&P384_OPS.common.q));
        let mut count = 0;
        for vector in SSWU_TEST_VECTORS.iter() {
            let input = BigUint::parse_bytes(vector[0].as_bytes(), 10).unwrap();
            let exp_x = BigUint::parse_bytes(vector[1].as_bytes(), 10).unwrap();
            let exp_y = BigUint::parse_bytes(vector[2].as_bytes(), 10).unwrap();
            let u_vec = vec!(input);
            let point = h2c.sswu(&u_vec).to_unencoded();
            assert_eq!(point.x.to_bytes_be(), exp_x.to_bytes_be(), "x test for count: {:?}", count);
            assert_eq!(point.y.to_bytes_be(), exp_y.to_bytes_be(), "y test for count: {:?}", count);
            count = count+1;
        }
    }

    pub const HASH_TO_BASE_VECTORS: [[&str; 2]; 5] = [
        ["", "15670280948239665018787050025088822552903093865230238970017602952833555416398748331082295637805213707088989441755988"],
		["1", "1942715482632358166165565369095283869513634648389774012602448122359464835733690346035199729746417427046377204715303"],
		["asdf", "24507112164256266255100924053603326775213507976390981967792131453083876194411216719447408537203841824718570787142464"],
		["test", "6409376039185531560017287982748544597515854411296193693488280424481644496093326544690902528863962436268623496771541"],
		["random", "16247250678686872222869936093984092594492729196895879130498408114251281419554923530849483086336127849429159109128818"],
    ];

    pub const SSWU_TEST_VECTORS: [[&str; 3]; 5] = [
        [
           "15670280948239665018787050025088822552903093865230238970017602952833555416398748331082295637805213707088989441755988",
           "15043091655123589139476535520316853145074562564067200072853707836963164937518115020044315814573473606362869394777187",
           "33136250779564189967647894388148954739171786982148515795299597591669909884906353483262749333579115340263403769866626",
        ],
        [
           "1942715482632358166165565369095283869513634648389774012602448122359464835733690346035199729746417427046377204715303",
           "31712666608794813838450831245768352608061820731219254600083599907999316691595120493791689938289282078353090933837041",
           "4206609551883326717841767788616124592725060605241985514692257065399455065867170452124896082541316335165985229730507",
        ],
        [
           "24507112164256266255100924053603326775213507976390981967792131453083876194411216719447408537203841824718570787142464",
           "293447988360561042611832928522597727479089496568847551940003813796772318506727270172476083341418873730785454701568",
           "9761653435465566913614766398945376337690717880421211083207566458447975999387790669048364346983316135762944196549898",
        ],
        [
           "6409376039185531560017287982748544597515854411296193693488280424481644496093326544690902528863962436268623496771541",
           "31475109408547543147199457632396492796169708514999370150255421041761202109773477769740788427961884005032653705307760",
           "33864641941997256225600777383976921381186308220482560482046046771370099718859942454304724377098962437161537347141181",
        ],
        [
           "16247250678686872222869936093984092594492729196895879130498408114251281419554923530849483086336127849429159109128818",
           "37310690097326955526874904412484957930185253899573931562503850700495595096900370544619577329930240463109907820476327",
           "1835534565261005396339419959321852444925359938134835077886551094760311758797304939482177839199023849240736081211984",
        ],
    ];
}

