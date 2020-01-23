//! h2c module
use crate::ring_ecc;
use ring_ecc::ec::CurveID;
use ring_ecc::ec::suite_b::ops::Elem;
use ring_ecc::ec::suite_b::ops::PrivateKeyOps as CurveOps;
use ring_ecc::arithmetic::montgomery::R;

use num::{BigUint,One,Zero};
use std::io::Error;

use super::{AffinePoint,P256,P384,Encoded};
use super::utils;
use super::errors;
use hkdf_sha512 as hkdf;

pub const SSWU_RO: &'static str = "SSWU-RO";

pub struct HashToCurve {
    mapping_name: &'static str,
	dst: String,
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
    pub fn new(id: CurveID, ops: &'static CurveOps, modulus: BigUint) -> Result<Self, Error> {
        let sqrt_exp = (&modulus+BigUint::from(1_u64))/BigUint::from(4_u64); // (p+1)/4
        let is_sq_exp = (&modulus-BigUint::from(1_u64))/BigUint::from(2_u64); // (p-1)/2
        match id {
            P256 => Ok(
                Self {
                    mapping_name: SSWU_RO,
                    dst: format!("VOPRF-P256-SHA512-{}-", SSWU_RO),
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
                }
            ),
            P384 => Ok(
                Self {
                    mapping_name: SSWU_RO,
                    dst: format!("VOPRF-P384-SHA512-{}-", SSWU_RO),
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
                }
            ),
            _ => Err(errors::unsupported())
        }
    }

    pub fn full(&self, alpha: &[u8]) -> AffinePoint<Encoded> {
        let map_to_curve = match self.mapping_name {
            "SSWU-RO" => Ok(|msg| self.sswu(msg)),
            _ => Err(errors::unsupported())
        }.unwrap();

        // See https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-05#section-6.6.1
        let u0 = self.hash_to_base(alpha, 0);
        let u1 = self.hash_to_base(alpha, 1);
        let Q0 = (map_to_curve)(&u0).to_jacobian();
        let Q1 = (map_to_curve)(&u1).to_jacobian();
        let R = Q0.add(&Q1).to_affine();

        // clear cofactor and return
        R.scalar_mul(&self.h_eff.to_bytes_be()).to_affine()
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
        let u = utils::biguint_to_elem_unenc(self.id, cops, &u_vec[0]).unwrap();

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
        let mut point = AffinePoint::new(self.id).unwrap();
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
        utils::minus_elem(ops.common, p,
                utils::biguint_to_elem_unenc(id, ops.common, &BigUint::from(minus_z)).unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ring_ecc::ec::suite_b::ops::p256::PRIVATE_KEY_OPS as P256_OPS;
    use ring_ecc::ec::suite_b::ops::p384::PRIVATE_KEY_OPS as P384_OPS;

    #[test]
    fn h2b_test() {
        for i in 0..2 {
            let h2c = match i {
                0 => HashToCurve::new(P256, &P256_OPS, utils::get_modulus_as_biguint(&P256_OPS.common)).unwrap(),
                1 => HashToCurve::new(P384, &P384_OPS, utils::get_modulus_as_biguint(&P384_OPS.common)).unwrap(),
                _ => panic!("bad i value"),
            };
            let mut count = 0;
            for vector in HASH_TO_BASE_VECTORS[i].iter() {
                let input = vector[0].as_bytes();
                let expected = BigUint::parse_bytes(vector[1].as_bytes(), 10).unwrap();
                let u_vec = h2c.hash_to_base(input, 0);
                assert_eq!(u_vec[0].to_bytes_be(), expected.to_bytes_be(), "h2b test for count: {:?} and curve: {:?}", count, h2c.id);
                count = count+1;
            }
        }
    }

    #[test]
    fn sswu_test() {
        for i in 0..2 {
            let h2c = match i {
                0 => HashToCurve::new(P256, &P256_OPS, utils::get_modulus_as_biguint(&P256_OPS.common)).unwrap(),
                1 => HashToCurve::new(P384, &P384_OPS, utils::get_modulus_as_biguint(&P384_OPS.common)).unwrap(),
                _ => panic!("bad i value"),
            };
            let mut count = 0;
            for vector in SSWU_TEST_VECTORS[i].iter() {
                let input = BigUint::parse_bytes(vector[0].as_bytes(), 10).unwrap();
                let exp_x = BigUint::parse_bytes(vector[1].as_bytes(), 10).unwrap();
                let exp_y = BigUint::parse_bytes(vector[2].as_bytes(), 10).unwrap();
                let u_vec = vec!(input);
                let point = h2c.sswu(&u_vec).to_unencoded();
                assert_eq!(point.x.to_bytes_be(), exp_x.to_bytes_be(), "x test for count: {:?} and curve: {:?}", count, h2c.id);
                assert_eq!(point.y.to_bytes_be(), exp_y.to_bytes_be(), "y test for count: {:?} and curve: {:?}", count, h2c.id);
                count = count+1;
            }
        }
    }

    #[test]
    fn full_test() {
        for i in 0..2 {
            let h2c = match i {
                0 => HashToCurve::new(P256, &P256_OPS, utils::get_modulus_as_biguint(&P256_OPS.common)).unwrap(),
                1 => HashToCurve::new(P384, &P384_OPS, utils::get_modulus_as_biguint(&P384_OPS.common)).unwrap(),
                _ => panic!("bad i value"),
            };
            let mut count = 0;
            for vector in FULL_HASH_TO_CURVE_VECTORS[i].iter() {
                let input = vector[0].as_bytes();
                let exp_x = BigUint::parse_bytes(vector[1].as_bytes(), 10).unwrap();
                let exp_y = BigUint::parse_bytes(vector[2].as_bytes(), 10).unwrap();
                let point = h2c.full(&input).to_unencoded();
                assert_eq!(point.x.to_bytes_be(), exp_x.to_bytes_be(), "x test for count: {:?} and curve: {:?}", count, h2c.id);
                assert_eq!(point.y.to_bytes_be(), exp_y.to_bytes_be(), "y test for count: {:?} and curve: {:?}", count, h2c.id);
                count = count+1;
            }
        }
    }


    // test vectors taken from
    // https://github.com/alxdavids/voprf-poc/blob/master/go/oprf/groups/ecgroup/h2c_test.go
    // (extra values derived for P-256).
    pub const HASH_TO_BASE_VECTORS: [[[&str; 2]; 5]; 2] = [
        [
            ["", "83535524130228921029437730219861701397353033315370087929938533023961338081610"],
            ["1", "24411817093339714607733915344869069814441911555609622374798957331470094523772"],
            ["asdf", "8517025996406018755637445955328551045057963231050648689430172083427787837670"],
            ["test", "42162176653912293023698044969243458081097858987743610568212464212738256993440"],
            ["random", "113589217908445930439596756311623491830232165921163730708362535829147959509856"],
        ],
        [
            ["", "15670280948239665018787050025088822552903093865230238970017602952833555416398748331082295637805213707088989441755988"],
            ["1", "1942715482632358166165565369095283869513634648389774012602448122359464835733690346035199729746417427046377204715303"],
            ["asdf", "24507112164256266255100924053603326775213507976390981967792131453083876194411216719447408537203841824718570787142464"],
            ["test", "6409376039185531560017287982748544597515854411296193693488280424481644496093326544690902528863962436268623496771541"],
            ["random", "16247250678686872222869936093984092594492729196895879130498408114251281419554923530849483086336127849429159109128818"],
        ]
    ];

    pub const SSWU_TEST_VECTORS: [[[&str; 3]; 5]; 2] = [
        [
            [
                "83535524130228921029437730219861701397353033315370087929938533023961338081610",
                "33246369943938658919323543545344077424579626066209059196315953968634806929589",
                "100579838009521953273512687994874777262635951008101990025114162974610192755778",
            ],
            [
                "24411817093339714607733915344869069814441911555609622374798957331470094523772",
                "712822355573976207315803288855643124830159649949246312338231401360734881035",
                "98802960212814468832734368047768081275713452913337463805082198558691662244642",
            ],
            [
                "8517025996406018755637445955328551045057963231050648689430172083427787837670",
                "47771283741417257855150510251356366045897902940740905499958561787800536897362",
                "54857595827007157103895576434859775056387098586181970767982835231065980002014",
            ],
            [
                "42162176653912293023698044969243458081097858987743610568212464212738256993440",
                "97348771341464391851577114286562352968120415428113837413159401842391823739253",
                "55393388684603922203260277905842915317789685517381484860940291161280052818602",
            ],
            [
                "113589217908445930439596756311623491830232165921163730708362535829147959509856",
                "25300982619573982394573922539431048395655419926103069164639074525314863918571",
                "8414461630540610217258268207762326490024997528417722370873173974260624106998",
            ],
        ],
        [
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
        ]
    ];

    pub const FULL_HASH_TO_CURVE_VECTORS: [[[&str; 3]; 5]; 2] = [
        [
            [
                "",
                "111016463477202659154256522758589203825459632420837451417967806756967750448452",
                "28947031354042851680314885070966712160557022699642034220096303355952942045989",
            ],
            [
                "1",
                "65568821604399219050273597479032438298872891216217300137665502250729949286199",
                "35147046718126604274424775790514260295309474447552727563844085072752648122008",
            ],
            [
                "asdf",
                "83230899211791055962927082806375996630304273504516139258674378179692817967234",
                "64838049334905204846189823100547327949475774292485779227251852745279450057765",
            ],
            [
                "test",
                "12450411665069793449868791705989898101059662143044179649910729165308035459023",
                "6602000077096978149184571662728829235076665216666459971093621024723181946870",
            ],
            [
                "random",
                "41990908553190152500965584160536237278250888305135761756480916601056399581307",
                "22340702353097352810037562436046093638816158965232419116385011882788505302027",
            ],
        ],
        [
            [
                "",
                "30080611775067838193475075004665419527937570396653956651519246592569896222441582047156381322632437363661635355059005",
                "20783652428854690810060531204648743925284619218538801076205938644463325455616369725402399433833851021039801860251878",
            ],
            [
                "1",
                "29650639659274268559136011553864194418207682311050323428173462440594796529912091771908609365933993237437866304383610",
                "3123044785607009045040711490412599434775424958141229760770582294918664212503090438091817468980366885107838379509098",
            ],
            [
                "asdf",
                "29969588127226151911382588418021312873012227179044443716367955445066566752849478826037970129940763289625714821443011",
                "17410069451102133321720859095615324374699361853698493986383537650837194987993565478405818082113600644209841551176018",
            ],
            [
                "test",
                "35545509722549146939660727050796900115452941653073989167788838920390302482004128874997252970331147295344750824226579",
                "27687865874587861560504941570144773748765286782363087596575730017098117334752417380625066587341878521688930304472033",
            ],
            [
                "random",
                "21634107956511686237571364665733337762160319624271853087423499351943896659075117271938533968539011259822501269661449",
                "24988434919453800168740599843788684084233292688651079542544278337911292329783936136946570334754766549649989066107853",
            ],
        ]
    ];
}

