use std::ops::{Add, AddAssign, Mul};

use linear_isomorphic::prelude::*;

#[derive(Default, Debug, Clone, Copy, PartialEq)]
pub struct Quadric<S, V>
where
    S: RealField + Mul<V, Output = V>,
    V: InnerSpace<S>,
{
    pub nxx: S,
    pub nyy: S,
    pub nzz: S,
    pub nxy: S,
    pub nxz: S,
    pub nyz: S,

    pub dn: V,
    pub d2: S,
    pub area: S,
}

impl<S, V> Add for Quadric<S, V>
where
    S: RealField + Mul<V, Output = V>,
    V: InnerSpace<S>,
{
    type Output = Self;

    fn add(mut self, other: Self) -> Self::Output {
        self.nxx += other.nxx;
        self.nyy += other.nyy;
        self.nzz += other.nzz;
        self.nxy += other.nxy;
        self.nxz += other.nxz;
        self.nyz += other.nyz;

        self.dn += other.dn;
        self.d2 += other.d2;
        self.area += other.area;

        self
    }
}

impl<S, V> AddAssign for Quadric<S, V>
where
    S: RealField + Mul<V, Output = V>,
    V: InnerSpace<S>,
{
    fn add_assign(&mut self, other: Self) {
        self.nxx += other.nxx;
        self.nyy += other.nyy;
        self.nzz += other.nzz;
        self.nxy += other.nxy;
        self.nxz += other.nxz;
        self.nyz += other.nyz;

        self.dn += other.dn;
        self.d2 += other.d2;
        self.area += other.area;
    }
}

impl<S, V> Quadric<S, V>
where
    S: RealField + Mul<V, Output = V>,
    V: InnerSpace<S>,
{
    pub fn add_border_quadric(&mut self, border: BorderQuadric<S, V>, p: V) {
        let dist = -border.n.dot(&p);

        self.nxx += border.nxx;
        self.nyy += border.nyy;
        self.nzz += border.nzz;
        self.nxy += border.nxy;
        self.nxz += border.nxz;
        self.nyz += border.nyz;

        let weighted_dist = border.total_weight * dist;

        self.dn += -border.total_weight * p.clone() - weighted_dist * border.n;
        self.d2 += border.total_weight * p.norm_squared() - weighted_dist * weighted_dist;
    }

    pub fn from_plane(p: V, n: V, weight: S) -> Self {
        let pn = p.dot(&n);
        Self {
            nxx: n[0] * n[0],
            nyy: n[1] * n[1],
            nzz: n[2] * n[2],
            nxy: n[0] * n[1],
            nxz: n[0] * n[2],
            nyz: n[1] * n[2],
            dn: n * pn,
            d2: pn * pn,
            area: weight,
        }
    }

    pub fn from_tri(v0: V, v1: V, v2: V) -> Self {
        let v01 = v1 - v0.clone();
        let v02 = v2 - v0.clone();

        let mut n = v01.cross(&v02);
        let mag = n.norm();
        let area = mag * S::from(0.5).unwrap();

        if is_degenerate_divisor(area) {
            return Self::default();
        }

        n = n.normalized();

        let distance = n.dot(&v0);

        Self {
            nxx: n[0] * n[0],
            nyy: n[1] * n[1],
            nzz: n[2] * n[2],
            nxy: n[0] * n[1],
            nxz: n[0] * n[2],
            nyz: n[1] * n[2],
            dn: n * distance,
            d2: distance * distance,
            area,
        }
    }

    pub fn error(&self, point: V) -> S {
        let mut v1 = V::default();
        let mut v2 = V::default();
        let mut v3 = V::default();
        v1.set_subset(&[self.nxx, self.nxy, self.nxz]);
        v2.set_subset(&[self.nxy, self.nyy, self.nyz]);
        v3.set_subset(&[self.nxz, self.nyz, self.nzz]);

        let x = point.dot(&v1);
        let y = point.dot(&v2);
        let z = point.dot(&v3);

        let mut v = V::default();
        v.set_subset(&[x, y, z]);
        let v_av = point.dot(&v);
        let b_tv = -point.dot(&self.dn);

        let q = v_av + S::from(2.0).unwrap() * b_tv + self.d2;

        let zero = S::from(0.0).unwrap();
        if is_degenerate(q) {
            zero
        } else {
            q.max(zero) * self.area
        }
    }

    pub fn optimize(&self) -> Option<(S, V)> {
        let a = self.nxx;
        let b = self.nxy;
        let c = self.nxz;
        let d = self.nyy;
        let e = self.nyz;
        let f = self.nzz;
        let r0 = self.dn[0];
        let r1 = self.dn[1];
        let r2 = self.dn[2];

        let ad = a * d;
        let ae = a * e;
        let af = a * f;
        let bc = b * c;
        let be = b * e;
        let bf = b * f;
        let df = d * f;
        let ce = c * e;
        let cd = c * d;

        let be_cd = be - cd;
        let bc_ae = bc - ae;
        let ce_bf = ce - bf;

        let det = a * df + S::from(2.0).unwrap() * b * ce - ae * e - bf * b - cd * c;
        if is_degenerate_divisor(det) {
            return None;
        }
        let denom = det;
        let nom0 = r0 * (df - e * e) + r1 * ce_bf + r2 * be_cd;
        let nom1 = r0 * ce_bf + r1 * (af - c * c) + r2 * bc_ae;
        let nom2 = r0 * be_cd + r1 * bc_ae + r2 * (ad - b * b);

        let mut point = V::default();
        point.set_subset(&[nom0, nom1, nom2]);
        point[0] = point[0] / denom;
        point[1] = point[1] / denom;
        point[2] = point[2] / denom;

        Some((self.error(point.clone()), point))
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq)]
pub struct BorderQuadric<S, V>
where
    S: RealField + Mul<V, Output = V>,
    V: InnerSpace<S>,
{
    pub nxx: S,
    pub nyy: S,
    pub nzz: S,
    pub nxy: S,
    pub nxz: S,
    pub nyz: S,

    pub n: V,
    // the weight of the border quadric plus the length of the edge
    pub total_weight: S,
}

impl<S, V> BorderQuadric<S, V>
where
    S: RealField + Mul<V, Output = V>,
    V: InnerSpace<S>,
{
    pub fn from_edge(weight: S, v0: V, v1: V, tri_normal: V) -> Self {
        let v01 = v1 - v0;
        let mut n = tri_normal.cross(&v01);
        let mag = n.norm();

        if is_degenerate_divisor(mag) {
            return Self::default();
        }

        n = n * (S::from(1.0).unwrap() / mag);

        let total_weight = weight * mag;

        Self {
            nxx: total_weight - n[0] * n[0] * total_weight,
            nyy: total_weight - n[1] * n[1] * total_weight,
            nzz: total_weight - n[2] * n[2] * total_weight,
            nxy: n[0] * n[1] * -total_weight,
            nxz: n[0] * n[2] * -total_weight,
            nyz: n[1] * n[2] * -total_weight,
            n,
            total_weight,
        }
    }
}

fn is_degenerate<S: RealField>(f: S) -> bool {
    f.is_nan() || f.is_infinite()
}
fn is_degenerate_divisor<S: RealField>(f: S) -> bool {
    is_degenerate(f) || f.abs() < S::from(1e-8).unwrap()
}
