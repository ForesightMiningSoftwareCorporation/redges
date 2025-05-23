//! Wedge data structure. Used for representing attribute vectors on the wedges surrounding
//! a vertex. Currently only useful for mesh simplification with attributes.
use std::collections::{BTreeMap, BTreeSet};

use linear_isomorphic::{InnerSpace, RealField};
use nalgebra::{ComplexField, DMatrix, DVector, Vector4};

use crate::{
    container_trait::{
        FaceAttributeGetter, FaceData, RedgeContainers, VertData, VertexAttributeGetter,
    },
    edge_handle::EdgeHandle,
    face_handle::{FaceHandle, FaceMetrics},
    vert_handle::VertHandle,
    FaceId, VertId,
};

/// A wedge is a collection of face corners incident on the same vertex, who all share the same
/// values.
#[derive(Debug, Clone)]
pub(crate) struct Wedge<S: RealField> {
    pub(crate) vertex: VertId,
    pub(crate) attributes: DVector<S>,
}

/// A wedge data structure allows us to collapse edges while preserving attribute discontinuities.
#[derive(Debug, Clone)]
pub(crate) struct WedgeDS<S: RealField> {
    pub(crate) wedges: Vec<Wedge<S>>,
    pub(crate) faces: BTreeMap<FaceId, BTreeMap<VertId, usize>>,
}

impl<S: RealField> Default for WedgeDS<S> {
    fn default() -> Self {
        Self::new()
    }
}

type WedgeQuadricOutput<S> = (
    DMatrix<S>,
    DVector<S>,
    S,
    DMatrix<S>,
    DVector<S>,
    S,
    Vec<(usize, usize)>,
    Vec<usize>,
);
impl<S: RealField> WedgeDS<S> {
    pub fn new() -> Self {
        Self {
            wedges: Vec::new(),
            faces: BTreeMap::new(),
        }
    }

    pub fn wedge_from_corner(&self, vert_id: VertId, face_id: FaceId) -> Option<Wedge<S>> {
        let list = self.faces.get(&face_id);
        list?;

        let list = list.unwrap();

        let wid = list.get(&vert_id);
        wid?;

        let wid = *wid.unwrap();

        Some(self.wedges[wid].clone())
    }

    pub fn wedge_id_from_corner(&self, vert_id: VertId, face_id: FaceId) -> Option<usize> {
        let list = self.faces.get(&face_id);
        list?;
        let list = list.unwrap();

        let wid = list.get(&vert_id);

        wid.copied()
    }

    pub fn wedge_quadric<R>(
        &self,
        edge: &EdgeHandle<R>,
        attribute_count: usize,
    ) -> WedgeQuadricOutput<S>
    where
        R: RedgeContainers,
        VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
        FaceData<R>: FaceAttributeGetter<S>,
        S: ComplexField,
    {
        // In non-manifold meshes, sometimes we generate weird edges that don't have incident faces.
        // but they could still connect to full faces whose quadrics we must optimize.
        let [f1, f2] = if edge.has_hedge() {
            let f1 = edge.hedge().face().id();
            let f2 = edge.hedge().radial_next().face().id();

            [f1, f2]
        } else {
            [FaceId::ABSENT, FaceId::ABSENT]
        };

        // Add the wedges of all faces touching the edge, except the two faces intersecting
        // at the edge.
        // TODO: what about the non-manifold case?
        let mut wedges = BTreeMap::new();
        for (vert, face) in edge
            .v1()
            .incident_faces()
            .map(|f| (edge.v1(), f))
            .chain(edge.v2().incident_faces().map(|f| (edge.v2(), f)))
            // Skip the two faces sharing the edge.
            .filter(|f| f.1.id() != f1 && f.1.id() != f2)
        {
            let (a, b, c) = face_geometric_quadric(&face);

            let wid = self.wedge_id_from_corner(vert.id(), face.id()).unwrap();
            let data = wedges.entry(wid).or_insert(Vec::new());

            let gradients = face_attribute_gradients(&face);
            let area = face.area();

            data.push(((a, b, c), gradients, area));
        }

        // For each side of the edge, see if wedges extend into the faces that will be removed.
        let merging_candidates = self.test_wedge_extensions(edge);

        let mut wedges_to_merge = Vec::new();
        // For each set of wedges that extend into the faces that will be deleted, merge them.
        for (w1, w2) in merging_candidates.into_iter().flatten() {
            // We probaly merged this wedge in a prior iteration.
            if !wedges.contains_key(&w1) || !wedges.contains_key(&w2) {
                continue;
            }
            // Move all data of wedge w2 into wedge w1, thus merging them.
            let wedges2 = wedges.get(&w2).unwrap().clone();
            let wedges1 = wedges.get_mut(&w1).unwrap();
            wedges1.extend(wedges2);
            wedges.remove(&w2);

            wedges_to_merge.push((w1, w2));
        }

        let wedge_count = wedges.len();
        // 3 because there are 3 coordinates in R^3.
        let dimension = 3 + wedge_count * attribute_count;
        let mut q = DMatrix::zeros(dimension, dimension);
        let mut final_b = DVector::zeros(dimension);
        let mut final_d = S::from(0.).unwrap();

        let mut geom_q = DMatrix::zeros(3, 3);
        let mut geom_b = DVector::zeros(3);
        let mut geom_d = S::from(0.).unwrap();

        // For each wedge, add its geometric component to the top left matrix.
        let mut wedge_order = Vec::new();
        for (i, (wid, wedge_data)) in wedges.into_iter().enumerate() {
            wedge_order.push(wid);
            for datum in wedge_data {
                // Add the geometric contributions of this face to the final quadric.
                let ((a, b, c), gradients, area) = datum;

                let mut top_left = q.view_mut((0, 0), (3, 3));
                top_left += a * area;
                geom_q += a * area;

                let mut top_three = final_b.rows_mut(0, 3);
                top_three += b * area;
                geom_b += b * area;

                final_d += c * area;
                geom_d += c * area;

                // Add the attribute information of this wedge to each section independently.
                // That is, each face contributes to the top left matrix and top three rows of the b vector.
                // but the attributes of each wedge are independent and thus assigned their own k regions
                // in both the matrix and the vector.
                let attrib_count = gradients.len();
                for (j, (g, d)) in gradients.into_iter().enumerate() {
                    let mut top_left = q.view_mut((0, 0), (3, 3));
                    top_left += g * g.transpose() * area;

                    let mut top_three = final_b.rows_mut(0, 3);
                    top_three += g * d * area;

                    let mut sub_attributes = final_b.rows_mut(3 + i * attrib_count, attrib_count);
                    sub_attributes[j] -= d * area;

                    final_d += d * d * area;

                    let offset = 3 + i * attrib_count;
                    let mut grad_row = q.fixed_view_mut::<1, 3>(offset + j, 0);
                    grad_row -= g.transpose() * area;

                    let mut grad_col = q.fixed_view_mut::<3, 1>(0, offset + j);
                    grad_col -= g * area;

                    q[(offset + j, offset + j)] += area;
                }
            }
        }

        (
            q,
            final_b,
            final_d,
            geom_q,
            geom_b,
            geom_d,
            wedges_to_merge,
            wedge_order,
        )
    }

    pub fn test_wedge_extensions<R>(&self, edge: &EdgeHandle<R>) -> [Option<(usize, usize)>; 2]
    where
        R: RedgeContainers,
        VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
        FaceData<R>: FaceAttributeGetter<S>,
        S: ComplexField,
    {
        // TODO: hack for non manifold weird cases, we need to refactor this function.
        if !edge.has_hedge() {
            return [None, None];
        }

        let mut result = [None, None];
        let f = edge.hedge().face();

        // If one endpoint is in the boundary and the other isn't,
        // then no overlap can exist on the wedges touching this edge.
        if edge.v1().is_in_boundary() != edge.v2().is_in_boundary() {
            return [None, None];
        }
        // TODO: what about non-manifold cases.

        // Test if the wedges on the left side extend into the left face.
        {
            let fa = edge.hedge().face_next().radial_next().face();
            let fb = edge.hedge().face_prev().radial_next().face();

            let va = edge.hedge().dest();
            let vb = edge.hedge().source();

            let wid_fa = self.wedge_id_from_corner(va.id(), f.id()).unwrap();
            let wid_fb = self.wedge_id_from_corner(vb.id(), f.id()).unwrap();

            let wid_a = self.wedge_id_from_corner(va.id(), fa.id()).unwrap();
            let wid_b = self.wedge_id_from_corner(vb.id(), fb.id()).unwrap();

            // There are two wedges that extend into this face.
            if wid_fa == wid_a && wid_b == wid_fb {
                result[0] = Some((wid_fa, wid_fb))
            }
        }

        // Test if the wedges on the right side extend into the right face.
        {
            let fa = edge.hedge().radial_next().face_next().radial_next().face();
            let fb = edge.hedge().radial_next().face_prev().radial_next().face();

            let va = edge.hedge().radial_next().dest();
            let vb = edge.hedge().radial_next().source();

            let wid_fa = self.wedge_id_from_corner(va.id(), f.id()).unwrap();
            let wid_fb = self.wedge_id_from_corner(vb.id(), f.id()).unwrap();

            let wid_a = self.wedge_id_from_corner(va.id(), fa.id()).unwrap();
            let wid_b = self.wedge_id_from_corner(vb.id(), fb.id()).unwrap();

            // There are two wedges that extend into this face.
            if wid_fa == wid_a && wid_b == wid_fb {
                result[1] = Some((wid_fa, wid_fb));
            }
        }

        result
    }

    // TODO: experimental and untested.
    /// Warning: untested.
    pub fn _test_wedge_extensions_non_manifold<R>(&self, edge: &EdgeHandle<R>) -> Vec<Vec<usize>>
    where
        R: RedgeContainers,
        VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
        FaceData<R>: FaceAttributeGetter<S>,
        S: ComplexField,
    {
        if !edge.has_hedge() {
            return Vec::new();
        }

        let mut result = Vec::new();
        for h in edge.hedge().radial_loop() {
            let h_prev = h.face_prev();
            let h_next = h.face_next();
            let fid = h.face().id();

            let wid1 = self.wedge_id_from_corner(h.source().id(), fid).unwrap();
            let wid2 = self.wedge_id_from_corner(h.dest().id(), fid).unwrap();

            let mut merge_group_1 = Vec::new();
            for face in h_prev
                .radial_loop()
                .map(|h| h.face())
                .filter(|f| f.id() != fid)
            {
                let wid = self
                    .wedge_id_from_corner(h.source().id(), face.id())
                    .unwrap();
                if wid == wid1 {
                    merge_group_1.push(wid);
                }
            }

            let mut merge_group_2 = Vec::new();
            for face in h_next
                .radial_loop()
                .map(|h| h.face())
                .filter(|f| f.id() != fid)
            {
                let wid = self.wedge_id_from_corner(h.dest().id(), face.id()).unwrap();
                if wid == wid2 {
                    merge_group_2.push(wid);
                }
            }

            if !merge_group_1.is_empty() {
                result.push(merge_group_1);
            }

            if !merge_group_2.is_empty() {
                result.push(merge_group_2);
            }
        }

        result
    }

    pub fn collapse_wedge<R>(
        &mut self,
        optimum: &DVector<S>,
        edge: &EdgeHandle<R>,
        wedges_to_merge: &Vec<(usize, usize)>,
        wedge_order: &[usize],
        attribute_count: usize,
    ) where
        R: RedgeContainers,
        VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
        FaceData<R>: FaceAttributeGetter<S>,
        S: ComplexField,
    {
        for (w1, w2) in wedges_to_merge {
            // Update all face corners pointing to w2, to now point to w1.
            let v1 = edge.v1();
            for face in v1.incident_faces() {
                let wid = self.wedge_id_from_corner(v1.id(), face.id()).unwrap();
                if wid == *w2 {
                    let verts = self.faces.get_mut(&face.id()).unwrap();
                    *verts.get_mut(&v1.id()).unwrap() = *w1;
                }
            }

            let v2 = edge.v2();
            for face in v2.incident_faces() {
                let wid = self.wedge_id_from_corner(v2.id(), face.id()).unwrap();
                if wid == *w2 {
                    let verts = self.faces.get_mut(&face.id()).unwrap();
                    *verts.get_mut(&v2.id()).unwrap() = *w1;
                }
            }

            // Invalidate all data in w2.
            self.wedges[*w2].vertex = VertId::ABSENT;
            self.wedges[*w2].attributes = DVector::default();
        }

        // Update the values assigned to the wedge, be careful to follow the order returned by the edge cost function.
        for (i, wid) in wedge_order.iter().enumerate() {
            let offset = 3 + i * attribute_count;
            let wedge = &mut self.wedges[*wid];
            for j in 0..attribute_count {
                wedge.attributes[j] = optimum[offset + j];
            }
        }
    }

    // TODO: This assumes manifold opology.
    pub fn insert_face_attributes<R>(&mut self, face: FaceHandle<R>)
    where
        R: RedgeContainers,
        VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
        FaceData<R>: FaceAttributeGetter<S>,
        S: ComplexField,
    {
        for vert_handle in face.vertices() {
            let n = face.data().attribute_count();
            if n == 0 {
                break;
            }

            let mut attribute = DVector::zeros(n);
            let vindex = face.data().inner_index(vert_handle.id());
            for i in 0..n {
                attribute[i] = face.data().attribute(vindex, i);
            }

            // Get the two faces adjacent to the current face, incident on the current vertex.
            let hedge = face
                .hedge()
                .face_loop()
                .find(|h| h.source().id() == vert_handle.id())
                .map(|h| h.id())
                .unwrap();
            let hedge = face.redge.hedge_handle(hedge);

            let f_prev = hedge.radial_next().face();
            let f_next = hedge.face_prev().radial_next().face();
            debug_assert!(f_prev.id() != face.id());
            debug_assert!(f_next.id() != face.id());
            debug_assert!(f_prev.id() != f_next.id());

            // See if the attribute of the prev face is contiguous with my attribute.
            let mut prev_wid = usize::MAX;
            if let Some(face_list) = self.faces.get(&f_prev.id()) {
                if let Some(wid) = face_list.get(&vert_handle.id()) {
                    let prev_attrib = self.wedges[*wid].attributes.clone();
                    let d = (&attribute - &prev_attrib).norm();

                    if S::from_real(d) < S::from(0.0001).unwrap() {
                        prev_wid = *wid;
                    }
                }
            }

            // See if the attribute of the next face is contiguous with my attribute.
            let mut next_wid = usize::MAX;
            if let Some(face_list) = self.faces.get(&f_next.id()) {
                if let Some(wid) = face_list.get(&vert_handle.id()) {
                    let next_attrib = self.wedges[*wid].attributes.clone();
                    let d = (&attribute - &next_attrib).norm();

                    if S::from_real(d) < S::from(0.0001).unwrap() {
                        next_wid = *wid;
                    }
                }
            }

            match (prev_wid != usize::MAX, next_wid != usize::MAX) {
                // Neither of my neighbours had a matching wedge, make a new one.
                (false, false) => {
                    let wid = self.wedges.len();
                    self.wedges.push(Wedge {
                        vertex: vert_handle.id(),
                        attributes: attribute,
                    });

                    let map = self.faces.entry(face.id()).or_default();
                    map.insert(vert_handle.id(), wid);
                }
                // Only my prev neighbour had a matching wedge, so assign myself to it.
                (true, false) => {
                    let map = self.faces.entry(face.id()).or_default();
                    map.insert(vert_handle.id(), prev_wid);
                }
                // Only my next neighbour had a matching wedge, so assign myself to it.
                (false, true) => {
                    let map = self.faces.entry(face.id()).or_default();
                    map.insert(vert_handle.id(), next_wid);
                }
                // Both neighbours had matching edges, assign and merge if necessary.
                (true, true) => {
                    // Assign to the prev wedge.
                    let map = self.faces.entry(face.id()).or_default();
                    map.insert(vert_handle.id(), prev_wid);

                    // My neighbours' wedges don't match each other, so merge them.
                    if prev_wid != next_wid {
                        for face in vert_handle.incident_faces() {
                            // Re-assign all faces pointing to `next_wid` to `prev_wid` instead.
                            if let Some(fwid) =
                                self.wedge_id_from_corner(vert_handle.id(), face.id())
                            {
                                if fwid == next_wid {
                                    let list = self.faces.get_mut(&face.id()).unwrap();
                                    let corner = list.get_mut(&vert_handle.id()).unwrap();
                                    *corner = prev_wid;
                                }
                            }
                        }

                        // Invalidate this wedge, but don't erase it because it would shift indices.
                        self.wedges[next_wid].vertex = VertId::ABSENT;
                        self.wedges[next_wid].attributes = DVector::default();
                    }
                }
            };
        }
    }

    pub fn vertex_wedges<R>(&self, vert: &VertHandle<R>) -> BTreeSet<usize>
    where
        R: RedgeContainers,
        VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
        FaceData<R>: FaceAttributeGetter<S>,
        S: ComplexField,
    {
        let mut wids = BTreeSet::new();
        for face in vert.incident_faces() {
            if let Some(wid) = self.wedge_id_from_corner(vert.id(), face.id()) {
                wids.insert(wid);
            }
        }

        wids
    }
}

/// Construct a quadric error metric matrix out from a face.
pub fn face_geometric_quadric<R: RedgeContainers, S>(
    face: &FaceHandle<R>,
) -> (nalgebra::Matrix3<S>, nalgebra::Vector3<S>, S)
where
    S: RealField + ComplexField,
    FaceData<R>: FaceAttributeGetter<S>,
    VertData<R>: InnerSpace<S>,
{
    let mut n = face.unit_normal();
    // Numeric hack to prevent nans.
    if !n[0].is_finite() || !n[1].is_finite() || !n[2].is_finite() {
        n[0] = S::from(1.).unwrap();
        n[1] = S::from(1.).unwrap();
        n[2] = S::from(1.).unwrap();
        n = n.normalized();
    }

    let n = nalgebra::Vector3::new(n[0], n[1], n[2]);

    let p = face.vertices().next().unwrap().data().clone();
    let p = nalgebra::Vector3::new(p[0], p[1], p[2]);

    let d = -p.dot(&n);

    (n * n.transpose(), n * d, d * d)
}

fn face_attribute_gradients<R: RedgeContainers, S>(
    face: &FaceHandle<R>,
) -> Vec<(nalgebra::Vector3<S>, S)>
where
    S: RealField + ComplexField,
    FaceData<R>: FaceAttributeGetter<S>,
    VertData<R>: InnerSpace<S>,
{
    assert!(face.is_active());
    let vids = face.data().attribute_vertices();

    // Order is REALLY important. We should always use the face data ordering as the canonical ordering
    // (as opposed to the implicit orderign given by the redge topology). Be very careful when modifying this.
    let positions: Vec<_> = vids
        .iter()
        .map(|id| {
            let d = face.vertex_by_id(*id).unwrap().data().clone();
            [d[0], d[1], d[2]]
        })
        .collect();

    let mut n = face.unit_normal();
    // Numeric hack to prevent nans.
    if !n[0].is_finite() || !n[1].is_finite() || !n[2].is_finite() {
        n[0] = S::from(1.).unwrap();
        n[1] = S::from(1.).unwrap();
        n[2] = S::from(1.).unwrap();
        n = n.normalized();
    }

    let n = nalgebra::Vector3::new(n[0], n[1], n[2]);

    let face_attribute_count = face.data().attribute_count();

    let mut face_attribs = [
        DVector::zeros(face_attribute_count),
        DVector::zeros(face_attribute_count),
        DVector::zeros(face_attribute_count),
    ];
    for (i, face_attrib) in face_attribs.iter_mut().enumerate() {
        for k in 0..face_attribute_count {
            face_attrib[k] = face.data().attribute(i, k);
        }
    }

    let mut gradients = Vec::new();
    for k in 0..face_attribute_count {
        let (g, d) = attribute_gradient(&positions, &n, &face_attribs, k);

        gradients.push((g, d));
    }

    gradients
}

fn attribute_gradient<S>(
    positions: &[[S; 3]],
    normal: &nalgebra::Vector3<S>,
    attributes: &[DVector<S>; 3],
    atrribute_id: usize,
) -> (nalgebra::Vector3<S>, S)
where
    S: RealField + ComplexField,
{
    let mut mat = nalgebra::Matrix4::<S>::zeros();
    let mut b = nalgebra::Vector4::<S>::zeros();
    for (i, attribs) in attributes.iter().enumerate() {
        debug_assert!(i < 3);

        mat[(i, 0)] = positions[i][0];
        mat[(i, 1)] = positions[i][1];
        mat[(i, 2)] = positions[i][2];
        mat[(i, 3)] = S::from(1.).unwrap();

        // TODO: Clamp is for uvs, and normals, should not be done in general.
        b[i] = attribs[atrribute_id].clamp(S::from(-1.0).unwrap(), S::from(1.0).unwrap());
    }

    mat[(3, 0)] = normal[0];
    mat[(3, 1)] = normal[1];
    mat[(3, 2)] = normal[2];

    // Trivial solution.
    if b == Vector4::zeros() {
        return (b.fixed_rows::<3>(0).into(), S::from(0.0).unwrap());
    }

    if let Some(inv) = mat.try_inverse() {
        let res = inv * b;
        (res.fixed_rows::<3>(0).into(), res[3])
    } else {
        // This is a hack. There are better ways to get a value from an undetermined system, however.
        (nalgebra::Vector3::default(), S::from(1.0).unwrap())
    }
}
