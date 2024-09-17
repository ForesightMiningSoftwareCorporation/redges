#import "@preview/cheq:0.1.0": checklist

#align(center, text(17pt)[
  *Redge Manual*
])

#show: checklist

= Introduction

This document explains the theory behind the radial-edge data structure as implemented for Foresight Spatial Labs code base. It is focused on the general understanding of the structure and does not document the implementation itself, for details on available functionality please refer to the source code.

= The data structure

The Redge consists of essentially four metadata structures, delineated as follows:

- `Vertex`:
  - Pointer to an edge 
- `Edge`:
  - Two pointers to the end point vertices.
  - Pointer to a half edge in a face incident on the edge.
  - Pointers for a circular-doubly-linked-list of the edges around the first endpoint.
  - Pointers for a circular-doubly-linked-list of the edges around the second endpoint.
- `Half Edge`: 
  - Pointer to the source vertex of the half edge.
  - Pointer to the edge parallel to the half edge.
  - Pointer to the next half edge in the radial cycle of half edges orbitting the edge.
  - Pointer to the previous half edge in the radial cycle of half edges orbitting the edge.
  - Pointer to the next half edge in the face loop.
  - Pointer to the previous half edge in the face loop.
  - Pointer to the face this half edge belongs to.
- `Face`: 
  - Pointer to a half edge in the loop for this face.

@redge_overview Shows the general structure of a radial Edge and how each element points to the others.

#figure(
  image("images/redge.svg", width: 80%),
  caption: [Overview of the Radial Edge data structure. (Not all pointers are shown to reduce visual noise.)],
) <redge_overview>

= Topological Transformations

The basic topological transformations are the following:

- Vertex Removal: Removes a vertex and all edges and faces incident on it.
- Edge Flip: Changes the orientation of an edge to the edge crossing it transversally across the quad made by its two incident faces.
- Edge Split: Fractions an edge and introduces new faces, edges and a new vertex.
- Edge Collapse: Makes an edge disappear along with all the faces incident onto it and merges the vertices of the endpoints.
- Face Split: Introduces a new vertex inside of a face and adds as many faces and edges as necessary to connect the mesh to it.
- Face Collapse: Replaces a face with a vertex and modifies the surrounding topology accordingly.

However, not all of the above are commonly used, so this document will focus only on those which are necessary for our current needs. More might be added as it becomes necessary.

== Edge Flip

An edge flip only makes sense in the context of a triangle, manifold mesh. It grabs the set made by the two faces incident on the edge, and changes the connectivity of the edge so that it now connects the transversal vertices, as seen in @edge_flip.

#figure(
  grid(
    columns: 2,     // 2 means 2 auto-sized columns
    gutter: 2mm,    // space between columns
      image("images/edge_flip.svg", width: 90%),
    ),
  caption: [States for an edge set before and after an edge flip. Edges are not labeled for shortness, by convention $h_i$ points to $e_i$.],
) <edge_flip>

Note that there are many cases where an edge flip is not possible. The edge should be manifold, it can't be a boundary edge, and additionally, the transversal edges should not be connected by any edge.

== Edge Collapse

Due to the genericity of the radial edge, an edge collapse is better decomposed into a few sub-operations that must be used together to achieve the final result.

=== Edge Removal

Ignoring any modifications done to incident faces (which may not exist), removing an edge consists on removing one of the two vertices, removing the edge from the set *and* from the two cycles of its endpoints, then updating the necessary pointers. As shown in @edge_edge_collapse.

#figure(
  grid(
    columns: 2,     // 2 means 2 auto-sized columns
    gutter: 2mm,    // space between columns
      image("images/edge_edge_collapse_before.svg", width: 40%),
      image("images/edge_edge_collapse_after.svg", width: 40%),
    ),
  caption: [States for an edge set before and after an edge collapse.],
) <edge_edge_collapse>

A checklist of what must be done for this element:

- [ ] $v_2$ must be removed.
- [ ] $e_0$ must be removed from the set.
- [ ] The vertex cyles $C_1$ and $C_2$ of $e_0$'s endpoints must be updated such that $e_0$ is removed from both and they now form a single doubly-linked-list.
- [ ] $v_1$ must now point to a different edge than $e_0$ if possible. If no edge exists, then it must point to null.
- [ ] Every edge that was incident to $v_2$ must now point to $v_1$.

Note that this leaves the Radial Edge in a broken state, all operations in this section must be applied to restore correctness.

=== Face Removal

As shown in @face_edge_collapse, when collapsing an edge, we will want to modify the connectivity of the face complices surrounding it. For each face incident on the edge, there will be one half edge that must be removed and thus pointers must be updated.

#figure(
  grid(
    columns: 2,     // 2 means 2 auto-sized columns
    gutter: 2mm,    // space between columns
      image("images/face_edge_collapse_before.svg", width: 70%),
      image("images/face_edge_collapse_after.svg", width: 70%),
    ),
  caption: [State for a face set before and after an edge collapse.],
) <face_edge_collapse>

This is a checklist of all elements that need to be updated within the same face complex and ignoring the edge itself.

- [ ] Face $f$ must now point to either $h_p$ or $h_n$.
- [ ] $h_0$ must be removed from its radial cycle.
- [ ] $h_0$ must be removed from the set of hedges.
- [ ] $h_p$ and $h_n$ must be attached so that $h_n$ comes after $h_p$ in the doubly linked list of the face loop. 
- [ ] $h_n$ must point to $v_1$ (this is not redundant, enforce it).

Note that this leaves the Radial Edge in a broken state, all operations in this section must be applied to restore correctness.

=== Degenerate Face Removal

In almost all cases, the combination of the above operations should be sufficient. However, if any of the faces incident on the edge is triangular, then we end up with degenarte faces with only two edges. These must be fixed by removing unnecessary elements and updating pointers.

#figure(
  grid(
    columns: 2,     // 2 means 2 auto-sized columns
    gutter: 2mm,    // space between columns
      image("images/degenerate_face.svg", width: 70%),
      image("images/fixed_face.svg", width: 70%),
    ),
  caption: [State for a degenerate face set before and after fixing it.],
) <degenerate_face>

@degenerate_face shows an overview of how the degenerate faces should be fixed. The checklist is:

- [ ] Remove the face $f$
- [ ] Remove all half edges associated with $f$.
- [ ] Remove one of the two edges.
- [ ] Grab all half edges orbitting the edge that was removed, insert them into the radial cycle of the surviving edge.
- [ ] Update the surviving edge to point to any of these half edges.
- [ ] Update the vertex cycles around the endpoints of the edges to remove the eliminated edge.
- [ ] Update $v_1$ and $v_2$ to point to the surviving edge.
- [ ] Update all half edges in the joint radial cycle to point to the surviving edge.