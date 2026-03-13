# Mesh Splitting Algorithm

This crate implements a method to split the meshes when the intersecting faces
are coplanar. This allows the geometry to match at the interfaces of the meshes,
which enables ray tracing for translucent materials reliably.


## Implementation

Features such as vertices, their normals and faces (polygons) generate from
these are kept in a list and indexed by their position in the list as shown
below. Hence,
matching vertices is done based on an identifier rather than comparing floating
points. This mitigates problems in operation with floating point errors,
avoiding openings in the mesh and computationally lighter to computer matches.

![Features global indexing for
construction](./imgs/Meshes_Set_Indexing.excalidraw.png)


- [Mesh splitting algorithm explained in
depth](https://git-pages.ecdf.ed.ac.uk/lidar-research/models/aetherus/simulation-nb/website/aetherus/remesh.html)

