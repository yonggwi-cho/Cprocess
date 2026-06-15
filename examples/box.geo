// Gmsh example geometry for cprocess: a 1 x 1 x 1.5 um silicon block,
// drawn in micrometres (load with: mesh gmsh file=box.msh scale=1um).
//
//   gmsh -3 -format msh41 box.geo -o box.msh
//
// Note: with OpenCASCADE the surface numbering below matches gmsh 4.x for
// Box(); if your version differs, check the surface ids in the GUI.
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1.5};
MeshSize{ PointsOf{ Volume{1}; } } = 0.08;

Physical Volume("sub") = {1};
Physical Surface("top") = {6};     // z = 1.5 face
Physical Surface("bottom") = {5};  // z = 0 face
