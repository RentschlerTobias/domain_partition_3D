Regression goldens: one small JSON of reference numbers per case (block count, inverted cells, min/p1
scaled Jacobian, Hausdorff), never meshes. Filled in step 1.4 as the modules migrate. Inputs that are
not versioned (dtOO meshes, HexMe) are looked up through HEXDP_DATA and the tests skip when it is unset.
