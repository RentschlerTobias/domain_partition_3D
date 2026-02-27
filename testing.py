from remesh import mesh_stl
from remesh import remsh_quasi_structured

input_path = './stl_files/Case09_post/LV_inner.stl'
output_path = './test.msh'

remsh_quasi_structured(input_stl=input_path, output_msh=output_path)


mesh_stl(
    path_to_stl="./stl_files/Case09_post/LV_inner.stl",
    output_path="LV_inner.msh",
    element_size=4
)
