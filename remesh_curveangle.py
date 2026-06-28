import subprocess
import sys
import os

OUT_DIR = "/root/repos/block_structured_meshing"

def run_gmsh_remesh(stl_file, name, out_msh, elem_size=1.0, timeout=30):
    script = f"""
import gmsh
import sys
import numpy as np

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.option.setNumber("Mesh.Algorithm", 11)

gmsh.merge("{stl_file}")

# classify surfaces with larger curve angle to merge small curves
gmsh.model.mesh.classifySurfaces(np.pi/2, True, True, np.pi/2)

print("After classify:")
print(f"  dim1: {{len(gmsh.model.getEntities(1))}}")
print(f"  dim2: {{len(gmsh.model.getEntities(2))}}")

# create geometry
gmsh.model.mesh.createGeometry()

print("After createGeometry:")
print(f"  dim1: {{len(gmsh.model.getEntities(1))}}")
print(f"  dim2: {{len(gmsh.model.getEntities(2))}}")

# set background field
field = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(field, "F", "{elem_size}")
gmsh.model.mesh.field.setAsBackgroundMesh(field)

# options
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", {elem_size})
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", {elem_size})
gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)
gmsh.option.setNumber("Mesh.QuadqsSizemapMethod", 0)

# generate mesh
try:
    gmsh.model.mesh.generate(2)
    
    # info
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2)
    n_quads = 0
    n_tris = 0
    for etype, etags, enodes in zip(elem_types, elem_tags, elem_node_tags):
        if etype == 3:
            n_quads = len(etags)
        elif etype == 2:
            n_tris = len(etags)
    
    print(f"RESULT: {{n_quads}} quads, {{n_tris}} tris")
    gmsh.write("{out_msh}")
    print(f"SAVED: {out_msh}")
except Exception as e:
    print(f"ERROR: {{e}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

gmsh.finalize()
"""
    
    tmp_script = f"{OUT_DIR}/tmp_remesh_{name}.py"
    with open(tmp_script, 'w') as f:
        f.write(script)
    
    env = os.environ.copy()
    env['PYTHONPATH'] = '/root/venv/lib/python3.12/site-packages'
    
    try:
        result = subprocess.run(
            ['/root/venv/bin/python3', tmp_script],
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env
        )
        print(f"[{name}] return code: {result.returncode}")
        print(f"[{name}] stdout:\n{result.stdout}")
        if result.stderr:
            print(f"[{name}] stderr:\n{result.stderr}")
    except subprocess.TimeoutExpired:
        print(f"[{name}] TIMEOUT after {timeout}s")
    finally:
        if os.path.exists(tmp_script):
            os.remove(tmp_script)

if __name__ == "__main__":
    run_gmsh_remesh(f"{OUT_DIR}/T1_9_hub_raw.stl", "hub", f"{OUT_DIR}/T1_9_hub_quad.msh", elem_size=1.0, timeout=30)
    run_gmsh_remesh(f"{OUT_DIR}/T1_9_shroud_raw.stl", "shroud", f"{OUT_DIR}/T1_9_shroud_quad.msh", elem_size=1.0, timeout=30)
