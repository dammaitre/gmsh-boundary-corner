import sys, os, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "api"))
os.environ["GMSH_LIB"] = os.path.join(os.path.dirname(__file__), "build", "libgmsh.so")

import gmsh

gui = "--gui" in sys.argv
gmsh.initialize(["gmsh"] if gui else ["gmsh", "-nopopup"])
gmsh.model.add("test_bc")

p1 = gmsh.model.geo.addPoint(0, 0, 0, 0.1)
p2 = gmsh.model.geo.addPoint(1, 0, 0, 0.1)
p3 = gmsh.model.geo.addPoint(1, 1, 0, 0.1)
p4 = gmsh.model.geo.addPoint(0, 1, 0, 0.1)
l1 = gmsh.model.geo.addLine(p1, p2)   # profil (y=0)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)
cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
sf = gmsh.model.geo.addPlaneSurface([cl])
gmsh.model.geo.synchronize()

bc = gmsh.model.mesh.field.add("BoundaryCorner")
gmsh.model.mesh.field.setNumbers(bc, "CurvesList",      [l1])
gmsh.model.mesh.field.setNumbers(bc, "AxisPoint",       [1.0, 0.0])
gmsh.model.mesh.field.setNumbers(bc, "StartPoint",      [0.7, 0.0])
gmsh.model.mesh.field.setNumber (bc, "Size",            0.01)
gmsh.model.mesh.field.setNumber (bc, "Ratio",           1.15)
gmsh.model.mesh.field.setNumber (bc, "NbLayers",        5)
gmsh.model.mesh.field.setNumber (bc, "NbCornerColumns", 4)
gmsh.model.mesh.field.setAsBackgroundMesh(bc)

gmsh.option.setNumber("Mesh.BoundaryCornerField", bc)

gmsh.model.mesh.generate(2)

elem_types, elem_tags, _ = gmsh.model.mesh.getElements(2)
print(f"\n=== Résultat ===")
type_names = {2: "triangle", 3: "quadrangle"}
for t, tags in zip(elem_types, elem_tags):
    if len(tags) == 0:
        continue
    name = type_names.get(t, f"type{t}")
    print(f"{name} — n={len(tags)}")

gmsh.write("/tmp/test_bc.msh")
print("\nfichier écrit : /tmp/test_bc.msh")

if gui:
    subprocess.run(["/usr/bin/gmsh", "/tmp/test_bc.msh"])

gmsh.finalize()
