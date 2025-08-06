using Gmsh
function generate_mesh1d_mt(zgrid, fname)
    # Initialize gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("1DLine")

    # Your input grid points
    xgrid = zgrid # [0.0, 0.25, 0.5, 0.75, 1.0]

    # Add points to the geometry
    point_tags = [gmsh.model.geo.addPoint(x, 0.0, 0.0, 1e6) for x in xgrid]

    # Create lines between successive points
    line_tags = [
        gmsh.model.geo.addLine(point_tags[i], point_tags[i+1]) for
        i = 1:length(point_tags)-1
    ]

    # Create a line loop and curve physical group (optional)
    gmsh.model.geo.synchronize()

    # Tag the lines as a physical group (e.g., for BCs or domain marking)
    gmsh.model.addPhysicalGroup(1, line_tags, 1)
    gmsh.model.setPhysicalName(1, 1, "LineDomain")

    gmsh.model.addPhysicalGroup(0, [point_tags[1]], 101)  # left endpoint
    gmsh.model.setPhysicalName(0, 101, "tag1")

    gmsh.model.addPhysicalGroup(0, [point_tags[end]], 102)  # right endpoint
    gmsh.model.setPhysicalName(0, 102, "tag2")

    # Generate 1D mesh
    gmsh.model.mesh.generate(1)

    # Optional: write to file
    gmsh.write(fname)

    # Finalize
    gmsh.finalize()
end
