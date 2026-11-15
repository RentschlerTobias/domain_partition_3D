"""hexdp: hexahedral block-structure generation and domain partition.

Subpackages, in pipeline order:

    io          read/write msh, vtk, ovm, json
    core        Geometry (patch roles, feature curves), BlockStructure, lattice
    detect      detectors: rotation axis / sweep test, thin walls, blades, periodicity
    strategies  block-structure strategies: algohex, sweep, cgrid, blade_ring
    ops         block operations: collapse, ogrid/pillow, split
    discretise  cell counts per edge, grading, first cell height from y+
    fill        TFI, Winslow smoothing
    quality     scaled Jacobian, non-orthogonality, Hausdorff, volume ratios
    export      vtk, msh, OpenFOAM polyMesh, CAD
    interfaces  conformal / cyclicAMI / cyclicPeriodicAMI / mixingPlane
    pipeline    case YAML, runner, jobs
    mcp         MCP server exposing the tools to agents

Plugins of every kind are found through entry points, see `hexdp.registry`.
"""

__version__ = "0.1.0"
