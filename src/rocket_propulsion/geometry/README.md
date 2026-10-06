# Geometry

Parametric geometry that connects dimensions to the compressible-flow core.

`nozzle_contour.py` generates conical or explicitly preliminary Hermite C-D
profiles. `rao.py` supplies an endpoint-tangent Rao/TOP approximation and
`moc.py` constructs a planar characteristic net with a first-order
axisymmetric area mapping. Every model returns the same immutable
`NozzleContour` station contract.

`comparison.py` reports length, surface area, divergence loss, exit angle, and
mesh evidence on a common basis. `export.py` emits CSV, standalone SVG, ASCII
DXF, and closed revolved ASCII STL. These are engineering exchange artifacts,
not manufacturing drawings; MOC model fidelity and refinement residuals remain
part of the result.

