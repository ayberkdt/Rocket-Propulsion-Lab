# User guide

## Start and navigate

Run `scripts/run.ps1`, then open <http://127.0.0.1:8765>. The server is local;
no analysis input is sent to an external service. Every form uses SI units and
shows its unit beside the field or result. Use the skip link or the horizontal
module navigation with the keyboard.

## Recommended analysis chain

1. Use **Termokimya** to inspect temperature-dependent species or frozen
   mixture properties.
2. Use **Denge** for a local experimental O/F sweep. Configure
   `RPLAB_CEA_EXECUTABLE` when an independently installed NASA CEA executable is
   available; the repository intentionally does not distribute it.
3. Size the ideal nozzle, then inspect **Tasarım dışı nozul** over ambient
   pressure or altitude. Treat the Summerfield result as a screening warning,
   not a viscous separation solution.
4. Compare conical, preliminary, Rao, and MOC geometry. MOC results include
   boundary residuals and mesh-refinement evidence. Export geometry only after
   reviewing the stated model fidelity.
5. Apply the loss budget and Bartz thermal model. Gas/coolant properties are
   explicit; missing thermal data stops the calculation.
6. Optimize staging, close a preliminary cycle balance, and run seeded
   uncertainty. Save the cases as `.rplab.json` or a printable HTML report.

## Workspaces

The current schema is version 2. Opening a version-1 document migrates it in
memory; saving writes version 2. A workspace stores inputs, scalar results,
model versions, warnings, and dataset identifiers. It does not embed large
thermochemical tables. A save/open round trip is deterministic.

## Fidelity labels

- **Perfect gas / frozen:** constant composition and, where stated, constant
  `gamma`.
- **Local Gibbs experimental:** element-balanced ideal-gas equilibrium with a
  limited species database. It is not a replacement for CEA validation.
- **Empirical screen / preliminary:** a correlation or lumped efficiency with a
  declared application band.
- **Experimental cooling:** bulk single-phase energy balance without wall
  conduction, boiling, cavitation, or turbomachinery maps.

The tool is for education and preliminary design, not flight certification,
hazardous-material procedures, manufacturing release, or human-rated decisions.

