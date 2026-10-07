"""Thruster torque about a migrating centre of mass.

A four-thruster aft cluster fires while oxidizer drains; the centre of mass
moves and the torque from one stronger thruster changes with it.
Numbers are illustrative.
"""

from rocket_propulsion.propulsion.burns import (
    PropellantMassModel,
    PropellantStation,
    ThrusterLayout,
    ThrusterMount,
)

layout = ThrusterLayout(
    tuple(
        ThrusterMount(name, position, (0.0, 0.0, 1.0))
        for name, position in (
            ("T1", (0.4, 0.0, -1.2)),
            ("T2", (-0.4, 0.0, -1.2)),
            ("T3", (0.0, 0.4, -1.2)),
            ("T4", (0.0, -0.4, -1.2)),
        )
    )
)
mass_model = PropellantMassModel(
    dry_mass_kg=800.0,
    dry_center_of_mass_m=(0.0, 0.0, 0.4),
    dry_inertia_kg_m2=((400.0, 0.0, 0.0), (0.0, 420.0, 0.0), (0.0, 0.0, 300.0)),
    stations=(
        PropellantStation("oxidizer_mass_kg", (0.05, 0.0, -0.2)),
        PropellantStation("fuel_mass_kg", (-0.05, 0.0, 0.1)),
    ),
)
thrusts = {"T1": 22.0, "T2": 20.0, "T3": 20.0, "T4": 20.0}

for oxidizer in (300.0, 150.0, 20.0):
    properties = mass_model.evaluate({"oxidizer_mass_kg": oxidizer, "fuel_mass_kg": 120.0})
    result = layout.resultant(thrusts, properties.center_of_mass_m)
    print(
        f"oxidizer {oxidizer:6.1f} kg  CoM z = {properties.center_of_mass_m[2]: .4f} m  "
        f"torque = {tuple(round(value, 4) for value in result.torque_n_m)} N m"
    )
print("axial efficiency:", layout.axial_efficiency(thrusts, (0.0, 0.0, 1.0)))

# Allocate thrust plus a pitch torque within 0..25 N per thruster (thrusters only push).
allocation = layout.allocate(
    force_n=(0.0, 0.0, 60.0),
    torque_n_m=(0.0, 1.5, 0.0),
    center_of_mass_m=properties.center_of_mass_m,
    maximum_thrust_n=25.0,
)
print("allocation:", {name: round(value, 3) for name, value in allocation.thrusts_n})
print("exact:", allocation.exact, "at maximum:", allocation.thrusters_at_maximum)
