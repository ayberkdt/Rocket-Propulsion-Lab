"""Minimal Python workflow using the public domain functions."""

from rocket_propulsion.compressible import calculate_off_design_nozzle
from rocket_propulsion.geometry import generate_rao_contour
from rocket_propulsion.propulsion import calculate_loss_budget

contour = generate_rao_contour(throat_area_m2=0.01, area_ratio=25.0)
nozzle = calculate_off_design_nozzle(
    chamber_pressure_pa=5_000_000.0,
    chamber_temperature_k=3500.0,
    throat_area_m2=contour.throat_area_m2,
    area_ratio=contour.area_ratio,
    ambient_pressure_pa=101_325.0,
    gamma=1.22,
    gas_constant_j_kg_k=355.0,
)
losses = calculate_loss_budget(
    ideal_thrust_n=nozzle.thrust_n,
    ideal_specific_impulse_s=nozzle.specific_impulse_s,
    divergence_efficiency=contour.estimated_divergence_efficiency or 1.0,
)

print(f"regime={nozzle.regime}")
print(f"ideal thrust={nozzle.thrust_n:.1f} N")
print(f"delivered thrust={losses.delivered_thrust_n:.1f} N")

