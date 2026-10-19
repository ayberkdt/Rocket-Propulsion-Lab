"""Precision-controlled rocket-equation closure without binary intermediate overflow.

References
----------
NASA-ROCKET: sequential rocket-equation mass ratios; rearrangement derived here.
NASA-7009B: numerical verification, declared limits and uncertainty reporting.
"""

from decimal import Context, Decimal, DecimalException, localcontext
from math import fsum, isfinite
from sys import float_info

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import DomainError


def _loss_fraction(exponent: Decimal) -> Decimal:
    """Evaluate 1-exp(-s), retaining tiny positive s with an alternating series.

    References
    ----------
    NASA-ROCKET: exponential mass ratio; Taylor rearrangement derived here.
    """
    if exponent == 0:
        return Decimal(0)
    if exponent >= Decimal("0.5"):
        return 1 - (-exponent).exp()
    term = total = exponent
    index = 2
    while True:
        term *= -exponent / index
        updated = total + term
        if updated == total:
            return total
        total = updated
        index += 1


def _binary_rows(rows: tuple) -> tuple[tuple, Decimal]:
    """Convert only final channels; reject positive quantities rounded to zero.

    Closure residuals may round below the binary range; their Decimal strings
    remain in the numerical evidence. Conversion differences describe numeric
    quantization, not physical uncertainty.

    References
    ----------
    NASA-7009B: finite-precision numerical errors and reporting boundaries.
    """
    converted = []
    maximum_difference = Decimal(0)
    for name, value, unit in rows:
        binary = float(value)
        diagnostic = name in {"mass_closure_residual", "propellant_closure_residual"}
        if not isfinite(binary) or (value != 0 and binary == 0 and not diagnostic):
            raise DomainError(f"Mission {name} is outside the representable binary output domain.")
        if value != 0 and not diagnostic:
            maximum_difference = max(maximum_difference, abs(Decimal(binary) - value) / abs(value))
        converted.append((name, binary, unit))
    return tuple(converted), maximum_difference


def _at_precision(dry: float, residual: float, reserve: float,
                  burns: tuple[tuple[float, float, float, float, float], ...],
                  precision: int) -> tuple | None:
    """Evaluate a complete ledger; None requests more precision near feasibility.

    With B=exp(-sum(dv/(g0*Isp*eta))), P=[dry*(1-B)+reserve]/(B-f).
    This inverse-ratio form avoids an overflowing exp(+s) intermediate.
    A deliberately conservative decimal guard reserves ten significant digits
    for accumulated rounding. This is an implementation convergence policy,
    not a rigorous interval proof or a NASA-prescribed acceptance tolerance.

    References
    ----------
    NASA-ROCKET: mass ratio and ordered burns; closure algebra derived here.
    NASA-ISP: mass flow, on-time and impulse identities.
    NASA-7009B: numerical domain and precision/convergence evidence.
    """
    with localcontext(Context(prec=precision)):
        d, f, r = map(Decimal, (dry, residual, reserve))
        g0 = Decimal(STANDARD_GRAVITY_M_S2)
        inputs = tuple(tuple(Decimal(x) for x in burn) for burn in burns)
        exponents = tuple(dv / (g0 * isp * eta) for dv, _, isp, eta, _ in inputs)
        exponent = sum(exponents, Decimal(0))
        # Even with no retained propellant, final mass is at least dry mass.
        if exponent > Decimal(float_info.max).ln() - d.ln():
            raise DomainError("Mission wet mass exceeds the finite binary output domain.")
        survival = (-exponent).exp()
        denominator = survival - f
        roundoff_scale = Decimal(10) ** (10 - precision) * max(Decimal(1), exponent) * len(burns)
        if abs(denominator) <= (survival + f) * roundoff_scale:
            return None
        if denominator < 0:
            raise DomainError("No finite wet mass closes the requested residual and delta-v.")
        loaded = (d * _loss_fraction(exponent) + r) / denominator
        retained = f * loaded + r
        mass = initial = d + loaded
        total_used = total_impulse = total_on = total_elapsed = Decimal(0)
        burn_rows = []
        conversion_difference = Decimal(0)
        for (dv, force, isp, eta, duty), burn_exponent in zip(inputs, exponents, strict=True):
            exhaust = g0 * isp
            used = mass * _loss_fraction(burn_exponent)
            final_mass = mass * (-burn_exponent).exp()
            flow = force / exhaust
            on_time = used / flow
            elapsed = on_time / duty
            impulse = used * exhaust
            rows = (("initial_mass", mass, "kg"), ("final_mass", final_mass, "kg"),
                    ("consumed_propellant", used, "kg"), ("tank_mass_flow", flow, "kg/s"),
                    ("on_time", on_time, "s"), ("elapsed_time", elapsed, "s"),
                    ("scalar_impulse", impulse, "N*s"), ("axial_impulse", impulse * eta, "N*s"),
                    ("initial_axial_acceleration", force * eta / mass, "m/s^2"))
            converted, difference = _binary_rows(rows)
            conversion_difference = max(conversion_difference, difference)
            burn_rows.append(converted)
            mass = final_mass
            total_used += used
            total_impulse += impulse
            total_on += on_time
            total_elapsed += elapsed
        mass_error = mass - d - retained
        propellant_error = loaded - retained - total_used
        if (abs(mass_error) > max(mass, d, retained) * roundoff_scale
                or abs(propellant_error) > max(loaded, retained, total_used) * roundoff_scale):
            return None
        rows = (
            ("dry_mass", d, "kg"), ("wet_mass", initial, "kg"),
            ("loaded_propellant", loaded, "kg"), ("consumed_propellant", total_used, "kg"),
            ("retained_propellant", retained, "kg"), ("final_mass", mass, "kg"),
            ("total_delta_v", sum((b[0] for b in inputs), Decimal(0)), "m/s"),
            ("scalar_impulse", total_impulse, "N*s"), ("on_time", total_on, "s"),
            ("elapsed_time", total_elapsed, "s"),
            ("mass_closure_residual", mass_error, "kg"),
            ("propellant_closure_residual", propellant_error, "kg"),
        )
        converted, difference = _binary_rows(rows)
        binary = {name: value for name, value, _ in converted}
        binary_residuals = {
            "mass_closure_residual": fsum((binary["final_mass"], -binary["dry_mass"],
                                            -binary["retained_propellant"])),
            "propellant_closure_residual": fsum((binary["loaded_propellant"],
                                                  -binary["retained_propellant"],
                                                  -binary["consumed_propellant"])),
        }
        converted = tuple((name, binary_residuals.get(name, value), unit)
                          for name, value, unit in converted)
        diagnostics = (
            ("method", "inverse_survival_ratio_and_decimal_precision_agreement"),
            ("working_precision_digits", str(precision)),
            ("closure_denominator_decimal", str(denominator)),
            ("denominator_subtraction_condition_decimal", str((survival + f) / denominator)),
            ("mass_closure_residual_decimal_kg", str(mass_error)),
            ("propellant_closure_residual_decimal_kg", str(propellant_error)),
            ("maximum_binary_conversion_relative_difference_decimal",
             str(max(conversion_difference, difference))),
            ("interpretation", "numerical diagnostics only; physical uncertainty is not quantified"),
        )
        return converted, tuple(burn_rows), diagnostics


def solve_mass_budget(dry: float, residual: float, reserve: float,
                      burns: tuple[tuple[float, float, float, float, float], ...]) -> tuple:
    """Require identical physical binary channels at two resolved precisions.

    Inputs have already passed the public, evidence-bearing domain checks.
    Decimal precision is isolated from the caller's process-wide context.
    Failure to resolve a boundary by 640 digits is explicit, not feasibility.

    References
    ----------
    NASA-ROCKET: mass closure; NASA-7009B: numerical convergence and limits.
    """
    previous = None
    checked = []
    for precision in (80, 160, 320, 640):
        checked.append(str(precision))
        try:
            current = _at_precision(dry, residual, reserve, burns, precision)
        except DecimalException as error:
            raise DomainError("Mission arithmetic exceeds the declared decimal calculation domain.") from error
        if current is None:
            previous = None
            continue
        # Decimal closure residuals are bounded above in each evaluation. They
        # measure arithmetic roundoff and shrink with precision, so requiring
        # their binary representations to agree would prevent convergence of
        # an otherwise stable physical ledger.
        physical_summary = tuple(row for row in current[0] if not row[0].endswith("closure_residual"))
        signature = physical_summary, current[1]
        if previous is not None and signature == previous:
            return (*current[:2], (*current[2], ("checked_precision_digits", ",".join(checked))))
        previous = signature
    raise DomainError("Mission feasibility/output rounding did not resolve within 640 decimal digits.")
