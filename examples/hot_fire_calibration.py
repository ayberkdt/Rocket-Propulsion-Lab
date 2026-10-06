"""Calibrate residual thrust uncertainty from repeated hot-fire scales.

The numbers are synthetic software fixtures, not flight-engine evidence.

References
----------
NIST TN 1297: https://www.nist.gov/pml/nist-technical-note-1297
NIST/SEMATECH variance components:
https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
NIST Dataplot consensus mean:
https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
"""

from __future__ import annotations

from rocket_propulsion.propulsion.burns import (
    CalibrationScope,
    HotFireScaleMeasurement,
    calibrate_hot_fire_scale,
)


def main() -> None:
    """Print correction and residual variance components.

    References
    ----------
    NIST Dataplot consensus mean:
    https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
    """

    offsets = (-0.025, -0.015, -0.005, 0.005, 0.015, 0.025)
    observations = tuple(
        HotFireScaleMeasurement(
            engine_id=f"engine-{engine_index + 1}",
            run_id=f"run-{run_index + 1}",
            condition_id="rated-steady-state",
            value=1.008 + offset + deviation,
            standard_uncertainty=0.001,
            source="synthetic qualification campaign",
        )
        for engine_index, offset in enumerate(offsets)
        for run_index, deviation in enumerate((-0.004, 0.0, 0.004))
    )
    result = calibrate_hot_fire_scale(
        "thrust_scale",
        observations,
        source="synthetic qualification campaign",
        rationale="demonstrate engine/run/measurement variance separation",
    )
    population = result.residual_uncertain_parameter(
        scope=CalibrationScope.POPULATION
    )
    same_engine = result.residual_uncertain_parameter(
        scope=CalibrationScope.SAME_ENGINE
    )
    print(f"consensus correction: {result.consensus_scale:.9f}")
    print(
        "process standard deviations: "
        f"within={result.within_engine_process_variance**0.5:.9f}, "
        f"between={result.between_engine_process_variance**0.5:.9f}"
    )
    print(
        "residual relative sigma: "
        f"population={population.standard_deviation:.9f}, "
        f"same-engine={same_engine.standard_deviation:.9f}"
    )
    print(f"solver iterations: {result.evidence.iterations}")
    print(f"calibration hash: {result.calibration_hash}")
    for warning in result.warnings:
        print(f"warning: {warning}")


if __name__ == "__main__":
    main()

