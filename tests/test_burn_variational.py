"""Variational equations built from the analytic partials match finite differences.

A consumer that estimates propulsion drivers integrates the sensitivity matrix
``S = dy/dp`` alongside the state:

    dS/dt = (df/dy) S + df/dp

and applies the saltation jump ``S+ = S- + (f- - f+) dtau/dp`` at every event
whose time ``tau`` depends on a parameter (here: every profile time event
moves with ``ignition_time_bias_s``).  These tests integrate exactly that with
the bridge's and the electric model's partials and compare the final
sensitivities with central finite differences of whole propagations.
"""

import unittest

from rocket_propulsion.propulsion.burns import (
    FixedEfficiencyThrottle,
    PowerLimitedPropulsionModel,
    PropagatorPropulsionBridge,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
)

PARAMETERS = ("thrust_scale", "mass_flow_scale", "ignition_time_bias_s")


def jumpy_artifact() -> TabulatedBurnArtifact:
    """Start jump, ramp, mid-burn throttle jump and stop jump."""

    def streams(flow: float) -> tuple[tuple[str, float], ...]:
        return (("ox", 0.7 * flow), ("fu", 0.3 * flow))

    return TabulatedBurnArtifact(
        segments=(
            TabulatedBurnSegment(
                0.0, 2.0, 80.0, 160.0, 80.0, 160.0, 0.40, 0.70, streams(0.40), streams(0.70), "ramp"
            ),
            TabulatedBurnSegment(
                2.0,
                5.0,
                120.0,
                120.0,
                120.0,
                120.0,
                0.55,
                0.55,
                streams(0.55),
                streams(0.55),
                "throttled",
            ),
        ),
        source_id="variational",
        source_sha256="5" * 64,
        reference_ids=("NASA-20040000363",),
    )


class BridgeVariationalTests(unittest.TestCase):
    def setUp(self) -> None:
        self.bridge = PropagatorPropulsionBridge(
            jumpy_artifact(),
            stream_bindings=(StreamStateBinding("ox", "m_ox"), StreamStateBinding("fu", "m_fu")),
        )

    def rhs(self, t, y, p, side):
        value = self.bridge.evaluate_at_epoch(
            time_s=t,
            ignition_time_s=0.0,
            vehicle_mass_kg=y[1],
            direction=(1.0, 0.0, 0.0),
            additional_states={"m_ox": y[2], "m_fu": y[3]},
            parameter_overrides=dict(zip(PARAMETERS, p, strict=True)),
            side=side,
        )
        drains = dict(value.additional_state_derivatives)
        f = [value.acceleration_m_s2[0], value.mass_derivative_kg_s, drains["m_ox"], drains["m_fu"]]
        partials = value.partials
        dfdy = [[0.0] * 4 for _ in range(4)]
        dfdy[0][1] = partials.acceleration_wrt_mass[0]
        acc_p = dict(partials.acceleration_wrt_parameters)
        mass_p = dict(partials.mass_rate_wrt_parameters)
        state_p = {
            name: dict(rates) for name, rates in partials.additional_state_rates_wrt_parameters
        }
        dfdp = [
            [acc_p[name][0] for name in PARAMETERS],
            [mass_p[name] for name in PARAMETERS],
            [state_p[name]["m_ox"] for name in PARAMETERS],
            [state_p[name]["m_fu"] for name in PARAMETERS],
        ]
        return f, dfdy, dfdp

    def propagate(self, p, *, with_sensitivity: bool, step: float = 0.01):
        t, end = -0.5, 6.5
        y = [0.0, 120.0, 30.0, 12.0]
        s = [[0.0] * 3 for _ in range(4)]
        bias = p[2]
        events = sorted(
            e.relative_time_s + bias
            for e in self.bridge.event_surfaces()
            if e.relative_time_s is not None
        )

        def combined(tt, yy, ss, side):
            f, dfdy, dfdp = self.rhs(tt, yy, p, side)
            ds = [
                [sum(dfdy[i][k] * ss[k][j] for k in range(4)) + dfdp[i][j] for j in range(3)]
                for i in range(4)
            ]
            return f, ds

        while t < end:
            target = min(end, t + step)
            for event in events:
                if t < event < target:
                    target = event
                    break
            h = target - t

            def add(base, k, factor):
                return [a + factor * b for a, b in zip(base, k, strict=True)]

            def add_s(base, k, factor):
                return [add(row, krow, factor) for row, krow in zip(base, k, strict=True)]

            f1, s1 = combined(t, y, s, "right")
            f2, s2 = combined(t + h / 2, add(y, f1, h / 2), add_s(s, s1, h / 2), "right")
            f3, s3 = combined(t + h / 2, add(y, f2, h / 2), add_s(s, s2, h / 2), "right")
            f4, s4 = combined(target, add(y, f3, h), add_s(s, s3, h), "left")
            y = [
                a + h / 6 * (b + 2 * c + 2 * d + e)
                for a, b, c, d, e in zip(y, f1, f2, f3, f4, strict=True)
            ]
            s = [
                [a + h / 6 * (b + 2 * c + 2 * d + e) for a, b, c, d, e in zip(*rows, strict=True)]
                for rows in zip(s, s1, s2, s3, s4, strict=True)
            ]
            t = target
            if with_sensitivity and t in events:
                # Saltation: every profile event time moves one-for-one with the bias.
                jump = self.bridge.time_event_jump_at_epoch(
                    time_s=t,
                    ignition_time_s=0.0,
                    vehicle_mass_kg=y[1],
                    direction=(1.0, 0.0, 0.0),
                    additional_states={"m_ox": y[2], "m_fu": y[3]},
                    parameter_overrides=dict(zip(PARAMETERS, p, strict=True)),
                )
                rates = dict(jump.additional_state_derivatives)
                for i, value in enumerate(
                    (
                        jump.acceleration_m_s2[0],
                        jump.mass_derivative_kg_s,
                        rates["m_ox"],
                        rates["m_fu"],
                    )
                ):
                    s[i][2] += value
        return y, s

    def test_sensitivities_match_finite_differences(self) -> None:
        nominal = (1.0, 1.0, 0.0)
        _, sensitivity = self.propagate(nominal, with_sensitivity=True)
        steps = (1e-5, 1e-5, 1e-5)
        for column, h in enumerate(steps):
            plus = list(nominal)
            minus = list(nominal)
            plus[column] += h
            minus[column] -= h
            y_plus, _ = self.propagate(tuple(plus), with_sensitivity=False)
            y_minus, _ = self.propagate(tuple(minus), with_sensitivity=False)
            for row in range(4):
                finite = (y_plus[row] - y_minus[row]) / (2 * h)
                scale = max(1.0, abs(finite))
                self.assertAlmostEqual(
                    sensitivity[row][column] / scale,
                    finite / scale,
                    places=6,
                    msg=f"d{row}/d{PARAMETERS[column]}",
                )

    def test_bias_sensitivity_is_wrong_without_saltation(self) -> None:
        """Guards the test itself: dropping the jump term must be detectable."""

        nominal = (1.0, 1.0, 0.0)
        _, with_jump = self.propagate(nominal, with_sensitivity=True)
        _, without_jump = self.propagate(nominal, with_sensitivity=False)
        self.assertGreater(abs(with_jump[0][2] - without_jump[0][2]), 1e-3)

    def test_free_space_clock_translation_preserves_state_and_zero_bias_sensitivity(self) -> None:
        """Time-homogeneous free flight is invariant to shifting the complete burn."""
        reference, sensitivity = self.propagate((1.0, 1.0, 0.0), with_sensitivity=True)
        # Independent exact trapezoid ledger: 1.1 + 1.65 kg, with 70/30 streams.
        for actual, expected in zip(reference[1:], (117.25, 28.075, 11.175), strict=True):
            self.assertLess(abs(actual - expected), 1e-10)
        for row in sensitivity:
            self.assertLess(abs(row[2]), 1e-9)
        for bias in (-0.125, 0.125):
            shifted, _ = self.propagate((1.0, 1.0, bias), with_sensitivity=False)
            for actual, expected in zip(shifted, reference, strict=True):
                self.assertLess(abs(actual - expected), 1e-9)


class ElectricVariationalTests(unittest.TestCase):
    """Chain rule through the consumer's power model: P(t) = P0 + k t."""

    def test_power_and_scale_sensitivities_match_finite_differences(self) -> None:
        model = PowerLimitedPropulsionModel(FixedEfficiencyThrottle(0.6, 2500.0, 500.0, 5000.0))

        def propagate(p0, thrust_scale, with_sensitivity):
            y = [0.0, 400.0, 30.0]  # v, m, xenon
            s = [[0.0, 0.0] for _ in range(3)]  # d/dP0, d/dthrust_scale
            t, end, h = 0.0, 3000.0, 25.0
            while t < end - 1e-9:

                def rhs(tt, yy):
                    value = model.evaluate(
                        level_index=0,
                        available_power_w=p0 - 0.5 * tt,
                        vehicle_mass_kg=yy[1],
                        direction=(1.0, 0.0, 0.0),
                        additional_states={"xenon_mass_kg": yy[2]},
                        parameter_overrides={"thrust_scale": thrust_scale},
                    )
                    partials = value.partials
                    f = [
                        value.acceleration_m_s2[0],
                        value.mass_derivative_kg_s,
                        value.mass_derivative_kg_s,
                    ]
                    dfdy = [[0.0, partials.acceleration_wrt_mass[0], 0.0], [0.0] * 3, [0.0] * 3]
                    acc = dict(partials.acceleration_wrt_parameters)
                    flow = dict(partials.mass_rate_wrt_parameters)
                    dfdp = [
                        [partials.acceleration_wrt_available_power[0], acc["thrust_scale"][0]],
                        [partials.mass_rate_wrt_available_power, flow["thrust_scale"]],
                        [partials.mass_rate_wrt_available_power, flow["thrust_scale"]],
                    ]
                    return f, dfdy, dfdp

                def step(tt, yy, ss):
                    f, dfdy, dfdp = rhs(tt, yy)
                    ds = [
                        [
                            sum(dfdy[i][k] * ss[k][j] for k in range(3)) + dfdp[i][j]
                            for j in range(2)
                        ]
                        for i in range(3)
                    ]
                    return f, ds

                def add(base, k, factor):
                    return [a + factor * b for a, b in zip(base, k, strict=True)]

                f1, s1 = step(t, y, s)
                f2, s2 = step(
                    t + h / 2,
                    add(y, f1, h / 2),
                    [add(r, k, h / 2) for r, k in zip(s, s1, strict=True)],
                )
                f3, s3 = step(
                    t + h / 2,
                    add(y, f2, h / 2),
                    [add(r, k, h / 2) for r, k in zip(s, s2, strict=True)],
                )
                f4, s4 = step(
                    t + h, add(y, f3, h), [add(r, k, h) for r, k in zip(s, s3, strict=True)]
                )
                y = [
                    a + h / 6 * (b + 2 * c + 2 * d + e)
                    for a, b, c, d, e in zip(y, f1, f2, f3, f4, strict=True)
                ]
                s = [
                    [
                        a + h / 6 * (b + 2 * c + 2 * d + e)
                        for a, b, c, d, e in zip(*rows, strict=True)
                    ]
                    for rows in zip(s, s1, s2, s3, s4, strict=True)
                ]
                t += h
            return y, s

        _, sensitivity = propagate(3000.0, 1.0, True)
        for column, (h, args) in enumerate(
            ((1.0, lambda d: (3000.0 + d, 1.0)), (1e-5, lambda d: (3000.0, 1.0 + d)))
        ):
            plus, _ = propagate(*args(h), False)
            minus, _ = propagate(*args(-h), False)
            for row in range(3):
                finite = (plus[row] - minus[row]) / (2 * h)
                scale = max(abs(finite), 1e-12)
                self.assertAlmostEqual(sensitivity[row][column] / scale, finite / scale, places=5)


if __name__ == "__main__":
    unittest.main()
