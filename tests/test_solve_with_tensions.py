import logging
import webbrowser
from pathlib import Path

from pyroll.core import Profile, PassSequence, RollPass, Roll, CircularOvalGroove, Transport, RoundGroove


def test_solve(tmp_path: Path, caplog):
    caplog.set_level(logging.INFO, logger="pyroll")

    import pyroll.lippmann_mahrenholz_force_torque

    in_profile = Profile.round(
        diameter=30e-3,
        temperature=1200 + 273.15,
        strain=0,
        material=["C45", "steel"],
        flow_stress=100e6,
        density=7.5e3,
        specific_heat_capacity=690,
    )

    roll_pass_with_front_tension = RollPass(
        label="Oval I",
        roll=Roll(
            groove=CircularOvalGroove(
                depth=8e-3,
                r1=6e-3,
                r2=40e-3
            ),
            nominal_radius=160e-3,
            rotational_frequency=1
        ),
        gap=2e-3,
        back_tension=0,
        front_tension=10e6
    )

    roll_pass_with_back_tension = RollPass(
        label="Oval I",
        roll=Roll(
            groove=CircularOvalGroove(
                depth=8e-3,
                r1=6e-3,
                r2=40e-3
            ),
            nominal_radius=160e-3,
            rotational_frequency=1
        ),
        gap=2e-3,
        back_tension=10e6,
        front_tension=0
    )

    roll_pass_without_tension = RollPass(
        label="Oval I",
        roll=Roll(
            groove=CircularOvalGroove(
                depth=8e-3,
                r1=6e-3,
                r2=40e-3
            ),
            nominal_radius=160e-3,
            rotational_frequency=1
        ),
        gap=2e-3,
        back_tension=0,
        front_tension=0
    )

    roll_pass_with_back_tension.solve(in_profile)
    roll_pass_with_front_tension.solve(in_profile)
    roll_pass_without_tension.solve(in_profile)

    assert roll_pass_with_back_tension.roll.roll_torque > roll_pass_without_tension.roll.roll_torque
    assert roll_pass_without_tension.roll.roll_torque > roll_pass_with_front_tension.roll.roll_torque
    assert roll_pass_with_back_tension.entry_point < roll_pass_with_back_tension.roll.neutral_point < roll_pass_without_tension.exit_point

    assert roll_pass_with_back_tension.roll.neutral_point > roll_pass_without_tension.roll.neutral_point > roll_pass_with_front_tension.roll.neutral_point
