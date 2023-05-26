import logging
import webbrowser
from pathlib import Path

from pyroll.core import Profile, PassSequence, RollPass, Roll, CircularOvalGroove, Transport, RoundGroove


def test_solve(tmp_path: Path, caplog):
    caplog.set_level(logging.INFO, logger="pyroll")

    def test_solve(tmp_path: Path, caplog):
        caplog.set_level(logging.INFO, logger="pyroll")

        import pyroll.lippmann_mahrenholz_power_and_labour

        in_profile = Profile.round(
            diameter=30e-3,
            temperature=1200 + 273.15,
            strain=0,
            material=["C45", "steel"],
        )

        sequence = PassSequence([
            RollPass(
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
                front_tension=6e6,

            ),
            Transport(
                label="I => II",
                duration=1
            ),
            RollPass(
                label="Round II",
                roll=Roll(
                    groove=RoundGroove(
                        r1=1e-3,
                        r2=12.5e-3,
                        depth=11.5e-3
                    ),
                    nominal_radius=160e-3,
                    rotational_frequency=1
                ),
                gap=2e-3,
                back_tension=6e6,
                front_tension=0,

            ),
        ])

        try:
            sequence.solve(in_profile)
        finally:
            print("\nLog:")
            print(caplog.text)

        try:
            import pyroll.report

            report = pyroll.report.report(sequence)

            report_file = tmp_path / "report.html"
            report_file.write_text(report)
            print(report_file)
            webbrowser.open(report_file.as_uri())

        except ImportError:
            pass
