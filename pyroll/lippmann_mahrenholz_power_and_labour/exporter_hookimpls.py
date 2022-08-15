from pyroll.core import RollPass
from pyroll.ui import Exporter
from pyroll.utils import for_units


@Exporter.hookimpl
@for_units(RollPass)
def columns(unit: RollPass):
    return dict(
        inverse_efficiency=f"{unit.inverse_forming_efficiency:.2f}",
        equivalent_neutral_line_angle=f"{unit.equivalent_neutral_line_angle:.2f}"
    )
