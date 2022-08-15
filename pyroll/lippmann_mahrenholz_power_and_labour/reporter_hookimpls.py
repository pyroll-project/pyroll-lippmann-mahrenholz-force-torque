from pyroll.core import RollPass
from pyroll.ui import Reporter
from pyroll.utils import for_units


@Reporter.hookimpl
@for_units(RollPass)
def unit_properties(unit: RollPass):
    return dict(
        inverse_efficiency=f"{unit.inverse_forming_efficiency:.2f}",
        neutral_line_angle=f"{unit.neutral_line_angle:.2f}"
    )
