from pyroll.core import RollPass
from pyroll.ui import Exporter, Reporter

from . import roll_pass_hookspecs
from . import roll_pass_hookimpls
from . import exporter_hookimpls, reporter_hookimpls

RollPass.plugin_manager.add_hookspecs(roll_pass_hookspecs)
RollPass.plugin_manager.register(roll_pass_hookimpls)

RollPass.Roll.plugin_manager.register(roll_pass_hookimpls)

Exporter.plugin_manager.register(exporter_hookimpls)
Reporter.plugin_manager.register(reporter_hookimpls)
