from pyroll.core import RollPass


@RollPass.hookspec
def forward_tension(roll_pass: RollPass):
    """Forward tension applied to the roll pass."""


@RollPass.hookspec
def backward_tension(roll_pass: RollPass):
    """Backward tension applied to the roll pass."""


@RollPass.hookspec
def equivalent_reduction(roll_pass: RollPass):
    """Equivalent reduction of stock."""


@RollPass.hookspec
def neutral_line_angle(roll_pass: RollPass):
    """Equivalent roll angle of the neutral line."""


@RollPass.hookspec
def lippmann_mahrenholz_roll_torque_function(roll_pass: RollPass):
    """Lippmann and Mahrenholz equation for factor Q of roll torque calculation."""


@RollPass.hookspec
def inverse_forming_efficiency(roll_pass: RollPass):
    """Lippmann and Mahrenholz equation for inverse forming efficiency."""
