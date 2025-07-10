import logging
import numpy as np
from pyroll.core import SymmetricRollPass, Hook

VERSION = "3.0.0"

log = logging.getLogger(__name__)

SymmetricRollPass.inverse_forming_efficiency = Hook[float]()
"""Inverse forming efficiency defined by Lippmann and Mahrenholz for the roll force of the roll pass."""

SymmetricRollPass.roll_torque_loss_function = Hook[float]()
"""Loss function defined by Lippmann and Mahrenholz for the roll torque of the roll pass."""

SymmetricRollPass.Roll.relative_neutral_angle = Hook[float]()
"""Relative Neutral angle defined by Lippmann and Mahrenholz for the roll pass."""

@SymmetricRollPass.Roll.relative_neutral_angle
def relative_neutral_angle(self: SymmetricRollPass.Roll):
    rp = self.roll_pass
    mean_flow_stress = (rp.in_profile.flow_stress + 2 * rp.out_profile.flow_stress) / 3
    abs_rel_draught  = abs(rp.rel_draught)
    back_tension_model_definition = -rp.back_tension
    front_tension_model_definition = -rp.front_tension
    relative_tension = (back_tension_model_definition - front_tension_model_definition) / mean_flow_stress

    p1 = np.sqrt((1 - abs_rel_draught) / abs_rel_draught)
    p2 = 1 / 2 * np.sqrt(rp.out_profile.equivalent_height / self.working_radius) * (relative_tension + np.log(1 - abs_rel_draught))
    p3 = 1 / 2 * np.arctan(np.sqrt(abs_rel_draught / (1 - abs_rel_draught)))

    return p1 * np.tan(p2 + p3)

@SymmetricRollPass.Roll.neutral_angle
def neutral_angle(self: SymmetricRollPass.Roll):
    return self.entry_angle * self.relative_neutral_angle

@SymmetricRollPass.inverse_forming_efficiency
def inverse_forming_efficiency(self: SymmetricRollPass):
    if np.isclose(self.rel_draught, 0):
        return 1

    mean_flow_stress = (self.in_profile.flow_stress + 2 * self.out_profile.flow_stress) / 3
    abs_rel_draught = abs(self.rel_draught)
    back_tension_model_definition = - self.back_tension

    return back_tension_model_definition/ mean_flow_stress + 2 * np.sqrt(
        (1 - abs_rel_draught) / abs_rel_draught) * np.arctan(
        np.sqrt(abs_rel_draught / (1 - abs_rel_draught))) - 1 + np.sqrt(
        self.roll.working_radius / self.out_profile.equivalent_height) * np.sqrt(
        (1 - abs_rel_draught) / abs_rel_draught) * np.log(
        np.sqrt(1 - abs_rel_draught) / (1 - abs_rel_draught * (1 - self.roll.relative_neutral_angle ** 2)))


@SymmetricRollPass.DiskElement.deformation_resistance
def deformation_resistance(self: SymmetricRollPass.DiskElement):
    mean_flow_stress = (self.in_profile.flow_stress + 2 * self.out_profile.flow_stress) / 3
    return mean_flow_stress * self.roll_pass.inverse_forming_efficiency


@SymmetricRollPass.deformation_resistance
def deformation_resistance(self: SymmetricRollPass):
    mean_flow_stress = (self.in_profile.flow_stress + 2 * self.out_profile.flow_stress) / 3
    return mean_flow_stress * self.inverse_forming_efficiency


@SymmetricRollPass.roll_force
def roll_force(self: SymmetricRollPass):
    return self.roll.contact_area * self.deformation_resistance


@SymmetricRollPass.roll_torque_loss_function
def roll_torque_loss_function(self: SymmetricRollPass):
    if np.isclose(self.rel_draught, 0):
        return 1

    abs_rel_draught = abs(self.rel_draught)

    return np.sqrt(self.roll.working_radius / self.out_profile.equivalent_height) * np.sqrt(
        (1 - abs_rel_draught) / abs_rel_draught) * (1 / 2 - self.roll.relative_neutral_angle)


@SymmetricRollPass.Roll.roll_torque
def roll_torque(self: SymmetricRollPass.Roll):
    rp = self.roll_pass
    mean_flow_stress = (rp.in_profile.flow_stress + 2 * rp.out_profile.flow_stress) / 3
    mean_width = (rp.in_profile.equivalent_width + 2 * rp.out_profile.equivalent_width) / 3
    height_change = rp.in_profile.equivalent_height - rp.out_profile.equivalent_height

    return mean_width * self.working_radius * mean_flow_stress * height_change * rp.roll_torque_loss_function
