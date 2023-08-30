import logging
import numpy as np
from pyroll.core import RollPass, Roll, Hook

VERSION = "2.0.0.post0"

log = logging.getLogger(__name__)

RollPass.inverse_forming_efficiency = Hook[float]()
"""Inverse forming efficiency defined by Lippmann and Mahrenholz for the roll force of the roll pass."""

RollPass.roll_torque_loss_function = Hook[float]()
"""Loss function defined by Lippmann and Mahrenholz for the roll torque for the roll pass."""

Roll.entry_angle = Hook[float]()
"""Angle at which the profile enters the roll gap."""

Roll.neutral_angle = Hook[float]()
"""Neutral angle defined by Lippmann and Mahrenholz."""


@RollPass.front_tension
def front_tension(self: RollPass):
    return 0


@RollPass.back_tension
def back_tension(self: RollPass):
    return 0


@RollPass.Roll.entry_angle
def entry_angle(self: RollPass.Roll):
    return -np.arcsin(-self.contact_length / self.working_radius)


@RollPass.Roll.neutral_angle
def neutral_angle(self: RollPass.Roll):
    rp = self.roll_pass
    mean_flow_stress = (rp.in_profile.flow_stress + 2 * rp.out_profile.flow_stress) / 3
    abs_rel_drought = np.abs(rp.rel_draught)

    p1 = np.sqrt((1 - abs_rel_drought) / abs_rel_drought)
    p2 = 1 / 2 * np.sqrt(rp.out_profile.equivalent_height / self.working_radius)
    p3 = np.log(1 - abs_rel_drought)

    return self.entry_angle * np.tan(
        p2 * ((rp.back_tension - rp.front_tension) / mean_flow_stress + p3) + 0.5 * np.arctan(p1))


@RollPass.Roll.neutral_point
def neutral_point(self: RollPass.Roll):
    return -self.working_radius * np.sin(self.neutral_angle)


@RollPass.inverse_forming_efficiency
def inverse_forming_efficiency(self: RollPass):
    if np.isclose(self.rel_draught, 0):
        return 1

    mean_flow_stress = (self.in_profile.flow_stress + 2 * self.out_profile.flow_stress) / 3
    abs_rel_drought = np.abs(self.rel_draught)
    relative_neutral_angle = self.roll.neutral_angle / self.roll.entry_angle

    return self.back_tension / mean_flow_stress + 2 * np.sqrt(
        (1 - abs_rel_drought) / abs_rel_drought) * np.arctan(
        np.sqrt(abs_rel_drought / (1 - abs_rel_drought))) - 1 + np.sqrt(
        self.roll.working_radius / self.out_profile.equivalent_height) * np.sqrt(
        (1 - abs_rel_drought) / abs_rel_drought) * np.log(
        np.sqrt(1 - abs_rel_drought) / (1 - abs_rel_drought * (1 - relative_neutral_angle ** 2)))


@RollPass.DiskElement.deformation_resistance
def deformation_resistance(self: RollPass.DiskElement):
    mean_flow_stress = (self.in_profile.flow_stress + 2 * self.out_profile.flow_stress) / 2
    return mean_flow_stress * self.roll_pass.inverse_forming_efficiency


@RollPass.deformation_resistance
def deformation_resistance(self: RollPass):
    mean_flow_stress = (self.in_profile.flow_stress + 2 * self.out_profile.flow_stress) / 3
    return mean_flow_stress * self.inverse_forming_efficiency


@RollPass.roll_force
def roll_force(self: RollPass):
    return self.roll.contact_area * self.deformation_resistance


@RollPass.roll_torque_loss_function
def roll_torque_loss_function(self: RollPass):
    if np.isclose(self.rel_draught, 0):
        return 1

    relative_neutral_angle = self.roll.neutral_angle / self.roll.entry_angle
    abs_rel_drought = np.abs(self.rel_draught)

    return np.sqrt(self.roll.working_radius / self.in_profile.equivalent_height) * np.sqrt(
        (1 - abs_rel_drought) / abs_rel_drought) * (1 / 2 - relative_neutral_angle)


@RollPass.Roll.roll_torque
def roll_torque(self: RollPass.Roll):
    rp = self.roll_pass
    mean_flow_stress = (rp.in_profile.flow_stress + 2 * rp.out_profile.flow_stress) / 3
    mean_width = (rp.in_profile.equivalent_width + 2 * rp.out_profile.equivalent_width) / 3
    height_change = rp.in_profile.equivalent_height - rp.out_profile.equivalent_height

    return mean_width * self.working_radius * mean_flow_stress * height_change * self.roll_pass.roll_torque_loss_function
