import logging
import numpy as np
from pyroll.core import RollPass, Hook, root_hooks

VERSION = "2.0.0"

log = logging.getLogger(__name__)

RollPass.inverse_forming_efficiency = Hook[float]()
"""Inverse forming efficiency of the roll pass."""

RollPass.roll_torque_loss_function = Hook[float]()
"""Loss function defined by Lippmann and Mahrenholz for the roll torque for the roll pass."""


@RollPass.mean_front_tension
def mean_front_tension(self: RollPass):
    raise ValueError(
        "You must provide a mean front tension to use the pyroll-lippmann-mahrenholz-power-and-labour plugin!")


@RollPass.mean_back_tension
def mean_back_tension(self: RollPass):
    raise ValueError(
        "You must provide a mean front tension to use the pyroll-lippmann-mahrenholz-power-and-labour plugin!")


@RollPass.mean_neutral_plane_angle
def mean_neutral_plane_angle(self: RollPass):

    mean_flow_stress = (self.in_profile.flow_stress + 2 * self.out_profile.flow_stress) / 3
    abs_rel_drought = np.abs(self.rel_draught)

    return np.sqrt((1 - abs_rel_drought) / abs_rel_drought) * np.tan(
        1 / 2 * np.sqrt(self.out_profile.equivalent_height / self.roll.working_radius) * (
                (self.mean_back_tension - self.mean_front_tension) / mean_flow_stress + np.log(
            1 - abs_rel_drought)) + 1 / 2 * np.arctan(np.sqrt(abs_rel_drought / (1 - abs_rel_drought))))


@RollPass.inverse_forming_efficiency
def inverse_forming_efficiency(self: RollPass):
    if np.isclose(self.rel_draught, 0):
        return 1

    mean_flow_stress = (self.in_profile.flow_stress + 2 * self.out_profile.flow_stress) / 3
    abs_rel_drought = np.abs(self.rel_draught)

    return self.mean_back_tension / mean_flow_stress + 2 * np.sqrt(
        (1 - abs_rel_drought) / abs_rel_drought) * np.arctan(
        np.sqrt(abs_rel_drought / (1 - abs_rel_drought))) - 1 + np.sqrt(
        self.roll.working_radius / self.out_profile.equivalent_height) * np.sqrt(
        (1 - abs_rel_drought) / abs_rel_drought) * np.log(
        np.sqrt(1 - abs_rel_drought) / (1 - abs_rel_drought * (1 - self.mean_neutral_plane_angle ** 2)))


@RollPass.roll_force
def roll_force(self: RollPass):
    mean_flow_stress = (self.in_profile.flow_stress + 2 * self.out_profile.flow_stress) / 3
    return self.roll.contact_area * mean_flow_stress * self.inverse_forming_efficiency


@RollPass.roll_torque_loss_function
def roll_torque_loss_function(self: RollPass):

    if np.isclose(self.rel_draught, 0):
        return 1

    abs_rel_drought = np.abs(self.rel_draught)
    return np.sqrt(self.roll.working_radius / self.in_profile.equivalent_height) * np.sqrt(
        (1 - abs_rel_drought) / abs_rel_drought) * (1 / 2 - self.mean_neutral_plane_angle)


@RollPass.Roll.roll_torque
def roll_torque(self: RollPass.Roll):
    mean_flow_stress = (self.roll_pass.in_profile.flow_stress + 2 * self.roll_pass.out_profile.flow_stress) / 3
    return self.roll_pass.out_profile.equivalent_width * mean_flow_stress * self.contact_length ** 2 * self.roll_pass.roll_torque_loss_function
