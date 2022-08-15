import logging
import numpy as np
from pyroll.core import RollPass

log = logging.getLogger(__name__)


@RollPass.hookimpl
def forward_tension(roll_pass: RollPass):
    return 0


@RollPass.hookimpl
def backward_tension(roll_pass: RollPass):
    return 0


@RollPass.hookimpl
def equivalent_reduction(roll_pass: RollPass):
    return np.abs((roll_pass.out_profile.equivalent_rectangle.height - roll_pass.in_profile.equivalent_rectangle.height) / roll_pass.in_profile.equivalent_rectangle.height)


@RollPass.hookimpl
def neutral_line_angle(roll_pass: RollPass):
    return np.sqrt((1 - roll_pass.equivalent_reduction) / roll_pass.equivalent_reduction) * np.tan(
        1 / 2 * np.sqrt(roll_pass.out_profile.equivalent_rectangle.height / roll_pass.roll.working_radius) * (
                (roll_pass.backward_tension - roll_pass.forward_tension) / roll_pass.mean_flow_stress + np.log(
            1 - roll_pass.equivalent_reduction)) + 1 / 2 *
        np.arctan(np.sqrt(roll_pass.equivalent_reduction / (1 - roll_pass.equivalent_reduction))))


@RollPass.hookimpl
def inverse_forming_efficiency(roll_pass: RollPass):
    return roll_pass.backward_tension / roll_pass.mean_flow_stress + 2 * \
           np.sqrt((1 - roll_pass.equivalent_reduction) / roll_pass.equivalent_reduction) * \
           np.arctan(np.sqrt(roll_pass.equivalent_reduction / (1 - roll_pass.equivalent_reduction))) \
           - 1 + np.sqrt(roll_pass.roll.working_radius / roll_pass.out_profile.equivalent_rectangle.height) \
           * np.sqrt((1 - roll_pass.equivalent_reduction) / roll_pass.equivalent_reduction) * \
           np.log(np.sqrt(1 - roll_pass.equivalent_reduction) / (
                   1 - roll_pass.equivalent_reduction * (1 - roll_pass.neutral_line_angle ** 2)))


@RollPass.hookimpl
def roll_force(roll_pass: RollPass):
    return roll_pass.roll.contact_area * roll_pass.mean_flow_stress * roll_pass.inverse_forming_efficiency


@RollPass.hookimpl
def lippmann_mahrenholz_roll_torque_function(roll_pass: RollPass):
    return np.sqrt(roll_pass.roll.working_radius / roll_pass.in_profile.equivalent_rectangle.height) * np.sqrt(
        (1 - roll_pass.equivalent_reduction) / roll_pass.equivalent_reduction) * (1 / 2 - roll_pass.neutral_line_angle)


@RollPass.Roll.hookimpl
def roll_torque(roll_pass: RollPass, roll: RollPass.Roll):
    return roll_pass.out_profile.equivalent_rectangle.width * roll_pass.mean_flow_stress * roll_pass.roll.contact_length ** 2 * roll_pass.lippmann_mahrenholz_roll_torque_function
