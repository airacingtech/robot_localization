/*
 * Copyright (c) 2025 AIRacing Tech
 */
#include "robot_localization/omnidirectional_motion_model.hpp"

#include <cmath>

#include "robot_localization/filter_common.hpp"

namespace robot_localization
{

void OmnidirectionalMotionModel::computeStateTransition(
  const Eigen::VectorXd & state,
  double delta_sec,
  Eigen::MatrixXd & transfer_function,
  Eigen::MatrixXd & transfer_function_jacobian)
{
  transfer_function.setIdentity();
  transfer_function_jacobian.setIdentity();

  double roll = state(StateMemberRoll);
  double pitch = state(StateMemberPitch);
  double yaw = state(StateMemberYaw);
  double x_vel = state(StateMemberVx);
  double y_vel = state(StateMemberVy);
  double z_vel = state(StateMemberVz);
  double pitch_vel = state(StateMemberVpitch);
  double yaw_vel = state(StateMemberVyaw);
  double x_acc = state(StateMemberAx);
  double y_acc = state(StateMemberAy);
  double z_acc = state(StateMemberAz);

  double sp = ::sin(pitch);
  double cp = ::cos(pitch);
  double cpi = 1.0 / cp;
  double tp = sp * cpi;

  double sr = ::sin(roll);
  double cr = ::cos(roll);

  double sy = ::sin(yaw);
  double cy = ::cos(yaw);

  // Transfer function
  transfer_function(StateMemberX, StateMemberVx) = cy * cp * delta_sec;
  transfer_function(StateMemberX, StateMemberVy) =
    (cy * sp * sr - sy * cr) * delta_sec;
  transfer_function(StateMemberX, StateMemberVz) =
    (cy * sp * cr + sy * sr) * delta_sec;
  transfer_function(StateMemberX, StateMemberAx) =
    0.5 * transfer_function(StateMemberX, StateMemberVx) * delta_sec;
  transfer_function(StateMemberX, StateMemberAy) =
    0.5 * transfer_function(StateMemberX, StateMemberVy) * delta_sec;
  transfer_function(StateMemberX, StateMemberAz) =
    0.5 * transfer_function(StateMemberX, StateMemberVz) * delta_sec;
  transfer_function(StateMemberY, StateMemberVx) = sy * cp * delta_sec;
  transfer_function(StateMemberY, StateMemberVy) =
    (sy * sp * sr + cy * cr) * delta_sec;
  transfer_function(StateMemberY, StateMemberVz) =
    (sy * sp * cr - cy * sr) * delta_sec;
  transfer_function(StateMemberY, StateMemberAx) =
    0.5 * transfer_function(StateMemberY, StateMemberVx) * delta_sec;
  transfer_function(StateMemberY, StateMemberAy) =
    0.5 * transfer_function(StateMemberY, StateMemberVy) * delta_sec;
  transfer_function(StateMemberY, StateMemberAz) =
    0.5 * transfer_function(StateMemberY, StateMemberVz) * delta_sec;
  transfer_function(StateMemberZ, StateMemberVx) = -sp * delta_sec;
  transfer_function(StateMemberZ, StateMemberVy) = cp * sr * delta_sec;
  transfer_function(StateMemberZ, StateMemberVz) = cp * cr * delta_sec;
  transfer_function(StateMemberZ, StateMemberAx) =
    0.5 * transfer_function(StateMemberZ, StateMemberVx) * delta_sec;
  transfer_function(StateMemberZ, StateMemberAy) =
    0.5 * transfer_function(StateMemberZ, StateMemberVy) * delta_sec;
  transfer_function(StateMemberZ, StateMemberAz) =
    0.5 * transfer_function(StateMemberZ, StateMemberVz) * delta_sec;
  transfer_function(StateMemberRoll, StateMemberVroll) = delta_sec;
  transfer_function(StateMemberRoll, StateMemberVpitch) = sr * tp * delta_sec;
  transfer_function(StateMemberRoll, StateMemberVyaw) = cr * tp * delta_sec;
  transfer_function(StateMemberPitch, StateMemberVpitch) = cr * delta_sec;
  transfer_function(StateMemberPitch, StateMemberVyaw) = -sr * delta_sec;
  transfer_function(StateMemberYaw, StateMemberVpitch) = sr * cpi * delta_sec;
  transfer_function(StateMemberYaw, StateMemberVyaw) = cr * cpi * delta_sec;
  transfer_function(StateMemberVx, StateMemberAx) = delta_sec;
  transfer_function(StateMemberVy, StateMemberAy) = delta_sec;
  transfer_function(StateMemberVz, StateMemberAz) = delta_sec;

  // Jacobian (analytical derivatives)
  double one_half_at_squared = 0.5 * delta_sec * delta_sec;

  double y_coeff = cy * sp * cr + sy * sr;
  double z_coeff = -cy * sp * sr + sy * cr;
  double dFx_dR = (y_coeff * y_vel + z_coeff * z_vel) * delta_sec +
    (y_coeff * y_acc + z_coeff * z_acc) * one_half_at_squared;
  double dFR_dR = 1.0 + (cr * tp * pitch_vel - sr * tp * yaw_vel) * delta_sec;

  double x_coeff = -cy * sp;
  y_coeff = cy * cp * sr;
  z_coeff = cy * cp * cr;
  double dFx_dP =
    (x_coeff * x_vel + y_coeff * y_vel + z_coeff * z_vel) * delta_sec +
    (x_coeff * x_acc + y_coeff * y_acc + z_coeff * z_acc) *
    one_half_at_squared;
  double dFR_dP =
    (cpi * cpi * sr * pitch_vel + cpi * cpi * cr * yaw_vel) * delta_sec;

  x_coeff = -sy * cp;
  y_coeff = -sy * sp * sr - cy * cr;
  z_coeff = -sy * sp * cr + cy * sr;
  double dFx_dY =
    (x_coeff * x_vel + y_coeff * y_vel + z_coeff * z_vel) * delta_sec +
    (x_coeff * x_acc + y_coeff * y_acc + z_coeff * z_acc) *
    one_half_at_squared;

  y_coeff = sy * sp * cr - cy * sr;
  z_coeff = -sy * sp * sr - cy * cr;
  double dFy_dR = (y_coeff * y_vel + z_coeff * z_vel) * delta_sec +
    (y_coeff * y_acc + z_coeff * z_acc) * one_half_at_squared;
  double dFP_dR = (-sr * pitch_vel - cr * yaw_vel) * delta_sec;

  x_coeff = -sy * sp;
  y_coeff = sy * cp * sr;
  z_coeff = sy * cp * cr;
  double dFy_dP =
    (x_coeff * x_vel + y_coeff * y_vel + z_coeff * z_vel) * delta_sec +
    (x_coeff * x_acc + y_coeff * y_acc + z_coeff * z_acc) *
    one_half_at_squared;

  x_coeff = cy * cp;
  y_coeff = cy * sp * sr - sy * cr;
  z_coeff = cy * sp * cr + sy * sr;
  double dFy_dY =
    (x_coeff * x_vel + y_coeff * y_vel + z_coeff * z_vel) * delta_sec +
    (x_coeff * x_acc + y_coeff * y_acc + z_coeff * z_acc) *
    one_half_at_squared;

  y_coeff = cp * cr;
  z_coeff = -cp * sr;
  double dFz_dR = (y_coeff * y_vel + z_coeff * z_vel) * delta_sec +
    (y_coeff * y_acc + z_coeff * z_acc) * one_half_at_squared;
  double dFY_dR = (cr * cpi * pitch_vel - sr * cpi * yaw_vel) * delta_sec;

  x_coeff = -cp;
  y_coeff = -sp * sr;
  z_coeff = -sp * cr;
  double dFz_dP =
    (x_coeff * x_vel + y_coeff * y_vel + z_coeff * z_vel) * delta_sec +
    (x_coeff * x_acc + y_coeff * y_acc + z_coeff * z_acc) *
    one_half_at_squared;
  double dFY_dP =
    (sr * tp * cpi * pitch_vel + cr * tp * cpi * yaw_vel) * delta_sec;

  transfer_function_jacobian = transfer_function;
  transfer_function_jacobian(StateMemberX, StateMemberRoll) = dFx_dR;
  transfer_function_jacobian(StateMemberX, StateMemberPitch) = dFx_dP;
  transfer_function_jacobian(StateMemberX, StateMemberYaw) = dFx_dY;
  transfer_function_jacobian(StateMemberY, StateMemberRoll) = dFy_dR;
  transfer_function_jacobian(StateMemberY, StateMemberPitch) = dFy_dP;
  transfer_function_jacobian(StateMemberY, StateMemberYaw) = dFy_dY;
  transfer_function_jacobian(StateMemberZ, StateMemberRoll) = dFz_dR;
  transfer_function_jacobian(StateMemberZ, StateMemberPitch) = dFz_dP;
  transfer_function_jacobian(StateMemberRoll, StateMemberRoll) = dFR_dR;
  transfer_function_jacobian(StateMemberRoll, StateMemberPitch) = dFR_dP;
  transfer_function_jacobian(StateMemberPitch, StateMemberRoll) = dFP_dR;
  transfer_function_jacobian(StateMemberYaw, StateMemberRoll) = dFY_dR;
  transfer_function_jacobian(StateMemberYaw, StateMemberPitch) = dFY_dP;
}

}  // namespace robot_localization
