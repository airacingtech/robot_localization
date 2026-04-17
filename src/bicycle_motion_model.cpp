/*
 * Copyright (c) 2025 AIRacing Tech
 * Bicycle motion model implementation
 */
#include "robot_localization/bicycle_motion_model.hpp"

#include <cmath>

#include "robot_localization/filter_common.hpp"

namespace robot_localization
{

void BicycleMotionModel::setParameters(double wheelbase_m, double steering_angle_rad)
{
  wheelbase_m_ = wheelbase_m;
  steering_angle_rad_ = steering_angle_rad;
}

void BicycleMotionModel::computeStateTransition(
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
  double vx = state(StateMemberVx);
  double vy = state(StateMemberVy);
  double z_vel = state(StateMemberVz);
  double pitch_vel = state(StateMemberVpitch);
  double yaw_vel = state(StateMemberVyaw);
  double x_acc = state(StateMemberAx);
  double y_acc = state(StateMemberAy);
  double z_acc = state(StateMemberAz);

  // Bicycle kinematics: position moves in direction of yaw
  // Linear velocity is longitudinal (x-direction in body frame)
  // Steering angle couples yaw rate to velocity: dyaw/dt = v * tan(delta) / L
  double v_longitudinal = std::sqrt(vx * vx + vy * vy);  // Total velocity magnitude
  double delta = steering_angle_rad_;
  double L = wheelbase_m_;

  double sp = ::sin(pitch);
  double cp = ::cos(pitch);
  double sr = ::sin(roll);
  double cr = ::cos(roll);
  double sy = ::sin(yaw);
  double cy = ::cos(yaw);

  double one_half_at_squared = 0.5 * delta_sec * delta_sec;

  // Position update: x and y move in direction of yaw with longitudinal velocity
  // x_new = x + v * cos(yaw) * dt
  // y_new = y + v * sin(yaw) * dt
  transfer_function(StateMemberX, StateMemberVx) = cy * delta_sec;
  transfer_function(StateMemberX, StateMemberVy) = -sy * delta_sec;
  transfer_function(StateMemberX, StateMemberAx) =
    0.5 * transfer_function(StateMemberX, StateMemberVx) * delta_sec;
  transfer_function(StateMemberX, StateMemberAy) =
    0.5 * transfer_function(StateMemberX, StateMemberVy) * delta_sec;

  transfer_function(StateMemberY, StateMemberVx) = sy * delta_sec;
  transfer_function(StateMemberY, StateMemberVy) = cy * delta_sec;
  transfer_function(StateMemberY, StateMemberAx) =
    0.5 * transfer_function(StateMemberY, StateMemberVx) * delta_sec;
  transfer_function(StateMemberY, StateMemberAy) =
    0.5 * transfer_function(StateMemberY, StateMemberVy) * delta_sec;

  // Z position (unchanged for planar motion)
  transfer_function(StateMemberZ, StateMemberVz) = delta_sec;
  transfer_function(StateMemberZ, StateMemberAz) =
    0.5 * transfer_function(StateMemberZ, StateMemberVz) * delta_sec;

  // Roll/Pitch/Yaw angles (3D Euler integration)
  transfer_function(StateMemberRoll, StateMemberVroll) = delta_sec;
  transfer_function(StateMemberRoll, StateMemberVpitch) = sr * sp * delta_sec / cp;
  transfer_function(StateMemberRoll, StateMemberVyaw) = cr * sp * delta_sec / cp;
  transfer_function(StateMemberPitch, StateMemberVpitch) = cr * delta_sec;
  transfer_function(StateMemberPitch, StateMemberVyaw) = -sr * delta_sec;

  // Yaw update: bicycle model coupling with steering
  // dyaw/dt = v * tan(delta) / L
  // yaw_new = yaw + v * tan(delta) / L * dt
  double tan_delta = std::tan(delta);
  double yaw_rate_from_steering = v_longitudinal * tan_delta / L;

  transfer_function(StateMemberYaw, StateMemberVpitch) = sr * delta_sec / cp;
  transfer_function(StateMemberYaw, StateMemberVyaw) = cr * delta_sec;
  transfer_function(StateMemberYaw, StateMemberVx) = tan_delta / L * delta_sec;
  // vy should NOT affect yaw rate in bicycle model - only longitudinal velocity matters
  transfer_function(StateMemberYaw, StateMemberVy) = 0.0;

  // Velocity updates
  transfer_function(StateMemberVx, StateMemberAx) = delta_sec;
  transfer_function(StateMemberVy, StateMemberAy) = delta_sec;
  transfer_function(StateMemberVz, StateMemberAz) = delta_sec;

  // Jacobian: mostly same as transfer function for linear terms
  transfer_function_jacobian = transfer_function;

  // Add nonlinear derivatives for yaw coupling
  // d(yaw)/d(steering) = v / L (constant, but handled by steering subscriber)
  // d(yaw_pred)/d(vx) = tan(delta) / L * dt (already in transfer function)
  // d(x)/d(yaw) depends on velocity
  double dx_dyaw = -sy * v_longitudinal * delta_sec;
  double dy_dyaw = cy * v_longitudinal * delta_sec;
  transfer_function_jacobian(StateMemberX, StateMemberYaw) = dx_dyaw;
  transfer_function_jacobian(StateMemberY, StateMemberYaw) = dy_dyaw;
}

}  // namespace robot_localization
