/*
 * Copyright (c) 2025 AIRacing Tech
 * Bicycle (Ackermann) motion model for wheeled vehicles
 */
#ifndef ROBOT_LOCALIZATION__BICYCLE_MOTION_MODEL_HPP_
#define ROBOT_LOCALIZATION__BICYCLE_MOTION_MODEL_HPP_

#include "robot_localization/motion_model.hpp"

namespace robot_localization
{

/**
 * @brief Kinematic bicycle motion model
 *
 * Models a vehicle with front-wheel steering (Ackermann geometry).
 * Heading and velocity are coupled: vehicle moves in direction it's pointing.
 * Assumes: no slip, steering angle input.
 *
 * State: [x, y, z, roll, pitch, yaw, vx, vy, vz, vroll, vpitch, vyaw, ax, ay, az]
 * Input: steering angle (delta), measured separately
 */
class BicycleMotionModel : public MotionModel
{
public:
  /**
   * @brief Set bicycle model parameters
   * @param[in] wheelbase_m Distance between front and rear axles (meters)
   * @param[in] steering_angle_rad Current front wheel steering angle (radians)
   */
  void setParameters(double wheelbase_m, double steering_angle_rad);

  void computeStateTransition(
    const Eigen::VectorXd & state,
    double delta_sec,
    Eigen::MatrixXd & transfer_function,
    Eigen::MatrixXd & transfer_function_jacobian) override;

private:
  double wheelbase_m_ = 2.6;           // Default wheelbase
  double steering_angle_rad_ = 0.0;    // Current steering angle
};

}  // namespace robot_localization

#endif  // ROBOT_LOCALIZATION__BICYCLE_MOTION_MODEL_HPP_
