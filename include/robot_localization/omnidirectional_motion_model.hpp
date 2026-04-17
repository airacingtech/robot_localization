/*
 * Copyright (c) 2025 AIRacing Tech
 * Omnidirectional motion model (original RL model)
 */
#ifndef ROBOT_LOCALIZATION__OMNIDIRECTIONAL_MOTION_MODEL_HPP_
#define ROBOT_LOCALIZATION__OMNIDIRECTIONAL_MOTION_MODEL_HPP_

#include "robot_localization/motion_model.hpp"

namespace robot_localization
{

/**
 * @brief Omnidirectional motion model
 *
 * Assumes 3D motion with velocities independent of heading.
 * Velocity can be in any direction regardless of which way vehicle is pointing.
 * Appropriate for drones, wheeled robots with omnidirectional wheels.
 */
class OmnidirectionalMotionModel : public MotionModel
{
public:
  void computeStateTransition(
    const Eigen::VectorXd & state,
    double delta_sec,
    Eigen::MatrixXd & transfer_function,
    Eigen::MatrixXd & transfer_function_jacobian) override;
};

}  // namespace robot_localization

#endif  // ROBOT_LOCALIZATION__OMNIDIRECTIONAL_MOTION_MODEL_HPP_
