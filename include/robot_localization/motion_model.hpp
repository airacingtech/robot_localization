/*
 * Copyright (c) 2025 AIRacing Tech
 * Motion model abstraction for configurable state transition models
 */
#ifndef ROBOT_LOCALIZATION__MOTION_MODEL_HPP_
#define ROBOT_LOCALIZATION__MOTION_MODEL_HPP_

#include "Eigen/Dense"

namespace robot_localization
{

/**
 * @brief Abstract base class for motion models
 *
 * Defines the interface for computing state transition matrices (F) and
 * their Jacobians (H) for different vehicle motion models.
 */
class MotionModel
{
public:
  virtual ~MotionModel() = default;

  /**
   * @brief Compute the state transition matrix F and its Jacobian
   *
   * @param[in] state Current state vector (15 elements)
   * @param[in] delta_sec Time delta in seconds
   * @param[out] transfer_function State transition matrix F (15x15)
   * @param[out] transfer_function_jacobian Jacobian of F (15x15)
   */
  virtual void computeStateTransition(
    const Eigen::VectorXd & state,
    double delta_sec,
    Eigen::MatrixXd & transfer_function,
    Eigen::MatrixXd & transfer_function_jacobian) = 0;
};

}  // namespace robot_localization

#endif  // ROBOT_LOCALIZATION__MOTION_MODEL_HPP_
