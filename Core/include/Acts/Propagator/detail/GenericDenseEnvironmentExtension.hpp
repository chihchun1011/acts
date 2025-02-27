// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <array>
#include <functional>

namespace Acts {
namespace detail {

/// @brief Evaluater of the k_i's and elements of the transport matrix
/// D of the RKN4 stepping. This implementation involves energy loss due to
/// ioninisation, bremsstrahlung, pair production and photonuclear interaction
/// in the propagation and the jacobian. These effects will only occur if the
/// propagation is in a TrackingVolume with attached material.
/// @note This it templated on the scalar type because of the autodiff plugin.
template <typename scalar_t>
struct GenericDenseEnvironmentExtension {
  using Scalar = scalar_t;
  /// @brief Vector3 replacement for the custom scalar type
  using ThisVector3 = Eigen::Matrix<Scalar, 3, 1>;

  /// Momentum at a certain point
  Scalar currentMomentum = 0.;
  /// Particles momentum at k1
  Scalar initialMomentum = 0.;
  /// Material that will be passed
  /// TODO : Might not be needed anymore
  Material material;
  /// Derivatives dLambda''dlambda at each sub-step point
  std::array<Scalar, 4> dLdl{};
  /// q/p at each sub-step
  std::array<Scalar, 4> qop{};
  /// Derivatives dPds at each sub-step
  std::array<Scalar, 4> dPds{};
  /// Derivative d(dEds)d(q/p) evaluated at the initial point
  Scalar dgdqopValue = 0.;
  /// Derivative dEds at the initial point
  Scalar g = 0.;
  /// k_i equivalent for the time propagation
  std::array<Scalar, 4> tKi{};
  /// Lambda''_i
  std::array<Scalar, 4> Lambdappi{};
  /// Energy at each sub-step
  std::array<Scalar, 4> energy{};

  /// @brief Default constructor
  GenericDenseEnvironmentExtension() = default;

  /// @brief Control function if the step evaluation would be valid
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @return Boolean flag if the step would be valid
  template <typename propagator_state_t, typename stepper_t>
  int bid(const propagator_state_t& state, const stepper_t& stepper) const {
    // Check for valid particle properties
    if (stepper.charge(state.stepping) == 0. || state.options.mass == 0. ||
        stepper.momentum(state.stepping) < state.options.momentumCutOff) {
      return 0;
    }

    // Check existence of a volume with material
    if (!state.navigation.currentVolume ||
        !state.navigation.currentVolume->volumeMaterial()) {
      return 0;
    }
    return 2;
  }

  /// @brief Evaluater of the k_i's of the RKN4. For the case of i = 0 this
  /// step sets up member parameters, too.
  ///
  /// @tparam stepper_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [out] knew Next k_i that is evaluated
  /// @param [out] kQoP k_i elements of the momenta
  /// @param [in] bField B-Field at the evaluation position
  /// @param [in] i Index of the k_i, i = [0, 3]
  /// @param [in] h Step size (= 0. ^ 0.5 * StepSize ^ StepSize)
  /// @param [in] kprev Evaluated k_{i - 1}
  /// @return Boolean flag if the calculation is valid
  template <typename propagator_state_t, typename stepper_t>
  bool k(const propagator_state_t& state, const stepper_t& stepper,
         ThisVector3& knew, const Vector3& bField, std::array<Scalar, 4>& kQoP,
         const int i = 0, const double h = 0.,
         const ThisVector3& kprev = ThisVector3::Zero()) {
    // i = 0 is used for setup and evaluation of k
    if (i == 0) {
      // Set up container for energy loss
      auto volumeMaterial = state.navigation.currentVolume->volumeMaterial();
      ThisVector3 position = stepper.position(state.stepping);
      material = (volumeMaterial->material(position.template cast<double>()));
      initialMomentum = stepper.momentum(state.stepping);
      currentMomentum = initialMomentum;
      qop[0] = stepper.charge(state.stepping) / initialMomentum;
      initializeEnergyLoss(state);
      // Evaluate k
      knew = qop[0] * stepper.direction(state.stepping).cross(bField);
      // Evaluate k for the time propagation
      Lambdappi[0] =
          -qop[0] * qop[0] * qop[0] * g * energy[0] /
          (stepper.charge(state.stepping) * stepper.charge(state.stepping));
      //~ tKi[0] = std::hypot(1, state.options.mass / initialMomentum);
      using std::hypot;
      tKi[0] = hypot(1, state.options.mass * qop[0]);
      kQoP[0] = Lambdappi[0];
    } else {
      // Update parameters and check for momentum condition
      updateEnergyLoss(state.options.mass, h, state.stepping, stepper, i);
      if (currentMomentum < state.options.momentumCutOff) {
        return false;
      }
      // Evaluate k
      knew = qop[i] *
             (stepper.direction(state.stepping) + h * kprev).cross(bField);
      // Evaluate k_i for the time propagation
      auto qopNew = qop[0] + h * Lambdappi[i - 1];
      Lambdappi[i] =
          -qopNew * qopNew * qopNew * g * energy[i] /
          (stepper.charge(state.stepping) * stepper.charge(state.stepping) *
           UnitConstants::C * UnitConstants::C);
      using std::hypot;
      tKi[i] = hypot(1, state.options.mass * qopNew);
      kQoP[i] = Lambdappi[i];
    }
    return true;
  }

  /// @brief After a RKN4 step was accepted by the stepper this method has an
  /// additional veto on the quality of the step. The veto lies in evaluation
  /// of the energy loss and the therewith constrained to keep the momentum
  /// after the step in reasonable values.
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [in] h Step size
  /// @return Boolean flag if the calculation is valid
  template <typename propagator_state_t, typename stepper_t>
  bool finalize(propagator_state_t& state, const stepper_t& stepper,
                const double h) const {
    // Evaluate the new momentum
    auto newMomentum =
        stepper.momentum(state.stepping) +
        (h / 6.) * (dPds[0] + 2. * (dPds[1] + dPds[2]) + dPds[3]);

    // Break propagation if momentum becomes below cut-off
    if (newMomentum < state.options.momentumCutOff) {
      return false;
    }

    // Add derivative dlambda/ds = Lambda''
    using std::sqrt;
    state.stepping.derivative(7) =
        -sqrt(state.options.mass * state.options.mass +
              newMomentum * newMomentum) *
        g / (newMomentum * newMomentum * newMomentum);

    // Update momentum
    state.stepping.pars[eFreeQOverP] =
        stepper.charge(state.stepping) / newMomentum;
    // Add derivative dt/ds = 1/(beta * c) = sqrt(m^2 * p^{-2} + c^{-2})
    using std::hypot;
    state.stepping.derivative(3) = hypot(1, state.options.mass / newMomentum);
    // Update time
    state.stepping.pars[eFreeTime] +=
        (h / 6.) * (tKi[0] + 2. * (tKi[1] + tKi[2]) + tKi[3]);

    return true;
  }

  /// @brief After a RKN4 step was accepted by the stepper this method has an
  /// additional veto on the quality of the step. The veto lies in the
  /// evaluation
  /// of the energy loss, the therewith constrained to keep the momentum
  /// after the step in reasonable values and the evaluation of the transport
  /// matrix.
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [in] h Step size
  /// @param [out] D Transport matrix
  /// @return Boolean flag if the calculation is valid
  template <typename propagator_state_t, typename stepper_t>
  bool finalize(propagator_state_t& state, const stepper_t& stepper,
                const double h, FreeMatrix& D) const {
    return finalize(state, stepper, h) && transportMatrix(state, stepper, h, D);
  }

 private:
  /// @brief Evaluates the transport matrix D for the jacobian
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [in] h Step size
  /// @param [out] D Transport matrix
  /// @return Boolean flag if evaluation is valid
  template <typename propagator_state_t, typename stepper_t>
  bool transportMatrix(propagator_state_t& state, const stepper_t& stepper,
                       const double h, FreeMatrix& D) const {
    /// The calculations are based on ATL-SOFT-PUB-2009-002. The update of the
    /// Jacobian matrix is requires only the calculation of eq. 17 and 18.
    /// Since the terms of eq. 18 are currently 0, this matrix is not needed
    /// in the calculation. The matrix A from eq. 17 consists out of 3
    /// different parts. The first one is given by the upper left 3x3 matrix
    /// that are calculated by dFdT and dGdT. The second is given by the top 3
    /// lines of the rightmost column. This is calculated by dFdL and dGdL.
    /// The remaining non-zero term is calculated directly. The naming of the
    /// variables is explained in eq. 11 and are directly related to the
    /// initial problem in eq. 7.
    /// The evaluation is based on propagating the parameters T and lambda
    /// (including g(lambda) and E(lambda)) as given in eq. 16 and evaluating
    /// the derivations for matrix A.
    /// @note The translation for u_{n+1} in eq. 7 is in this case a
    /// 3-dimensional vector without a dependency of Lambda or lambda neither in
    /// u_n nor in u_n'. The second and fourth eq. in eq. 14 have the constant
    /// offset matrices h * Id and Id respectively. This involves that the
    /// constant offset does not exist for rectangular matrix dFdu' (due to the
    /// missing Lambda part) and only exists for dGdu' in dlambda/dlambda.

    auto& sd = state.stepping.stepData;
    auto dir = stepper.direction(state.stepping);

    D = FreeMatrix::Identity();
    const double half_h = h * 0.5;

    // This sets the reference to the sub matrices
    // dFdx is already initialised as (3x3) zero
    auto dFdT = D.block<3, 3>(0, 4);
    auto dFdL = D.block<3, 1>(0, 7);
    // dGdx is already initialised as (3x3) identity
    auto dGdT = D.block<3, 3>(4, 4);
    auto dGdL = D.block<3, 1>(4, 7);

    ActsMatrix<3, 3> dk1dT = ActsMatrix<3, 3>::Zero();
    ActsMatrix<3, 3> dk2dT = ActsMatrix<3, 3>::Identity();
    ActsMatrix<3, 3> dk3dT = ActsMatrix<3, 3>::Identity();
    ActsMatrix<3, 3> dk4dT = ActsMatrix<3, 3>::Identity();

    Vector3 dk1dL = Vector3::Zero();
    Vector3 dk2dL = Vector3::Zero();
    Vector3 dk3dL = Vector3::Zero();
    Vector3 dk4dL = Vector3::Zero();

    /// Propagation of derivatives of dLambda''dlambda at each sub-step
    std::array<double, 4> jdL{};

    // Evaluation of the rightmost column without the last term.
    jdL[0] = dLdl[0];
    dk1dL = dir.cross(sd.B_first);

    jdL[1] = dLdl[1] * (1. + half_h * jdL[0]);
    dk2dL = (1. + half_h * jdL[0]) * (dir + half_h * sd.k1).cross(sd.B_middle) +
            qop[1] * half_h * dk1dL.cross(sd.B_middle);

    jdL[2] = dLdl[2] * (1. + half_h * jdL[1]);
    dk3dL = (1. + half_h * jdL[1]) * (dir + half_h * sd.k2).cross(sd.B_middle) +
            qop[2] * half_h * dk2dL.cross(sd.B_middle);

    jdL[3] = dLdl[3] * (1. + h * jdL[2]);
    dk4dL = (1. + h * jdL[2]) * (dir + h * sd.k3).cross(sd.B_last) +
            qop[3] * h * dk3dL.cross(sd.B_last);

    dk1dT(0, 1) = sd.B_first.z();
    dk1dT(0, 2) = -sd.B_first.y();
    dk1dT(1, 0) = -sd.B_first.z();
    dk1dT(1, 2) = sd.B_first.x();
    dk1dT(2, 0) = sd.B_first.y();
    dk1dT(2, 1) = -sd.B_first.x();
    dk1dT *= qop[0];

    dk2dT += half_h * dk1dT;
    dk2dT = qop[1] * VectorHelpers::cross(dk2dT, sd.B_middle);

    dk3dT += half_h * dk2dT;
    dk3dT = qop[2] * VectorHelpers::cross(dk3dT, sd.B_middle);

    dk4dT += h * dk3dT;
    dk4dT = qop[3] * VectorHelpers::cross(dk4dT, sd.B_last);

    dFdT.setIdentity();
    dFdT += h / 6. * (dk1dT + dk2dT + dk3dT);
    dFdT *= h;

    dFdL = h * h / 6. * (dk1dL + dk2dL + dk3dL);

    dGdT += h / 6. * (dk1dT + 2. * (dk2dT + dk3dT) + dk4dT);

    dGdL = h / 6. * (dk1dL + 2. * (dk2dL + dk3dL) + dk4dL);

    // Evaluation of the dLambda''/dlambda term
    D(7, 7) += (h / 6.) * (jdL[0] + 2. * (jdL[1] + jdL[2]) + jdL[3]);

    // The following comment lines refer to the application of the time being
    // treated as a position. Since t and qop are treated independently for now,
    // this just serves as entry point for building their relation
    //~ double dtpp1dl = -state.options.mass * state.options.mass * qop[0] *
    //~ qop[0] *
    //~ (3. * g + qop[0] * dgdqop(energy[0], state.options.mass,
    //~ state.options.absPdgCode,
    //~ state.options.meanEnergyLoss));

    double dtp1dl = qop[0] * state.options.mass * state.options.mass /
                    std::hypot(1, qop[0] * state.options.mass);
    double qopNew = qop[0] + half_h * Lambdappi[0];

    //~ double dtpp2dl = -state.options.mass * state.options.mass * qopNew *
    //~ qopNew *
    //~ (3. * g * (1. + half_h * jdL[0]) +
    //~ qopNew * dgdqop(energy[1], state.options.mass,
    //~ state.options.absPdgCode,
    //~ state.options.meanEnergyLoss));

    double dtp2dl = qopNew * state.options.mass * state.options.mass /
                    std::hypot(1, qopNew * state.options.mass);
    qopNew = qop[0] + half_h * Lambdappi[1];

    //~ double dtpp3dl = -state.options.mass * state.options.mass * qopNew *
    //~ qopNew *
    //~ (3. * g * (1. + half_h * jdL[1]) +
    //~ qopNew * dgdqop(energy[2], state.options.mass,
    //~ state.options.absPdgCode,
    //~ state.options.meanEnergyLoss));

    double dtp3dl = qopNew * state.options.mass * state.options.mass /
                    std::hypot(1, qopNew * state.options.mass);
    qopNew = qop[0] + half_h * Lambdappi[2];
    double dtp4dl = qopNew * state.options.mass * state.options.mass /
                    std::hypot(1, qopNew * state.options.mass);

    //~ D(3, 7) = h * state.options.mass * state.options.mass * qop[0] /
    //~ std::hypot(1., state.options.mass * qop[0])
    //~ + h * h / 6. * (dtpp1dl + dtpp2dl + dtpp3dl);

    D(3, 7) = (h / 6.) * (dtp1dl + 2. * (dtp2dl + dtp3dl) + dtp4dl);
    return true;
  }

  /// @brief Initializer of all parameters related to a RKN4 step with energy
  /// loss of a particle in material
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @param [in] state Deliverer of configurations
  template <typename propagator_state_t>
  void initializeEnergyLoss(const propagator_state_t& state) {
    using std::hypot;
    energy[0] = hypot(initialMomentum, state.options.mass);
    // use unit length as thickness to compute the energy loss per unit length
    Acts::MaterialSlab slab(material, 1);
    // Use the same energy loss throughout the step.
    if (state.options.meanEnergyLoss) {
      g = -computeEnergyLossMean(slab, state.options.absPdgCode,
                                 state.options.mass,
                                 static_cast<double>(qop[0]));
    } else {
      // TODO using the unit path length is not quite right since the most
      //      probably energy loss is not independent from the path length.
      g = -computeEnergyLossMode(slab, state.options.absPdgCode,
                                 state.options.mass,
                                 static_cast<double>(qop[0]));
    }
    // Change of the momentum per path length
    // dPds = dPdE * dEds
    dPds[0] = g * energy[0] / initialMomentum;
    if (state.stepping.covTransport) {
      // Calculate the change of the energy loss per path length and
      // inverse momentum
      if (state.options.includeGgradient) {
        if (state.options.meanEnergyLoss) {
          dgdqopValue = deriveEnergyLossMeanQOverP(
              slab, state.options.absPdgCode, state.options.mass,
              static_cast<double>(qop[0]));
        } else {
          // TODO path length dependence; see above
          dgdqopValue = deriveEnergyLossModeQOverP(
              slab, state.options.absPdgCode, state.options.mass,
              static_cast<double>(qop[0]));
        }
      }
      // Calculate term for later error propagation
      dLdl[0] = (-qop[0] * qop[0] * g * energy[0] *
                     (3. - (initialMomentum * initialMomentum) /
                               (energy[0] * energy[0])) -
                 qop[0] * qop[0] * qop[0] * energy[0] * dgdqopValue);
    }
  }

  /// @brief Update of the kinematic parameters of the RKN4 sub-steps after
  /// initialization with energy loss of a particle in material
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @tparam stepper_t Type of the stepper
  /// @param [in] h Stepped distance of the sub-step (1-3)
  /// @param [in] state State of the stepper
  /// @param [in] i Index of the sub-step (1-3)
  template <typename stepper_state_t, typename stepper_t>
  void updateEnergyLoss(const double mass, const double h,
                        const stepper_state_t& state, const stepper_t& stepper,
                        const int i) {
    // Update parameters related to a changed momentum
    currentMomentum = initialMomentum + h * dPds[i - 1];
    using std::sqrt;
    energy[i] = sqrt(currentMomentum * currentMomentum + mass * mass);
    dPds[i] = g * energy[i] / currentMomentum;
    qop[i] = stepper.charge(state) / currentMomentum;
    // Calculate term for later error propagation
    if (state.covTransport) {
      dLdl[i] = (-qop[i] * qop[i] * g * energy[i] *
                     (3. - (currentMomentum * currentMomentum) /
                               (energy[i] * energy[i])) -
                 qop[i] * qop[i] * qop[i] * energy[i] * dgdqopValue);
    }
  }
};

}  // namespace detail

}  // namespace Acts
