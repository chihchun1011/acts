// This file is part of the Acts project.
//
// Copyright (C) 2021-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFitting/Chi2FitterError.hpp"
#include "Acts/TrackFitting/detail/VoidChi2Components.hpp"
// TODO: generalize voidKalmanCalibrator?
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <functional>
#include <map>
#include <memory>

namespace Acts {

namespace Experimental {

/// Extension struct which holds delegates to customize the GX2F behavior
template <typename traj_t>
struct Chi2FitterExtensions {
  using TrackStateProxy = typename MultiTrajectory<traj_t>::TrackStateProxy;
  using ConstTrackStateProxy =
      typename MultiTrajectory<traj_t>::ConstTrackStateProxy;
  using Parameters = typename TrackStateProxy::Parameters;

  using Calibrator = Delegate<void(const GeometryContext&, TrackStateProxy)>;
  using OutlierFinder = Delegate<bool(ConstTrackStateProxy)>;

  /// The Calibrator is a dedicated calibration algorithm that allows
  /// to calibrate measurements using track information, this could be
  /// e.g. sagging for wires, module deformations, etc.
  Calibrator calibrator;

  /// Determines whether a measurement is supposed to be considered as an
  /// outlier
  OutlierFinder outlierFinder;

  /// Default constructor which connects the default void components
  Chi2FitterExtensions() {
    calibrator.template connect<&voidChi2Calibrator<traj_t>>();
    outlierFinder.template connect<&voidChi2OutlierFinder<traj_t>>();
  }
};

/// Combined options for the GX2F fitter.
/// @tparam traj_t The trajectory type
template <typename traj_t>
struct Chi2FitterOptions {
  /// PropagatorOptions with context.
  ///
  /// @param gctx The goemetry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param cctx The calibration context for this fit
  /// @param extensions_ The chi2 extensions
  /// @param logger_ The logger wrapper
  /// @param pOptions The plain propagator options
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param nIter Number of update steps to the parameters
  /// @param calcFinalChi2_ Whether to run additional propagation to calculate
  /// final chi2
  /// @param freeToBoundCorrection_ Correction for non-linearity effect during
  /// transform from free to bound
  Chi2FitterOptions(const GeometryContext& gctx,
                    const MagneticFieldContext& mctx,
                    std::reference_wrapper<const CalibrationContext> cctx,
                    Chi2FitterExtensions<traj_t> extensions_,
                    LoggerWrapper logger_,
                    const PropagatorPlainOptions& pOptions,
                    bool mScattering = false, bool eLoss = false, int nIter = 1,
                    bool calcFinalChi2_ = true,
                    const FreeToBoundCorrection& freeToBoundCorrection_ =
                        FreeToBoundCorrection(false))
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        extensions(std::move(extensions_)),
        propagatorPlainOptions(pOptions),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        nUpdates(nIter),
        calcFinalChi2(calcFinalChi2_),
        freeToBoundCorrection(freeToBoundCorrection_),
        logger(logger_) {}
  /// Contexts are required and the options must not be default-constructible.
  Chi2FitterOptions() = delete;

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  Chi2FitterExtensions<traj_t> extensions;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// Whether to consider multiple scattering
  bool multipleScattering = false;  // TODO: add later

  /// Whether to consider energy loss
  bool energyLoss = false;  // TODO: add later

  /// Number of iterations to improve chi2
  int nUpdates = 1;

  /// Whether to do an additional propagation step, just to get the latest chi2
  /// value
  bool calcFinalChi2 = true;

  /// Whether to include non-linear correction during global to local
  /// transformation
  FreeToBoundCorrection freeToBoundCorrection;

  /// Logger
  LoggerWrapper logger;
};

template <typename traj_t>
struct Chi2FitterResult {
  // Fitted states that the actor has handled.
  std::shared_ptr<traj_t> fittedStates;

  // This is the index of the 'tip' of the track stored in multitrajectory.
  // This correspond to the last measurement state in the multitrajectory.
  // Since this GX2F only stores one trajectory, it is unambiguous.
  // SIZE_MAX is the start of a trajectory.
  size_t lastMeasurementIndex = SIZE_MAX;

  // This is the index of the 'tip' of the states stored in multitrajectory.
  // This correspond to the last state in the multitrajectory.
  // Since this GX2F only stores one trajectory, it is unambiguous.
  // SIZE_MAX is the start of a trajectory.
  size_t lastTrackIndex = SIZE_MAX;

  // The optional Parameters at the provided surface
  std::optional<BoundTrackParameters> fittedParameters;

  // Counter for states with non-outlier measurements
  size_t measurementStates = 0;

  // Counter for measurements holes
  // A hole correspond to a surface with an associated detector element with no
  // associated measurement. Holes are only taken into account if they are
  // between the first and last measurements.
  size_t measurementHoles = 0;

  // Counter for handled states
  size_t processedStates = 0;

  // Indicator if track fitting has been done
  bool finished = false;

  // Measurement surfaces without hits
  std::vector<const Surface*> missedActiveSurfaces;

  // collectors
  std::vector<ActsScalar> collectorMeasurements;
  std::vector<ActsScalar> collectorCovariance;
  std::vector<ActsScalar> collectorResiduals;

  /// first derivative of chi2 wrt starting track parameters
  BoundVector collectorDerive1Chi2Sum = BoundVector::Zero();
  BoundMatrix collectorDerive2Chi2Sum = BoundMatrix::Zero();

  BoundMatrix jacobianFromStart = BoundMatrix::Identity();

  // chi2 fitter results
  ActsDynamicVector residuals;
  ActsDynamicMatrix covariance;
  ActsScalar chisquare = -1;
  std::vector<ActsScalar> chisquares;

  Result<void> result{Result<void>::success()};
};

/// Chi2 fitter implementation.
///
/// @tparam propagator_t Type of the propagation class
template <typename propagator_t, typename traj_t>
class Chi2Fitter {
  using Chi2Navigator = typename propagator_t::Navigator;

 public:
  Chi2Fitter(propagator_t pPropagator) : m_propagator(std::move(pPropagator)) {}

 private:
  /// The propgator for the transport and material update
  propagator_t m_propagator;

  /// @brief Propagator Actor plugin for the Chi2Fitter
  ///
  /// @tparam parameters_t The type of parameters used for "local" parameters.
  ///
  /// The Chi2Actor does not rely on the measurements to be
  /// sorted along the track.
  template <typename parameters_t>
  class Actor {
   public:
    /// Broadcast the result_type
    using result_type = Chi2FitterResult<traj_t>;

    /// Allows retrieving measurements for a surface
    const std::map<GeometryIdentifier,
                   std::reference_wrapper<const SourceLink>>*
        inputMeasurements = nullptr;

    /// Whether to consider multiple scattering.
    bool multipleScattering = false;  // TODO: add later

    /// Whether to consider energy loss.
    bool energyLoss = false;  // TODO: add later

    /// Whether to include non-linear correction during global to local
    /// transformation
    FreeToBoundCorrection freeToBoundCorrection;

    /// Extension struct
    Chi2FitterExtensions<traj_t> extensions;

    /// @brief Chi square actor operation
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    void operator()(propagator_state_t& state, const stepper_t& stepper,
                    result_type& result) const {
      const auto& logger = state.options.logger;

      if (result.finished) {
        return;
      }

      // Add the measurement surface as external surface to navigator.
      // We will try to hit those surface by ignoring boundary checks.
      if (result.processedStates == 0) {
        for (auto measurementIt = inputMeasurements->begin();
             measurementIt != inputMeasurements->end(); measurementIt++) {
          state.navigation.externalSurfaces.insert(
              std::pair<uint64_t, GeometryIdentifier>(
                  measurementIt->first.layer(), measurementIt->first));
        }
      }

      // wait for surface
      auto surface = state.navigation.currentSurface;
      if (surface != nullptr) {
        auto res = processSurface(surface, state, stepper, result);
        if (!res.ok()) {
          ACTS_ERROR("chi2 | Error in processSurface: " << res.error());
          result.result = res.error();
        }
      }

      // finalization
      if (not result.finished) {
        if (result.measurementStates == inputMeasurements->size() or
            (result.measurementStates > 0 and
             state.navigation.navigationBreak)) {
          result.missedActiveSurfaces.resize(result.measurementHoles);
          ACTS_VERBOSE("chi2 | Finalize...");
          result.finished = true;
        }
      }
    }

    /// @brief Chi2 actor operation: process surface
    ///
    /// @tparam propagator_state_t is the type of Propagator state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result The mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    Result<void> processSurface(const Surface* surface,
                                propagator_state_t& state,
                                const stepper_t& stepper,
                                result_type& result) const {
      const auto& logger = state.options.logger;

      // We need the full jacobianFromStart, so we'll need to calculate it no
      // matter if we have a measurement or not.

      // Transport the covariance to the surface
      stepper.transportCovarianceToBound(state.stepping, *surface,
                                         freeToBoundCorrection);

      // Update state and stepper with pre material effects
      materialInteractor(surface, state, stepper,
                         MaterialUpdateStage::PreUpdate);
      // TODO: do we need the materialInteractor before we access the
      // boundState? In the material-only case in the KF, the materialInteractor
      // is called with fullUpdate *after* retrieving the boundState.

      // Bind the transported state to the current surface
      auto res = stepper.boundState(state.stepping, *surface, false);
      if (!res.ok()) {
        return res.error();
      }
      auto& [boundParams, jacobian, pathLength] = *res;
      // jacobian is "the stepwise jacobian towards the bound state (from last
      // bound)"

      result.jacobianFromStart = jacobian * result.jacobianFromStart;

      // Try to find the surface in all measurement surfaces
      auto sourcelink_it = inputMeasurements->find(surface->geometryId());
      // inputMeasurements is a std::map<GeometryIdentifier, source_link_t>
      if (sourcelink_it != inputMeasurements->end()) {
        ACTS_VERBOSE("chi2 |    processSurface: Measurement surface "
                     << surface->geometryId() << " detected.");

        // add a full TrackState entry multi trajectory
        // (this allocates storage for all components, we will set them later)
        result.lastTrackIndex = result.fittedStates->addTrackState(
            ~(TrackStatePropMask::Smoothed | TrackStatePropMask::Filtered),
            result.lastTrackIndex);

        // now get track state proxy back
        auto trackStateProxy =
            result.fittedStates->getTrackState(result.lastTrackIndex);

        trackStateProxy.setReferenceSurface(surface->getSharedPtr());

        // assign the source link to the track state
        trackStateProxy.setUncalibrated(sourcelink_it->second);

        // Fill the track state
        trackStateProxy.predicted() = std::move(boundParams.parameters());
        if (boundParams.covariance().has_value()) {
          trackStateProxy.predictedCovariance() =
              std::move(*boundParams.covariance());
        }
        trackStateProxy.jacobian() = std::move(jacobian);
        trackStateProxy.pathLength() = std::move(pathLength);

        extensions.calibrator(state.geoContext, trackStateProxy);

        visit_measurement(trackStateProxy.calibratedSize(), [&](auto N) {
          constexpr size_t kMeasurementSize = decltype(N)::value;

          // simple projection matrix H_is, composed of 1 and 0, 2x6 or 1x6
          const ActsMatrix<kMeasurementSize, eBoundSize> proj =
              trackStateProxy.projector()
                  .template topLeftCorner<kMeasurementSize, eBoundSize>();

          const auto Hi =
              (proj * result.jacobianFromStart).eval();  // 2x6 or 1x6
          const auto localMeasurements =
              trackStateProxy
                  .template calibrated<kMeasurementSize>();  // 2x1 or 1x1

          const auto& covariance =
              trackStateProxy.template calibratedCovariance<
                  kMeasurementSize>();  // 2x2 or 1x1. Should
                                        // be diagonal.
          const auto covInv = covariance.inverse();

          auto residuals =
              localMeasurements - proj * trackStateProxy.predicted();

          // TODO: use detail::calculateResiduals? Theta/Phi?
          const auto derive1Chi2 =
              (-2 * Hi.transpose() * covInv * residuals).eval();
          const auto derive2Chi2 = (2 * Hi.transpose() * covInv * Hi).eval();
          result.collectorDerive1Chi2Sum += derive1Chi2;
          result.collectorDerive2Chi2Sum += derive2Chi2;

          for (int i = 0; i < localMeasurements.rows(); ++i) {
            result.collectorMeasurements.push_back(localMeasurements(i));
            result.collectorResiduals.push_back(residuals(i));
            result.collectorCovariance.push_back(covariance(i, i));
            // we assume measurements are not correlated
          }
        });

        //====================================

        // Get and set the type flags
        auto& typeFlags = trackStateProxy.typeFlags();
        typeFlags.set(TrackStateFlag::ParameterFlag);
        if (surface->surfaceMaterial() != nullptr) {
          typeFlags.set(TrackStateFlag::MaterialFlag);
        }

        if (not extensions.outlierFinder(trackStateProxy)) {
          typeFlags.set(TrackStateFlag::MeasurementFlag);
          ++result.measurementStates;
        } else {
          ACTS_VERBOSE("chi2 | Measurement is determined to be an outlier.");
          typeFlags.set(TrackStateFlag::OutlierFlag);
        }

        // We count the processed states
        ++result.processedStates;
        // Update the number of holes only when encoutering a measurement
        result.measurementHoles = result.missedActiveSurfaces.size();
        // Since we encountered a measurement update the lastMeasurementIndex to
        // the lastTrackIndex.
        result.lastMeasurementIndex = result.lastTrackIndex;

      } else if (surface->associatedDetectorElement() != nullptr ||
                 surface->surfaceMaterial() != nullptr) {
        // We only create track states here if there is already measurement
        // detected or if the surface has material (no holes before the first
        // measurement
        if (result.measurementStates > 0 ||
            surface->surfaceMaterial() != nullptr) {
          // No source links on surface, add either hole or passive material
          // TrackState entry multi trajectory. No storage allocation for
          // uncalibrated/calibrated measurement and filtered parameter
          result.lastTrackIndex = result.fittedStates->addTrackState(
              ~(TrackStatePropMask::Calibrated | TrackStatePropMask::Filtered |
                TrackStatePropMask::Smoothed),
              result.lastTrackIndex);

          // now get track state proxy back
          auto trackStateProxy =
              result.fittedStates->getTrackState(result.lastTrackIndex);

          // Set the surface
          trackStateProxy.setReferenceSurface(surface->getSharedPtr());

          // Set the track state flags
          auto& typeFlags = trackStateProxy.typeFlags();
          typeFlags.set(TrackStateFlag::ParameterFlag);
          if (surface->surfaceMaterial() != nullptr) {
            typeFlags.set(TrackStateFlag::MaterialFlag);
          }
          if (surface->associatedDetectorElement() != nullptr) {
            ACTS_VERBOSE("Detected hole on " << surface->geometryId());
            // If the surface is sensitive, set the hole type flag
            typeFlags.set(TrackStateFlag::HoleFlag);

            // Count the missed surface
            result.missedActiveSurfaces.push_back(surface);
          } else if (surface->surfaceMaterial() != nullptr) {
            ACTS_VERBOSE("Detected in-sensitive surface "
                         << surface->geometryId());
          }

          // note: the track state is already transported/bound to the surface.

          // Fill the track state
          trackStateProxy.predicted() = std::move(boundParams.parameters());
          if (boundParams.covariance().has_value()) {
            trackStateProxy.predictedCovariance() =
                std::move(*boundParams.covariance());
          }
          trackStateProxy.jacobian() = std::move(jacobian);
          trackStateProxy.pathLength() = std::move(pathLength);

          // We count the processed state
          ++result.processedStates;
        }
      }

      // Update state and stepper with post material effects
      materialInteractor(surface, state, stepper,
                         MaterialUpdateStage::PostUpdate);
      // TODO: is it more expensive to split the materialInteractor into
      // preUpdate and postUpdate vs fullUpdate? One could just fullUpdate in
      // case we have no measurement, like the KF?

      return Result<void>::success();
    }

    /// @brief Chi2Fitter actor operation : material interaction
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param surface The surface where the material interaction happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param updateStage The materal update stage
    ///
    template <typename propagator_state_t, typename stepper_t>
    void materialInteractor(const Surface* surface, propagator_state_t& state,
                            stepper_t& stepper,
                            const MaterialUpdateStage& updateStage =
                                MaterialUpdateStage::FullUpdate) const {
      const auto& logger = state.options.logger;
      // Indicator if having material
      bool hasMaterial = false;

      if (surface and surface->surfaceMaterial()) {
        // Prepare relevant input particle properties
        detail::PointwiseMaterialInteraction interaction(surface, state,
                                                         stepper);
        // Evaluate the material properties
        if (interaction.evaluateMaterialSlab(state, updateStage)) {
          // Surface has material at this stage
          hasMaterial = true;

          // Evaluate the material effects
          interaction.evaluatePointwiseMaterialInteraction(multipleScattering,
                                                           energyLoss);

          // Screen out material effects info
          ACTS_VERBOSE("Material effects on surface: "
                       << surface->geometryId()
                       << " at update stage: " << updateStage << " are :");
          ACTS_VERBOSE("eLoss = "
                       << interaction.Eloss << ", "
                       << "variancePhi = " << interaction.variancePhi << ", "
                       << "varianceTheta = " << interaction.varianceTheta
                       << ", "
                       << "varianceQoverP = " << interaction.varianceQoverP);

          // Update the state and stepper with material effects
          interaction.updateState(state, stepper);
        }
      }

      if (not hasMaterial) {
        // Screen out message
        ACTS_VERBOSE("No material effects on surface: " << surface->geometryId()
                                                        << " at update stage: "
                                                        << updateStage);
      }
    }
  };

  template <typename parameters_t>
  class Aborter {
   public:
    /// Broadcast the action_type
    using action_type = Actor<parameters_t>;

    template <typename propagator_state_t, typename stepper_t,
              typename result_t>
    bool operator()(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const result_t& result) const {
      // const auto& logger = state.options.logger;
      if (!result.result.ok() or result.finished) {
        return true;
      }
      return false;
    }
  };

 public:
  /// Fit implementation of the GX2F
  ///
  /// @tparam source_link_iterator_t Iterator type used to pass source links
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param it Begin iterator for the fittable uncalibrated measurements
  /// @param end End iterator for the fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param chi2FitterOptions Chi2FitterOptions steering the fit
  /// @param trajectory Input trajectory storage to append into
  /// @note The input measurements are given in the form of @c SourceLink s.
  /// It's the calibrators job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_iterator_t, typename start_parameters_t,
            typename parameters_t = BoundTrackParameters>
  Result<Chi2FitterResult<traj_t>> fit(
      source_link_iterator_t it, source_link_iterator_t end,
      const start_parameters_t& sParameters,
      const Chi2FitterOptions<traj_t>& chi2FitterOptions,
      std::shared_ptr<traj_t> trajectory = {}) const {
    const auto& logger = chi2FitterOptions.logger;

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_VERBOSE("chi2 | preparing " << std::distance(it, end)
                                     << " input measurements");
    std::map<GeometryIdentifier, std::reference_wrapper<const SourceLink>>
        inputMeasurements;

    for (; it != end; ++it) {
      const SourceLink& sl = *it;
      inputMeasurements.emplace(sl.geometryId(), sl);
    }

    // TODO: for now, we use STL objects to collect the information during
    // propagation. Use dynamic Matrix instead? Performance?

    // Create the ActionList and AbortList
    using Chi2Aborter = Aborter<parameters_t>;
    using Chi2Actor = Actor<parameters_t>;

    using Chi2Result = typename Chi2Actor::result_type;
    using Actors = ActionList<Chi2Actor>;
    using Aborters = AbortList<Chi2Aborter>;

    // the result object which will be returned. Overridden every iteration.
    Chi2Result c2r;

    for (int i = 0; i <= chi2FitterOptions.nUpdates; ++i) {
      // Create relevant options for the propagation options
      PropagatorOptions<Actors, Aborters> propOptions(
          chi2FitterOptions.geoContext, chi2FitterOptions.magFieldContext,
          logger);

      // Set the trivial propagator options
      propOptions.setPlainOptions(chi2FitterOptions.propagatorPlainOptions);

      // Catch the actor and set the measurements
      auto& chi2Actor = propOptions.actionList.template get<Chi2Actor>();
      chi2Actor.inputMeasurements = &inputMeasurements;
      chi2Actor.multipleScattering = chi2FitterOptions.multipleScattering;
      chi2Actor.energyLoss = chi2FitterOptions.energyLoss;
      chi2Actor.freeToBoundCorrection = chi2FitterOptions.freeToBoundCorrection;
      chi2Actor.extensions = chi2FitterOptions.extensions;

      typename propagator_t::template action_list_t_result_t<
          CurvilinearTrackParameters, Actors>
          inputResult;

      auto& r = inputResult.template get<Chi2FitterResult<traj_t>>();

      if (trajectory) {
        r.fittedStates = trajectory;
      } else {
        r.fittedStates = std::make_shared<traj_t>();
      }

      using paramType = typename std::conditional<
          std::is_same<start_parameters_t, parameters_t>::value,
          std::variant<start_parameters_t>,
          std::variant<start_parameters_t, parameters_t>>::type;
      paramType vParams = sParameters;
      // start_parameters_t and parameter_t can be the same

      auto result = m_propagator.template propagate(sParameters, propOptions,
                                                    std::move(inputResult));

      if (!result.ok()) {
        ACTS_ERROR("chi2 | it=" << i
                                << " | propapation failed: " << result.error());
        return result.error();
      }

      Chi2Result c2rCurrent =
          std::move(result.value().template get<Chi2Result>());

      /// It could happen that the fit ends in zero measurement states.
      /// The result gets meaningless so such case is regarded as fit failure.
      if (c2rCurrent.result.ok() and not c2rCurrent.measurementStates) {
        c2rCurrent.result = Result<void>(Chi2FitterError::NoMeasurementFound);
      }

      if (!c2rCurrent.result.ok()) {
        ACTS_ERROR("chi2 | it=" << i << " | Chi2Fitter failed: "
                                << c2rCurrent.result.error() << ", "
                                << c2rCurrent.result.error().message());
        return c2rCurrent.result.error();
      }

      c2rCurrent.residuals =
          Eigen::Map<ActsDynamicVector>(c2rCurrent.collectorResiduals.data(),
                                        c2rCurrent.collectorResiduals.size());

      ActsDynamicVector variance =
          Eigen::Map<ActsDynamicVector>(c2rCurrent.collectorCovariance.data(),
                                        c2rCurrent.collectorCovariance.size());
      c2rCurrent.covariance = variance.asDiagonal();

      c2rCurrent.chisquare = c2rCurrent.residuals.transpose() *
                             c2rCurrent.covariance.inverse() *
                             c2rCurrent.residuals;
      ACTS_VERBOSE("chi2 | it=" << i << " | χ² = " << c2rCurrent.chisquare);

      // copy over data from previous runs (namely chisquares vector)
      c2rCurrent.chisquares.reserve(c2r.chisquares.size() + 1);
      c2rCurrent.chisquares.insert(c2rCurrent.chisquares.end(),
                                   c2r.chisquares.begin(),
                                   c2r.chisquares.end());
      c2rCurrent.chisquares.push_back(c2rCurrent.chisquare);
      // TODO: is there a more elegant/optimal way to do this?

      c2r = std::move(c2rCurrent);

      if (i == chi2FitterOptions.nUpdates) {
        // don't update parameters in last iteration
        c2r.fittedParameters = std::visit(
            [](auto&& prevParams) {
              return BoundTrackParameters(
                  prevParams.referenceSurface().getSharedPtr(),
                  prevParams.parameters(), prevParams.covariance());
            },
            vParams);
        break;

        // TODO: verify if another step would be useful, e.g. by comparing the
        // delta to a desired minimum value
      }

      // calculate updates to parameters
      BoundVector delta_start_parameters =
          c2r.collectorDerive2Chi2Sum.colPivHouseholderQr().solve(
              c2r.collectorDerive1Chi2Sum);

      c2r.fittedParameters = std::visit(
          [delta_start_parameters, logger, i](auto&& prevParams) {
            BoundVector newParamsVec =
                prevParams.parameters() - delta_start_parameters;
            ACTS_VERBOSE("chi2 | it=" << i << " | updated parameters = "
                                      << newParamsVec.transpose());

            return BoundTrackParameters(
                prevParams.referenceSurface().getSharedPtr(), newParamsVec,
                prevParams.covariance());
          },
          vParams);

      vParams = c2r.fittedParameters.value();  // passed to next iteration
    }

    // Return the converted track
    return c2r;
  }
};

}  // namespace Experimental

}  // namespace Acts
