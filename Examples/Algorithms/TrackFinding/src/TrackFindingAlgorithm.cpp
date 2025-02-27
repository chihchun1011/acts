// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <memory>
#include <stdexcept>

#include <boost/histogram.hpp>

ActsExamples::TrackFindingAlgorithm::TrackFindingAlgorithm(
    Config config, Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("TrackFindingAlgorithm", level),
      m_cfg(std::move(config)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements input collection");
  }
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links input collection");
  }
  if (m_cfg.inputInitialTrackParameters.empty()) {
    throw std::invalid_argument(
        "Missing initial track parameters input collection");
  }
  if (m_cfg.outputTrajectories.empty()) {
    throw std::invalid_argument("Missing trajectories output collection");
  }
}

ActsExamples::ProcessCode ActsExamples::TrackFindingAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  const auto& measurements =
      ctx.eventStore.get<MeasurementContainer>(m_cfg.inputMeasurements);
  const auto& sourceLinks =
      ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);
  const auto& initialParameters = ctx.eventStore.get<TrackParametersContainer>(
      m_cfg.inputInitialTrackParameters);

  // Prepare the output data with MultiTrajectory
  TrajectoriesContainer trajectories;
  trajectories.reserve(initialParameters.size());

  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  Acts::PropagatorPlainOptions pOptions;
  pOptions.maxSteps = 10000;

  MeasurementCalibrator calibrator{measurements};
  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;
  Acts::MeasurementSelector measSel{m_cfg.measurementSelectorCfg};

  Acts::CombinatorialKalmanFilterExtensions<Acts::VectorMultiTrajectory>
      extensions;
  extensions.calibrator.connect<&MeasurementCalibrator::calibrate>(&calibrator);
  extensions.updater.connect<
      &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
      &kfUpdater);
  extensions.smoother.connect<
      &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
      &kfSmoother);
  extensions.measurementSelector
      .connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(
          &measSel);

  IndexSourceLinkAccessor slAccessor;
  slAccessor.container = &sourceLinks;
  Acts::SourceLinkAccessorDelegate<IndexSourceLinkAccessor::Iterator>
      slAccessorDelegate;
  slAccessorDelegate.connect<&IndexSourceLinkAccessor::range>(&slAccessor);

  // Set the CombinatorialKalmanFilter options
  ActsExamples::TrackFindingAlgorithm::TrackFinderOptions options(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext, slAccessorDelegate,
      extensions, Acts::LoggerWrapper{logger()}, pOptions, &(*pSurface));

  // Perform the track finding for all initial parameters
  ACTS_DEBUG("Invoke track finding with " << initialParameters.size()
                                          << " seeds.");

  auto mtj = std::make_shared<Acts::VectorMultiTrajectory>();
  auto results = (*m_cfg.findTracks)(initialParameters, options, mtj);

  // Compute shared hits from all the reconstructed tracks
  if (m_cfg.computeSharedHits) {
    computeSharedHits(sourceLinks, results);
  }

  // Loop over the track finding results for all initial parameters
  for (std::size_t iseed = 0; iseed < initialParameters.size(); ++iseed) {
    m_nTotalSeeds++;
    // The result for this seed
    auto& result = results[iseed];
    if (result.ok()) {
      // Get the track finding output object
      auto& trackFindingOutput = result.value();
      // Create a Trajectories result struct
      trajectories.emplace_back(trackFindingOutput.fittedStates,
                                trackFindingOutput.lastMeasurementIndices,
                                trackFindingOutput.fittedParameters);
    } else {
      ACTS_WARNING("Track finding failed for seed " << iseed << " with error"
                                                    << result.error());
      m_nFailedSeeds++;
      // Track finding failed. Add an empty result so the output container has
      // the same number of entries as the input.
      trajectories.push_back(Trajectories());
    }
  }

  ACTS_DEBUG("Finalized track finding with " << trajectories.size()
                                             << " track candidates.");

  m_memoryStatistics.local().hist += mtj->statistics().hist;

  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::TrackFindingAlgorithm::finalize()
    const {
  ACTS_INFO("TrackFindingAlgorithm statistics:");
  ACTS_INFO("- total seeds: " << m_nTotalSeeds);
  ACTS_INFO("- failed seeds: " << m_nFailedSeeds);
  ACTS_INFO("- failure ratio: " << static_cast<double>(m_nFailedSeeds) /
                                       m_nTotalSeeds);

  auto memoryStatistics =
      m_memoryStatistics.combine([](const auto& a, const auto& b) {
        Acts::VectorMultiTrajectory::Statistics c;
        c.hist = a.hist + b.hist;
        return c;
      });
  std::stringstream ss;
  memoryStatistics.toStream(ss);
  ACTS_DEBUG("Track State memory statistics (averaged):\n" << ss.str());
  return ProcessCode::SUCCESS;
}
