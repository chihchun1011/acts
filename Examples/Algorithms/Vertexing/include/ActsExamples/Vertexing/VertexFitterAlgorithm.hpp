// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <string>

namespace ActsExamples {

class VertexFitterAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Optional. Input track parameters collection
    std::string inputTrackParameters;
    /// Optional. Input trajectories container.
    std::string inputTrajectories;
    /// Input proto vertex collection
    std::string inputProtoVertices;
    /// Output vertex collection
    std::string outputVertices;
    /// The magnetic field
    std::shared_ptr<Acts::MagneticFieldProvider> bField;
    /// Constraint vertex fit bool
    bool doConstrainedFit = false;
    /// Vertex constraint position
    Acts::Vector4 constraintPos = Acts::Vector4(0, 0, 0, 0);
    /// Vertex constraint covariance matrix
    Acts::SymMatrix4 constraintCov =
        Acts::Vector4(3 * Acts::UnitConstants::mm * Acts::UnitConstants::mm,
                      3 * Acts::UnitConstants::mm * Acts::UnitConstants::mm,
                      10 * Acts::UnitConstants::mm * Acts::UnitConstants::mm,
                      1 * Acts::UnitConstants::ns * Acts::UnitConstants::ns)
            .asDiagonal();
  };

  VertexFitterAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  /// Fit the input vertices.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
