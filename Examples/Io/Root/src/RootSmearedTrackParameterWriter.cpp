// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootSmearedTrackParameterWriter.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <ios>
#include <iostream>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::RootSmearedTrackParameterWriter::RootSmearedTrackParameterWriter(
    const ActsExamples::RootSmearedTrackParameterWriter::Config& config,
    Acts::Logging::Level level)
    : TrackParameterWriter(config.inputTrackParameters,
                           "RootSmearedTrackParameterWriter", level),
      m_cfg(config) {
  
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    auto path = m_cfg.filePath;
    m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + path);
    }
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  } else {
    // The estimated track parameters
    m_outputTree->Branch("event_nr", &m_eventNr);
    m_outputTree->Branch("loc0", &m_loc0);
    m_outputTree->Branch("loc1", &m_loc1);
    m_outputTree->Branch("phi", &m_phi);
    m_outputTree->Branch("theta", &m_theta);
    m_outputTree->Branch("qop", &m_qop);
    m_outputTree->Branch("time", &m_time);
    m_outputTree->Branch("p", &m_p);
    m_outputTree->Branch("pt", &m_pt);
    m_outputTree->Branch("eta", &m_eta);
  }
}

ActsExamples::RootSmearedTrackParameterWriter::~RootSmearedTrackParameterWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootSmearedTrackParameterWriter::endRun() {
  if (m_outputFile != nullptr) {
    m_outputFile->cd();
    m_outputTree->Write();
    ACTS_INFO("Write estimated parameters from seed to tree '"
              << m_cfg.treeName << "' in '" << m_cfg.filePath << "'");
    m_outputFile->Close();
  }
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootSmearedTrackParameterWriter::writeT(
    const ActsExamples::AlgorithmContext& ctx,
    const TrackParametersContainer& trackParams) {

  if (m_outputFile == nullptr) {
    return ProcessCode::SUCCESS;
  }

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;

  ACTS_VERBOSE("Writing " << trackParams.size() << " track parameters");

  // Loop over the estimated track parameters
  for (size_t iparams = 0; iparams < trackParams.size(); ++iparams) {
    // The reference surface of the parameters, i.e. also the reference surface
    // of the first space point
    
    // The smeared parameters vector
    const auto params = trackParams[iparams].parameters();
    m_loc0 = params[Acts::eBoundLoc0];
    m_loc1 = params[Acts::eBoundLoc1];
    m_phi = params[Acts::eBoundPhi];
    m_theta = params[Acts::eBoundTheta];
    m_qop = params[Acts::eBoundQOverP];
    m_time = params[Acts::eBoundTime];
    m_p = std::abs(1.0 / m_qop);
    m_pt = m_p * std::sin(m_theta);
    m_eta = std::atanh(std::cos(m_theta));

    m_outputTree->Fill();
  }

  return ProcessCode::SUCCESS;
}