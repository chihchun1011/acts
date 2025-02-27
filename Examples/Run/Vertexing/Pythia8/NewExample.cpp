// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/MagneticFieldOptions.hpp"
#include "ActsExamples/Printers/ParticlesPrinter.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Options/ParticleSelectorOptions.hpp"
#include "ActsExamples/Options/ParticleSmearingOptions.hpp"
#include "ActsExamples/Options/Pythia8Options.hpp"
#include "ActsExamples/Options/VertexingOptions.hpp"
#include "ActsExamples/Reconstruction/ReconstructionBase.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/Printers/TrackParametersPrinter.hpp"
#include "ActsExamples/Io/Root/RootSmearedTrackParameterWriter.hpp"
#include "ActsExamples/Vertexing/AdaptiveMultiVertexFinderAlgorithm.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addPythia8Options(desc);
  Options::addParticleSelectorOptions(desc);
  Options::addVertexingOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  Options::addParticleSmearingOptions(desc);
  auto vars = Options::parse(desc, argc, argv);
  if (vars.empty()) {
    return EXIT_FAILURE;
  }

  // basic setup
  auto logLevel = Options::readLogLevel(vars);
  auto outputDir = ensureWritableDirectory(vars["output-dir"].as<std::string>());
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));
  Sequencer sequencer(Options::readSequencerConfig(vars));

  // Setup the magnetic field
  auto magneticField = Options::readMagneticField(vars);

  // setup event generator
  EventGenerator::Config evgen = Options::readPythia8Options(vars, logLevel);
  evgen.outputParticles = "particles_generated";
  evgen.randomNumbers = rnd;
  sequencer.addReader(std::make_shared<EventGenerator>(evgen, logLevel));

  // print particle_init
  // ParticlesPrinter::Config printInitialCfg;
  // printInitialCfg.inputParticles = evgen.outputParticles;
  // sequencer.addAlgorithm(
  //     std::make_shared<ParticlesPrinter>(printInitialCfg, logLevel));

  // Output particle_init
  RootParticleWriter::Config writeInitial;
  writeInitial.inputParticles = evgen.outputParticles;
  writeInitial.filePath = joinPaths(outputDir, "particles_initial.root");
  writeInitial.treeName = "particles";
  // std::cout<<writeInitial.treeName<<std::endl;
  sequencer.addWriter(
      std::make_shared<RootParticleWriter>(writeInitial, logLevel));

  // pre-select particles
  ParticleSelector::Config selectParticles =
      Options::readParticleSelectorConfig(vars);
  selectParticles.inputParticles = evgen.outputParticles;
  selectParticles.outputParticles = "particles_selected";
  // smearing only works with charge particles for now
  selectParticles.removeNeutral = true;
  selectParticles.absEtaMax = vars["vertexing-eta-max"].as<double>();
  selectParticles.rhoMax = vars["vertexing-rho-max"].as<double>() * 1_mm;
  selectParticles.ptMin = vars["vertexing-pt-min"].as<double>() * 1_MeV;
  sequencer.addAlgorithm(
      std::make_shared<ParticleSelector>(selectParticles, logLevel));
  // Run the particle smearing
  auto particleSmearingCfg = setupParticleSmearing(
      vars, sequencer, rnd, selectParticles.outputParticles);

  // print track parameter
  // TrackParametersPrinter::Config printTracks;
  // printTracks.inputTrackParameters = particleSmearingCfg.outputTrackParameters;
  // sequencer.addAlgorithm(
  // std::make_shared<TrackParametersPrinter>(printTracks, logLevel));

  //track parameter writer
  RootSmearedTrackParameterWriter::Config trackParamsWriterCfg;
  trackParamsWriterCfg.inputTrackParameters = particleSmearingCfg.outputTrackParameters;
  trackParamsWriterCfg.filePath = joinPaths(outputDir, "/smearedparams.root");
  trackParamsWriterCfg.treeName = "smearedparams";
  sequencer.addWriter(std::make_shared<RootSmearedTrackParameterWriter>(
      trackParamsWriterCfg, logLevel));

  // find vertices
  AdaptiveMultiVertexFinderAlgorithm::Config findVertices;
  findVertices.bField = magneticField;
  findVertices.inputTrackParameters = particleSmearingCfg.outputTrackParameters;
  findVertices.outputProtoVertices = "protovertices";
  sequencer.addAlgorithm(std::make_shared<AdaptiveMultiVertexFinderAlgorithm>(
      findVertices, logLevel));

  return sequencer.run();
}
