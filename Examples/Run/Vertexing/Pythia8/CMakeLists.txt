set(_common_libraries
  ActsCore
  ActsExamplesFramework
  ActsExamplesCommon
  ActsExamplesGenerators
  ActsExamplesGeneratorsPythia8
  ActsExamplesGeneratorsPythia8Options
  ActsExamplesTruthTracking
  ActsExamplesVertexing
  Boost::program_options)

add_executable(
  ActsExampleAdaptiveMultiVertexFinder
  AdaptiveMultiVertexFinderExample.cpp)
target_link_libraries(
  ActsExampleAdaptiveMultiVertexFinder
  PRIVATE ${_common_libraries} )

add_executable(
  New
  NewExample.cpp)
target_link_libraries(
  New
  PRIVATE ${_common_libraries} ActsExamplesIoRoot ActsExamplesPrinters)

add_executable(
  ActsExampleIterativeVertexFinder
  IterativeVertexFinderExample.cpp)
target_link_libraries(
  ActsExampleIterativeVertexFinder
  PRIVATE ${_common_libraries})
add_executable(
  ActsExampleVertexFitter
  VertexFitterExample.cpp)
target_link_libraries(
  ActsExampleVertexFitter
  PRIVATE ${_common_libraries})
add_executable(
  ActsTutorialVertexFinder
  TutorialVertexFinderExample.cpp)
target_link_libraries(
  ActsTutorialVertexFinder
  PRIVATE ${_common_libraries})

install(
  TARGETS
    New
    ActsExampleAdaptiveMultiVertexFinder
    ActsExampleIterativeVertexFinder
    ActsExampleVertexFitter
    ActsTutorialVertexFinder
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
