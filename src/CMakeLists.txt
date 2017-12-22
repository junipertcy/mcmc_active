cmake_minimum_required(VERSION 2.8)

add_executable(mcmc_active ActiveLearnMain.cpp Utility.cpp Graph.cpp Learner.cpp MCMC.cpp MCMCAlg.cpp TypeModel.cpp)

set_target_properties(mcmc_active PROPERTIES COMPILE_DEFINITIONS "OUTPUT_HISTORY=0")

if (Boost_FOUND)
    target_link_libraries(mcmc_active ${Boost_LIBRARIES})
endif (Boost_FOUND)
