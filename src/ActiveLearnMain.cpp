/* ~~~~~~~~~~~~~~~~ Notes ~~~~~~~~~~~~~~~~ */
//
// The code is written by XXX.
//
/* ~~~~~~~~~~~~~~~~ Includes ~~~~~~~~~~~~~~~~ */
// Boost
#include <boost/program_options.hpp>

// Program headers
#include "TypeModel.h"
#include "MCMCAlg.h"

namespace po = boost::program_options;

int main(int argc, char const *argv[]) {
    std::string graph_path;
    unsigned int runs;
    unsigned int numTypeInModel;
    unsigned int numPhase;
    unsigned int numTop;

    po::options_description description("Options");
    description.add_options()
            /* Options go here */
            ("input_graph,i", po::value<std::string>(&graph_path), "Path to the graph file (*.gml).")
            ("num_of_runs,r", po::value<unsigned int>(&runs)->default_value(3), "Number of independent runs.")
            ("num_of_groups,g", po::value<unsigned int>(&numTypeInModel)->default_value(3), "Number of groups in the model, which can be different from the number of types in the synthetic data.")
            ("num_of_phases,p", po::value<unsigned int>(&numPhase)->default_value(32), "Number of phases for active learning, each phase returns a fixed number of vertices for query. The number of vertices to query is indicated by the \"numTop\" variable.")
            ("num_nodes_to_query,q", po::value<unsigned int>(&numTop)->default_value(1), "Number of vertices to query in each phase.");
    po::variables_map var_map;
    po::store(po::parse_command_line(argc, argv, description), var_map);
    po::notify(var_map);

    if (var_map.count("help") > 0 || argc == 1) {
        std::cout << "MCMC algorithms for active learning.\n";
        std::cout << "Usage:\n"
                  << "  " + std::string(argv[0]) + " [--option_1=value] [--option_s2=value] ...\n";
        std::cout << description;
        return 0;
    }
    if (var_map.count("input_graph") == 0) {
        std::cout << "input_graph (*.gml file) is required (-i flag)\n";
        return 1;
    }

    // "modelType" indicates the type of the underlying block model. We have quite a lot block models available in this package.
    unsigned int modelType = 1;
    bool groupCorrected = false;

    // "frozentypes" indicates the types/groups that are fixed. For example, in the "sea food web" network, the "resource" type is fixed (not for learn).
    set<unsigned int> frozentypes;

    //* Active learning algorithm *//
    unsigned int numOptInit = 50; // number of initials for optimization
    unsigned int numLearnerInit = 50; // number of initials for active learning
    unsigned int numOptStep = 1000;// number of iterations in each initial for optimization
    unsigned int numLearnerStep = 1000; // number of iterations in each initial for active learning
    unsigned int learningMethod = 1; // indicates which active learner is used: 1-MutualInfo, 2-AvgAgree, 3-RandomLearner, 4-TopSeqLearner, 5-MaxDegree, 6-MaxBtwn.

    for (unsigned int i = 0; i < runs; ++i) {
        clog << "Run # " << i + 1 << " is running..." << endl;
        std::unique_ptr<Graph>graph = std::make_unique<Graph>(graph_path);
        std::unique_ptr<MCMCAlg>mcmc_algorithm = std::make_unique<MCMCAlg>(
                *graph, numTypeInModel, frozentypes, numOptInit, numOptStep, numLearnerInit,
                numLearnerStep, numPhase, numTop, learningMethod, modelType,
                groupCorrected);
        mcmc_algorithm->runMCMCAlg();

    }
    return 0;
};


