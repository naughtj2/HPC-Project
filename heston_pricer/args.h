/**
 * @file args.h
 * @brief Command-line options for model, product, and MC configuration.
 */
#ifndef PRICER_ARGS_H
#define PRICER_ARGS_H
#include <string>
#include <iostream>
#include "mc.h"
#include "options.h"
namespace pricer {
    /**
     * @brief Parsed CLI arguments and helpers.
    */
    struct Args {
        OptionSpec opt; 
        MCConfig cfg; 
        HestonSpec H; 
        Model model{Model::BS};

        /** @brief Map option name string to @ref OptionType. */
        static OptionType option_from_string(const std::string& s);
        
        /** @brief Parse command line and populate members. */
        void parse(int argc, char** argv);
    };
} 
#endif
