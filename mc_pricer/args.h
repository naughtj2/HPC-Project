/**
 * @file args.h
 * @brief Command-line argument parsing for option specs and Monte Carlo config.
 */
#ifndef PRICER_ARGS_H
#define PRICER_ARGS_H

#include <string>
#include <iostream>
#include "mc.h"
#include "options.h"

namespace pricer {

    /**
 * @brief Simple container for option specification and MC configuration with CLI parsing.
 */
    struct Args {
        OptionSpec opt;
        MCConfig cfg;
        /** @brief Map a string token to an @ref OptionType. */
        static OptionType option_from_string(const std::string& s);
        /** @brief Parse command-line arguments in place.
         *  @param argc Argument count
         *  @param argv Argument vector
         */
        void parse(int argc, char** argv);
    };

} 

#endif
