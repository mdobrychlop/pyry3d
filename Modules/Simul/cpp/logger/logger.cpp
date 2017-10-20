#include "logger.hpp"
#include "../extern.hpp"

#include <iostream>

void Logger::set_log(std::string filename){
    logfile.open(filename + "pyry.log", std::ios_base::app);
}

Logger::Logger(){

}

void Logger::write_to_logfile(int complex_index){
    logfile << "cx score for iteration ";
    logfile << complexes[complex_index].iteration_index;
    logfile << " ";
    logfile << complexes[complex_index].score;
    logfile << " components: restraints: ";
    logfile << (complexes[complex_index].restraints * restraints_penalty);
    logfile << " collisions: ";
    logfile << (complexes[complex_index].clashes    * clashes_penalty);
    logfile << " map filling: ";
    logfile << (complexes[complex_index].freespace  * freespace_penalty);
    logfile << " outbox atoms: ";
    logfile << (complexes[complex_index].outbox     * outbox_penalty);
    logfile << " density filling: ";
    logfile << (complexes[complex_index].density    * density_penalty);
    logfile << " symmetry penalty: ";
    logfile << (complexes[complex_index].symmetry   * symmetry_penalty);
    logfile << " chi2 penalty: ";
    logfile << (complexes[complex_index].chi2       * chi2_penalty);
    logfile << " rg penalty: ";
    logfile << (complexes[complex_index].rg           * rg_penalty);
    logfile << std::endl;
}   