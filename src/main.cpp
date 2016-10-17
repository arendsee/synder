#include "commands.h"
#include "version.h"

#include <string>
#include <iostream>

int Offsets::syn_start;
int Offsets::syn_stop;
int Offsets::in_start;
int Offsets::in_stop;
int Offsets::out_start;
int Offsets::out_stop;

int main(int argc, char *argv[])
{
    const option::Descriptor usage[] = {
        {
            UNKNOWN, 0, "", "", option::Arg::None,
            "synder - Explore genomes using a synteny map\n"
            "Usage:\n"
            "  synder [SUBCOMMAND] [OPTIONS]\n"
            "  synder [SUBCOMMAND] --help\n"
            "Subcommands:\n"
            "  search - predict search intervals\n"
            "  filter - remove links that disagree with the synteny map\n"
            "  map    - trace intervals across genomes\n"
            "  count  - count overlaps\n"
            "  dump   - print all blocks with contiguous set ids\n"
            "Arguments:"
        },
        {
            HELP, 0, "h", "help", option::Arg::None,
            "  -h, --help \tPrint usage and exit."
        },
        {
            VERSION, 0, "v", "version", option::Arg::None,
            "  -v, --version \tPrint version and exit."
        },
        {0,0,0,0,0,0}
    };

    SHIFT 

    BOILERPLATE

    if(options[VERSION]){
        std::cout << SYNDER_VERSION << std::endl;
        return EXIT_SUCCESS;
    }

    if(argv[0][0] == '-'){
        std::cerr << "Please provide a subcommand\n";  
        return EXIT_FAILURE;
    }

    SHIFT 
    std::string subcommand(parse.nonOption(0));

    if( subcommand == "search" ) return !subcommand_search ( argc, argv );
    if( subcommand == "filter" ) return !subcommand_filter ( argc, argv );
    if( subcommand == "map"    ) return !subcommand_map    ( argc, argv );
    if( subcommand == "count"  ) return !subcommand_count  ( argc, argv );
    if( subcommand == "dump"   ) return !subcommand_dump   ( argc, argv );

    std::cerr << "Failed to parse arguments, subcommand '" << parse.nonOption(0) << "' is not defined\n\n";
    option::printUsage(std::cout, usage);

    return EXIT_FAILURE;
}
