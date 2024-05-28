#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "lix.h"


void usage() {
    fprintf(stderr, "Usage: lix [OPTIONS]\n");
    fprintf(stderr, "    Build lix index from tbi index\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -i, --input FILE       Input tbi index file\n");
    fprintf(stderr, "    -o, --output FILE      Output lix index file\n");
    fprintf(stderr, "    -m, --min-shift INT    Set minimum interval size for lix indices to 2^INT [14]\n");
    fprintf(stderr, "    -h, --help             Print this help message and exit\n");
}

int main(int argc, char *argv[]) {
    if (argc == 1) {
        usage();
        return 0;
    }

    // Define your command line options here
    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"min-shift", required_argument, 0, 'm'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    // Variables to store command line option values
    char *input_file = NULL;
    char *output_file = NULL;
    int min_shift = 14;

    // Parse command line options
    int option;
    while ((option =
        getopt_long(argc, argv, "i:o:m:h", long_options, NULL)) != -1)
    {
        switch (option) {
            case 'i':
                input_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'm':
                // Your code to parse min-shift option goes here
                min_shift = atoi(optarg);
                break;
            case 'h':
                usage();
                return 0;
            default:
                fprintf(stderr, "Error: Invalid option\n");
                usage();
                exit(EXIT_FAILURE);
        }
    }

    // Check if required options are provided
    if (input_file == NULL || output_file == NULL) {
        fprintf(stderr, "Error: Input and output files must be specified\n");
        exit(EXIT_FAILURE);
    }

    int ret = lix_build(input_file, min_shift, output_file);
    if (ret != 0) {
        fprintf(stderr, "Error: Failed to build lix index, return code is %d\n", ret);
        exit(EXIT_FAILURE);
    }

    return 0;
}
