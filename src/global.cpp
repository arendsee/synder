#include "global.h"

void check_in_offset(size_t start, size_t stop){
    if(start == 0 && Offsets::in_start > 0){
        fprintf(stderr, "0 starts are illegal in in 1-based input\n");
        fprintf(stderr, "Are you sure your input is 1-based?\n");
        exit(EXIT_FAILURE);
    }
    if(stop == 0 && Offsets::in_stop > 0){
        fprintf(stderr, "0 stops are illegal in in 1-based input\n");
        fprintf(stderr, "Are you sure your input is 1-based?\n");
        exit(EXIT_FAILURE);
    }
}
