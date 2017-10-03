#include "types.h"

std::vector<long> to_one_base(std::vector<long> x){
    for(long i = 0; i < x.size(); i++){
        x[i]++;
    }
    return x;
}
