#include <stdio.h>
#include <stdbool.h>

#include "test.h"
#include "itree/iv.h"
#include "itree/search.h"

bool test_all(){
    return test_iv() &&
           test_search();
}
