#include "performance.h"

// TODO: Test if this works
static long long int nbr_double_flops = 0;

void __reset_flop_counter() {
#ifdef ENABLE_PERF_COUNTER
    nbr_double_flops = 0;
#endif
}

void __record_double_flops(const int amount) {
#ifdef ENABLE_PERF_COUNTER
    nbr_double_flops += amount;
#endif
}

long long int __nbr_double_flops(void) { return nbr_double_flops; }
