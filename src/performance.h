#ifndef PERFORMANCE_H
#define PERFORMANCE_H

#define ENABLE_PERF_COUNTER true

inline long long int nbr_double_flops = 0;

inline void __reset_flop_counter() {
#if ENABLE_PERF_COUNTER
  nbr_double_flops = 0;
#endif
}

inline void __record_double_flops(const int amount) {
#if ENABLE_PERF_COUNTER
  nbr_double_flops += amount;
#endif
}

#endif  // PERFORMANCE_H
