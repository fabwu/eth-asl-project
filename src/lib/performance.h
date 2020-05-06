#ifndef PERFORMANCE_H
#define PERFORMANCE_H

#define ENABLE_PERF_COUNTER

void __reset_flop_counter();

void __record_double_flops(const int amount);

long long int __nbr_double_flops(void);

#endif  // PERFORMANCE_H
