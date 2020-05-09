#ifndef PERFORMANCE_H
#define PERFORMANCE_H

// XXX undef this to disable perf counter
#define ENABLE_PERF_COUNTER

#ifdef ENABLE_PERF_COUNTER
extern long long int nbr_double_flops;
#define __reset_flop_counter() (nbr_double_flops = 0)
#define __record_double_flops(amount) (nbr_double_flops += amount)
#define __nbr_double_flops() (nbr_double_flops)
#else
#define __reset_flop_counter() ((void)0)
#define __record_double_flops(amount) ((void)0)
#define __nbr_double_flops() (0)
#endif

#endif  // PERFORMANCE_H
