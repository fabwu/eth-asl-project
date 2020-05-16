#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "../lib/tsc_x86.h"
#include "../lib/types.h"

double rtd_generic_with_rot90(const struct image_t *image,
                              const struct image_t *domain_block,
                              const struct block_t *range_block) {
    const int dbs = domain_block->size;
    double rtd_sum1 = 0;
    double rtd_sum2 = 0;

    size_t idx_rb1 = range_block->rel_y * image->size + range_block->rel_x;
   
    int init = dbs*dbs - dbs; 

    for (size_t i = 0; i < dbs; i++) {
        size_t idx_db1 = init + i;
        for (size_t j = 0; j < dbs; j += 2) {
            int idx_db2 = idx_db1 - dbs; 
            int idx_rb2 = idx_rb1 + 1;
            double ri1 = image->data[idx_rb1];
            double ri2 = image->data[idx_rb2];
            double di1 = domain_block->data[idx_db1];
            double di2 = domain_block->data[idx_db2];
            
            rtd_sum1 = fma(ri1, di1, rtd_sum1);
            rtd_sum2 = fma(ri2, di2, rtd_sum2);
            
            idx_rb1 = idx_rb1 + 2;
            idx_db1 = idx_db2 - dbs;
        }
        idx_rb1 += image->size - dbs;
    }


    return rtd_sum1 + rtd_sum2;
}
#define CYCLES_REQUIRED 1e8

int main(int argc, char const *argv[]) {
    myInt64 i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = 2;
    double ret;

    struct image_t img = make_image(256, 1);
    struct image_t db = make_image(16, 1);
    struct block_t rb =  make_block(3*16, 3*16,16,16);
    
    while(1) {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i) {
            ret = rtd_generic_with_rot90(&img, &db, &rb);
        }
        cycles = stop_tsc(start);

        //printf("%lld %lld\n", cycles, num_runs);
        if(cycles >= CYCLES_REQUIRED) break;

        if(num_runs == 0) {
            printf("OVERFLOW!!!");
            break;
        }

        num_runs *= 2;
    }
    start = start_tsc();

    for (i = 0; i < num_runs; ++i) {
        ret = rtd_generic_with_rot90(&img, &db, &rb);
    }

    cycles = stop_tsc(start)/num_runs;
    
    ret = rtd_generic_with_rot90(&img, &db, &rb);
    if(ret != 3924013) {
        printf("ERROR!");
    } else {
        printf("ret: %f cycles: %lld\n", ret, cycles);
    }
}
