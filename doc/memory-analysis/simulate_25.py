from cachesimulator import *
import itertools


def simulate_quadtree_depth(curr_quadtree_depth, image_size, rbs, var_image, var_prepared_domain_blocks, var_db_sums,
                            var_db_sums_squared, var_prepared_range_block, ms: MemorySystem):
    nbr_dbs = 2 ** (2 * curr_quadtree_depth - 2)
    nbr_rbs = 2 ** (2 * curr_quadtree_depth)
    dbs = int(image_size / (2 ** (curr_quadtree_depth - 1)))

    def BLOCK_CORD_REL_X(block_id, block_size):
        return (block_id % (int(image_size / block_size))) * block_size

    def BLOCK_CORD_REL_Y(block_id, block_size):
        return int((block_id / (image_size / block_size))) * block_size

    def BLOCK_CORD_X(block_id, block_size):
        return (block_id % (image_size / block_size))

    def BLOCK_CORD_Y(block_id, block_size):
        return int((block_id / (image_size / block_size)))

    def sim_prep_db_scale_block(db_idx, block_rel_y, block_rel_x):
        scaled_image_idx = 0
        original_image_idx = block_rel_y * image_size + block_rel_x

        for y in range(rbs):
            for x in range(rbs):
                ms.read_var(var_image[original_image_idx])
                ms.read_var(var_image[original_image_idx + 1])
                ms.read_var(var_image[original_image_idx + image_size])
                ms.read_var(var_image[original_image_idx + image_size + 1])
                ms.write_var(var_prepared_domain_blocks[db_idx * rbs * rbs])

                scaled_image_idx += 1
                original_image_idx += 2
            original_image_idx += image_size

    def sim_prepare_domain_blocks_norotation():
        for i in range(nbr_dbs):
            sim_prep_db_scale_block(i, BLOCK_CORD_REL_X(i, dbs), BLOCK_CORD_REL_Y(i, dbs))

            for j in range(rbs * rbs):
                ms.read_var(var_prepared_domain_blocks[i * rbs * rbs + j])

            ms.write_var(var_db_sums[i])
            ms.write_var(var_db_sums_squared[i])

    def sim_load_block(rb_idx):
        rel_x = BLOCK_CORD_REL_X(rb_idx, rbs)
        rel_y = BLOCK_CORD_REL_Y(rb_idx, rbs)

        idx = 0
        idx_in_image = rel_y * image_size + rel_x
        for i in range(rbs):
            for j in range(rbs):
                ms.read_var(var_image[idx_in_image])
                ms.write_var(var_prepared_range_block[idx])
                idx += 1
                idx_in_image += 1
            idx_in_image += image_size - rbs

    sim_prepare_domain_blocks_norotation()
    for idx_rb in range(nbr_rbs):
        sim_load_block(idx_rb)

        a = BLOCK_CORD_REL_Y(idx_rb, rbs)
        b = BLOCK_CORD_REL_X(idx_rb, rbs)
        rtd_start_rb = a * image_size + b

        for idx_db in range(nbr_dbs):
            idx_prep_db = idx_db * rbs * rbs

            rtd_idx_rb = rtd_start_rb
            dbs_i = 0
            for i in range(rbs):
                dbs_j = 0
                for j in range(0, rbs, 2):
                    idx_0_db1 = dbs_i + j
                    idx_0_db2 = idx_0_db1 + 1
                    idx_90_db1 = rbs * rbs - dbs_j - rbs + i
                    idx_90_db2 = idx_90_db1 - rbs
                    idx_180_db1 = rbs * rbs - dbs_i - j - 1
                    idx_180_db2 = idx_180_db1 - 1
                    idx_270_db1 = dbs_j + rbs - i - 1
                    idx_270_db2 = idx_270_db1 + rbs

                    idx_rb2 = rtd_idx_rb + 1

                    ms.read_var(var_image[rtd_idx_rb])
                    ms.read_var(var_image[idx_rb2])

                    ms.read_var(var_prepared_domain_blocks[idx_prep_db + idx_0_db1])
                    ms.read_var(var_prepared_domain_blocks[idx_prep_db + idx_0_db2])
                    ms.read_var(var_prepared_domain_blocks[idx_prep_db + idx_90_db1])
                    ms.read_var(var_prepared_domain_blocks[idx_prep_db + idx_90_db2])
                    ms.read_var(var_prepared_domain_blocks[idx_prep_db + idx_180_db1])
                    ms.read_var(var_prepared_domain_blocks[idx_prep_db + idx_180_db2])
                    ms.read_var(var_prepared_domain_blocks[idx_prep_db + idx_270_db1])
                    ms.read_var(var_prepared_domain_blocks[idx_prep_db + idx_270_db2])

                    rtd_idx_rb += 2
                    dbs_j = dbs_j + rbs + rbs
            rtd_idx_rb += image_size - rbs
            dbs_i += rbs


def simulate_full_depth(image_size, MAX_QUADTREE_DEPTH=8):
    """
        OUR HARDWARE
    """
    ms = MemorySystem()
    ms.add_cache(Cache(total_size=32768, block_size=64, associativity=8))
    ms.add_cache(Cache(total_size=262144, block_size=64, associativity=4))
    ms.add_cache(Cache(total_size=8388608, block_size=64, associativity=16))

    """
        CONSTANTS
    """
    UPPER_BOUND_DOMAIN_BLOCKS = 4 ** MAX_QUADTREE_DEPTH
    UPPER_BOUND_RANGE_BLOCKS = 4 ** (MAX_QUADTREE_DEPTH + 1)

    var_full_image = ms.create_doubles(image_size * image_size)
    var_prepared_domain_blocks = ms.create_doubles(int(image_size * image_size / 4))
    var_db_sums = ms.create_doubles(UPPER_BOUND_DOMAIN_BLOCKS)
    var_db_sums_squared = ms.create_doubles(UPPER_BOUND_DOMAIN_BLOCKS)
    var_prepared_range_block = ms.create_doubles(int(image_size * image_size / 4))

    results = {}
    for curr_quadtree_depth in range(1, MAX_QUADTREE_DEPTH):
        rbs = int(image_size / (2 ** curr_quadtree_depth))
        if 2 ** curr_quadtree_depth < size and rbs >= 4:
            simulate_quadtree_depth(curr_quadtree_depth, image_size, rbs, var_full_image, var_prepared_domain_blocks,
                                    var_db_sums, var_db_sums_squared, var_prepared_range_block, ms)

            # Capture state
            results[curr_quadtree_depth] = ms.caches

    return results


for size in [64, 128, 256, 512]:
    print("IMAGE SIZE: %s x %s" % (size, size))
    results = simulate_full_depth(size, MAX_QUADTREE_DEPTH=8)

    print("depth;l1misses;l1accesses;l1mr;l2misses;l2accesses;l2mr;l3misses;l3accesses;l3mr")
    for depth in results:
        l1miss = results[depth][0].stats.misses
        l1acc = results[depth][0].stats.accesses
        l2miss = results[depth][1].stats.misses
        l2acc = results[depth][1].stats.accesses
        l3miss = results[depth][2].stats.misses
        l3acc = results[depth][2].stats.accesses
        print("%d;%d;%d;%f;%d;%d;%f;%d;%d;%f" % (
            depth, l1miss, l1acc, l1miss / l1acc, l2miss, l2acc, l2miss / l2acc, l3miss, l3acc, l3miss/l3acc))
