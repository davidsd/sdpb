# Algorithm

## Define shared memory window sizes and split factors.

Shared memory windows (`double` elements):
1. `L_X_inv`: `size = sum(dim^2) * num_primes(trmm)`
2. `L_Y`: same size as `L_X_inv`
3. `L_X_inv F_p`: `size = sum(dim^2) * num_primes(trmm) * #(p)`
4. `(L_X_inv F_p) L_Y`: reuse the window for `L_X_inv F` 
5. `P` (syrk input) 
6. `S` (syrk output)

where:
- `dim` is a block dimension, `sum(dim^2)` is a sum over all blocks on a node.
- `#(p)` is how many `F_p`'s we process simultaneously (`p=1..primal_dimension`).
  For example: if `#(p)=3`, we process `F_{1..3}` and then `F_{4..6}`, etc.
- `num_primes(trmm)` depends on precision and the maximum number of additions `max(dim^2) * #(p)`
- `num_primes(syrk)` depends on precision and the maximum number of additions in syrk, which is `height(P) = sum(dim^2)`.


Other temporary allocations (`BigFloat` elements):

1. `W = (F_p F_{p+1}...)` (stacked horizontally): `size = sum(dim^2) * #(p)`. It's a `vector<Block_Diagonal_Matrix>`, where vector index is `p`.
2. `L_X_inv F_p`. Now we think of them as stacked vertically, but in fact we reuse the same `vector<Block_Diagonal_Matrix>`. We fill it directly from residues.
3. `(L_X_inv F_p) L_Y`: in-place, reuse the same structure.
4. `P`: fill from `(L_X_inv F_p) L_Y` by putting all elements from each `p` into a single column.
As long as we don't split `S` window (syrk output), we need all elements of `P`.
This means `size = sum(dim^2) * primal_dimension`.
**TODO:** we don't need the same block structure, let's assign the same number of rows to each CPU.
**TODO:** in extreme case, we might fill only parts of `P` and pay extra cost for repeating the same trmm calls. Let's not think about it now.
5. `S`: `size = primal_dimension^2`. This is not temporary, by the way.

Suppose that we have some free memory.

Constant size: `P`, `S`
Variable size: `W`, residues for `L_X_inv`, `L_Y`, `W`, `P`, `S`.
(`L_X_inv` and `L_Y` sizes vary due to `num_primes`)

Straightforward way:
```
num_primes(syrk) = num_primes(sum(dim^2))
for split_factor(S): 1..primal_dimension:
  for split_factor(trmm): 1..primal_dimension:
    delta_p = ceil(primal_dimension / split_factor(trmm))
    #(W) = sim(dim^2) * delta_p 
    num_primes(trmm) = num_primes(max(dim^2) * delta_p)
    #(residues L_X_inv) = sim(dim^2) * num_primes(trmm) 
    #(residues L_Y) = sim(dim^2) * num_primes(trmm) 
    #(residues W) = #(W) * num_primes(trmm)
    for slpit_factor(P): 1..primal_dimension
      #(residues P) = num_primes(syrk) * sum(dim^2) / split_factor(P)
      if #(W)*bigfloat + (#r(L_X_inv) + #r(L_Y) + max(#r(W), #r(P)))*double < mem_limit:
        return split factors
```

This is up to `O(primal_dimension)^3`, let's optimize:
```
num_primes(syrk) = num_primes(sum(dim^2))
min_primes(trmm) = num_primes(max(dim^2))
max_primes(trmm) = num_primes(max(dim^2) * primal_dimension)

total_trmm_bytes = sum(dim^2) * (delta_p * bigfloat + (delta_p + 2) * num_primes(trmm) * double) <= mem_limit

delta_p_min = floor((mem_limit/sum(dim^2) - 2 * max_primes(trmm)) / (bigfloat + max_primes(trmm))
delta_p_max = floor((mem_limit/sum(dim^2) - 2 * min_primes(trmm)) / (bigfloat + min_primes(trmm))

for split_factor(S): 1..primal_dimension:
  S_window_width = ceil(primal_dimension/split_factor(S)) 
  #r(S) = S_window_width^2 * num_primes(syrk)
  min #r(P) = 2 * S_window_width * num_primes(syrk)
  if #r(S) + min #r(P) > min(mem_limit, shmem_limit):
    continue
  
  for split_factor(trmm): primal_dimension/delta_p_max .. primal_dimension/delta_p_min:
    delta_p = ceil(primal_dimension / split_factor(trmm))
    if total_trmm_bytes(delta_p) <= mem_limit && shmem_trmm_bytes(delta_p) < shmem_limit:
      define split_factor(P) to fit shmem_limit
      return split factors
```


## Steps

Before solver iterations: determine split factors etc. via
`get_compute_S_config()`.
Allocate shared memory buffer(s) with size(s) specified by this `Compute_S_Config`.

Within each step - `initialize_P()`:
1. Compute L_X_inv via El::Trsm().
2. Normalize and shift (=multiply by 2^N) rows of L_X_inv and columns of L_Y.
3. Compute residues for L_X_inv and L_Y, put them into two `Vector_Matrix_Residues_Window`'s.
4. Allocate `Block_Matrix P`.
5. Loop over split_factor(trmm):
   1. Copy F_p for several p's to `vector<Block_Diagonal_Matrix> G`
   2. Normalize and shift **columns** of G (each p, each block separately).
   3. Compute G residues, put into window. The matrices with the same block index but different p index are stacked **horizontally** as a single matrix. So we can store everything in `Vector_Matrix_Residues_Window` again. Maybe we need views for each submatrix (to make restoring from residues easier). 
   4. Compute G := L_X_inv G (via BLAS), restore from residues into `vector<Block_Diagonal_Matrix> G`. 
   5. Remove normalization for `G` (using `L_X_inv` row norms and `G` column norms).
   6. Normalize and shift **rows** of `G`.
   7. Compute G residues. Put into `Vector_Matrix_Residues_Window` (another?). The matrices with the same block index but different p index are stacked **vertically** as a single matrix.
   8. Compute G := G L_Y (via BLAS), restore from residues into `G`. 
   9. Remove normalization (using `G` row norms and `L_Y` column norms)
   10. Reshape each G_p to the p-th column of P.
   **TODO:** can we restore from residues directly into P and then remove normalization? It's less convenient and doesn't save much memory. But saves some communication.

Then syrk P (same as before).
   
## First touch

**TODO:** how should we do a first touch, when shared window is reused? Shall we just use NUMA nodes instead of physical nodes to ensure memory locality?
It's probably easier and better to use NUMA nodes instead of physical nodes.

Disadvantages of NUMA nodes:
- If we have large blocks that need more cores than a single NUMA node has, we cannot do a good load balancing.
- reduce-scatter S time is proportional to the number of nodes. Sending data is usually fast, but serializing and deserializing takes time.

What if we work with physical nodes, but want to optimize memory access?
Assume that we have 2 NUMA nodes. Let's allocate the first half of the window on the first node, and the second half on the second node.
Now let's take all BLAS jobs (for trmm and syrk) and assign each one to a NUMA node (e.g. all jobs with the first prime belong to the first node, and all jobs with the last prime belong to the second node). After that, distribute all jobs from a single NUMA node among all cores from this node.

If we reuse the same buffer for trmm and for syrk, it could be tricky to align them, in particular when we have two input windows for syrk. In that case, one has to put the first half of primes for both windows into the first hasl of the underlying buffer. To keep constant `prime_stride`, we have to interleave two windows: e.g.:

```[(win1, prime1)(win2, prime1)(win1, prime2)(win2, prime2)...]```

It's easier to start with separate buffers for each window.
NB: adjust `compute_S_config()` and memory estimates in it.
