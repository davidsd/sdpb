I'll have a more in-depth write-up, but I just wanted to explain at a high level quickly how the Karatsuba algorithm works, and how it's used here. 

Karatsuba multiplication is an algorithm that reduces the complexity of multiplication. Suppose I have two numbers,
*A* and *B*, whose product *C* I wish to compute. *A* can be expressed as a two-digit number with digits *A_0* and *A_1*, i.e. 
*A* = *A_0* + *A_1* *b*, where *b* is our base. *B* can be expressed in the same way. Then *C* = *C_0* + *C_1* *b* + *C_2* *b^2*. 
Gradeschool arithmetic tells us that *C_0* = *A_0* *B_0*, *C_1* = *A_0* *B_1* + *A_1* *B_0*, and *C_2* = *A_1* *B_1*. In other words,
there are four 1-digit multiplications we must compute in order to calculate one 2-digit multiplication. Anatoly Karatsuba noticed that
*C_1* = (*A_0* + *A_1)(*B_0* + *B_1*) - *C_0* - *C_2*: if you already know two of the products, you can determine the rest with a few
additions and only one more multiplication, for three multiplications in total. If *A* and *B* are more than 2-digit numbers,
this procedure can be done recursively, to get a complexity of O(N^log2(3)). 

In our case, these numbers are actually matrices, and there are billions of elements in play. Therefore, memory is a major limiter: we
cannot simply do the above Karatsuba trick in one line because we can't just store the matrices in temporary variables. To that end, 
my Karatsuba algorithm implementation is roughly as follows:

1. C_0 = A_0 B_0
2. A_0 += A_1
3. B_0 += B_1
4. C_1 = A_0 B_0
5. A_0 -= A_1
6. B_0 -= B_1
7. C_1 -= C_0
8. C_2 = A_1 B_1
9. C_1 -= C_2

However, there are more complications to do with memory. Because each multiplication result is independent and has to be treated as such,
if we were to just caculate all the possible multiplications first and then add and shuffle as above, we would need O(N^log2(3)) memory,
which is problematic when dealing with problems as large as SDPB often does. With that said, the final amount of memory that the result
is stored in is only O(N); much of the results from above merge and overlap with each other. The actual explanation of what that overlap
precisely is warrants visualizations, so for the purposes of this quick note it's enough to note that it exists. It's the most memory 
efficient for us to do the overlapping and merging as soon as possible. This gives us a memory complexit of O(N). Our implementation becomes:

1. C_0 = A_0 B_0
2. A_0 += A_1
3. B_0 += B_1
4. C_1 = A_0 B_0
5. A_0 -= A_1
6. B_0 -= B_1
7. C_1 -= C_0
9. Deal with overlap
8. C_2 = A_1 B_1
9. C_1 -= C_2
10. Deal with overlap

There are further complications due to the fact that we want to work at a certain working precision. In the context of the above Karatsuba,
that amounts to saying that we wish to calculate C_0 and C_1 but not C_2. In that case, it's advantageous to just directly calculate the grade
school way. Hence, there are two kinds of functions in the code: `gradeschool` and `karatsuba`. The choice of one or the other is just based
on the knowlege of whether C_2 is nonzero or not.

Beyond the `gradeschool` and `karatsuba` classification, there are two more categorizations. Functions which list only A and C are symmetric
cases of the above algorithms, if B = A^T; the algorithms themselves are roughly similar to the above, but are done carefully so as to
minimize computations given the symmetries. Functions are listed as either `_cpu`, `_gpu`, or without either classifier. The vanilla version
is a purely host-side version of the algorithm. The `_gpu` version is a purely device-side version of the algorithm. The `_cpu` version starts
on the host-side, but as soon as in its recursion it reaches a point where the problem can fit on the GPU's memory, it copies over to the
device and calls the `_gpu` version.

With all of that said, the only code that anyone *should* see except those mucking about in the Karatsuba itself are the helpful public 
wrapper functions that hide away the recursion.

A note on recursion: while an iterative way would have been possible, the recursion was the best choice because while it may
add a slight amount of overhead, it made the development much simpler. Debugging (mostly) was greatly simplified, and it allowed modularity
like the above `_cpu` functions. Also, the depth of recursion goes as log(N), so the overhead of recursion is negligible compared to the
massive multiplications that happen.
