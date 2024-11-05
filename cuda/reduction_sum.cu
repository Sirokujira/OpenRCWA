/*
reduction_sum.cu

tid : thread number (=0,...)
n   : array size (<=1024)
s   : array (shared memory)
sum : sum (output)
*/
__device__
void reduction_sum(int tid, int n, real_t *s, real_t *sum)
{
	__syncthreads();

	for (int stride = (n + 1) >> 1; n > 1; stride = (stride + 1) >> 1, n = (n + 1) >> 1) {
		if (tid + stride < n) {
			s[tid] += s[tid + stride];
		}
		__syncthreads();
	}

	if (tid == 0) {
		*sum = s[0];
	}
}
