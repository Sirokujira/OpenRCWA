/*
rank2device.cu
*/

// rank -> device number
int rank2device(int comm_rank, int nhost, const int ndevice[])
{
	int device = 0;

	if (nhost <= 1) {
		// single node
		device = comm_rank;
	}
	else {
		// cluster
		device = -1;
		int rank = -1;
		for (int ihost = 0; ihost < nhost; ihost++) {
			for (int idevice = 0; idevice < ndevice[ihost]; idevice++) {
				if (++rank == comm_rank) {
					device = idevice;
					break;
				}
			}
			if (device >= 0) {
				break;
			}
		}
		if (device < 0) device = 0;
	}

	int num_device;
	cudaGetDeviceCount(&num_device);
	if (device >= num_device) device = num_device - 1;

	return device;
}
