#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <cuda_runtime.h>

class GPUManager
{
    private:
        int numgpus;
        std::vector<__global__ void*()> 
    public:
        GPUManager(int numGPUs);

        void passToDevice(int devID, void* pointer, size_t size);
        void

};