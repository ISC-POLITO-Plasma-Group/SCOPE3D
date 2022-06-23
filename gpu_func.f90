module gpu_func
      use cudafor  
      use cusparse
      type(cusparseHandle):: handle
      logical use_gpu
      integer(8):: pBufferSizeInBytes
      character(1), device, allocatable:: pBuffer(:)
      integer(8):: pBufferSizeInBytes1d
      character(1), device, allocatable:: pBuffer1d(:)
      
end module
