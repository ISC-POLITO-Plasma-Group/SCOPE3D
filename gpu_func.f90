module gpu_func
      use cudafor  
      use cusparse

      type(cusparseHandle) :: handle
      logical use_gpu
      integer(8) :: pBufferSizeInBytes
      character(1), device, allocatable :: pBuffer(:)
      integer(8) :: pBufferSizeInBytes1d
      character(1), device, allocatable :: pBuffer1d(:)
      integer(kind=cuda_stream_kind) :: streamd1, streamd2, streamdd1, streamdd2, streamddd1
      integer :: stridd1, stridd2, striddd1, striddd2, stridddd1 
end module
