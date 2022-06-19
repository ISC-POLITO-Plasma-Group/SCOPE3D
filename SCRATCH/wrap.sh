set -ex
if [[ $OMPI_COMM_WORLD_RANK == 0 ]]; then
	nsys profile --trace=nvtx,cuda,openacc --capture-range=nvtx --nvtx-capture=timestep --env-var=NSYS_NVTX_PROFILER_REGISTER_ONLY=0  $@
       	
else
	$@
fi

