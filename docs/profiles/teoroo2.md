# The Teoroo2 Cluster

!!! note "Note for scheduling GPU jobs on Teoroo2" 
    
    The official nextflow does not support scheduling GPU resources on a local
    executor. On Teoroo, a custom nextflow build is avialble. Use `/sw/nf` 
    instead of `nextflow`, and set the GPUs available to a workflow with 
    `CUDA_VISIALBE_DEVICES`, namely:
    
    ``` bash
    export CUDA_VISIBLE_DEVICES=1,2,3,4
    /sw/nextflow main.nf -entry h2o_demo -profile teoroo2
    ```

 The Teoroo2 profiel use containerized runtimes as the [standard profile], the
 difference is that local copies of the images are available and used directly.
 Local singularity images in sif format are stored in the `/sw/pinnacle` folder.
 All processes are configured to run on single thread, and the pinn-based jobs
 are configured to run with one GPU.

## Profile

```groovy
--8<-- "profiles/teoroo2.config"
```
 
 [standard profile]: ../overview/#Standard_profile
