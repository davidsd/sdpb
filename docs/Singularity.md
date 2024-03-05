# Singularity

Singularity
(since 2021 split into two forks, [SingularityCE](https://www.sylabs.io/docs/) and [Apptainer](https://apptainer.org/))
is the preferred method for running containers on HPC systems.
For laptop and desktop Linux machines, it is often a more user-friendly way of running containers than Docker.
Singularity can use existing Docker images.

## Install Singularity
On Debian and Ubuntu, Singularity is available in the default
repository in the `singularity-container` package.  So installing
Singularity should be as simple as running

    sudo apt-get install singularity-container

On Redhat and Centos, Singularity is available through EPEL.
The Singularity website has
detailed [installation instructions](https://www.sylabs.io/guides/3.2/user-guide/installation.html#install-the-centos-rhel-package-using-yum).

On HPC systems, if Singularity is available, it is probably installed
as a module.  Run

    module av

and look for `singularity`.

## Convert from Docker image
Once Singularity is installed, the next step is to download and
convert the image

    singularity pull docker://bootstrapcollaboration/sdpb:master

Depending on your version of Singularity, this will create a file
named `sdpb-master.simg` or `sdpb_master.sif`.

If you want to build Singularity image from current sources, you should build Docker image first
and push it to [local registry](https://docs.docker.com/registry/):

    # Build Docker image
    docker build . --tag sdpb
    
    # Run local registry
    docker run -d -p 5000:5000 --name registry registry:2.8
    
    # Push to local regisrty
    docker image tag sdpb localhost:5000/sdpb:master
    docker push localhost:5000/sdpb:master
    
    # Pull from local registry and convert to Singularity image
    singularity pull --no-https docker://localhost:5000/sdpb:master
    
    # Stop local registry and remove all data
    docker container stop registry && docker container rm -v registry

## Run Singularity image
Singularity should automatically mount your home directory.  If your
Singularity image is named `sdpb-2.7.0.sif`, you can invoke the SDPB
programs by prepending the command with

    singularity exec sdpb-2.7.0.sif

So to convert the JSON input file at `/home/user/input.json`, run the command

    singularity exec sdpb-2.7.0.sif mpirun -n 4 pmp2sdp --precision 1024 -i /home/user/input.json -o /home/user/input

This uses 4 cores when running pmp2sdp. You can change that number to
match your own machine.

To find a primal-dual solution, 

    singularity exec sdpb-2.7.0.sif mpirun -n 4 sdpb --precision=1024 -s /home/user/input

In theory, Singularity can be used to run jobs across multiple nodes.
We have not been able to make that work yet.  So for large, multi-node
jobs, it is still recommended to build from source using the
instructions in [Install.md](../Install.md).
