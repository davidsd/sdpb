# Singularity

[Singularity](https://www.sylabs.io/docs/) is the preferred method for
running containers on HPC systems.  For laptop and desktop Linux
machines, it is often a more user friendly way of running containers
than Docker.  Singularity can use existing Docker images.

On Debian and Ubuntu, Singularity is available in the default
repository in the `singularity-container` package.  So installing
Singularity should be as simple as running

    sudo apt-get install singularity-container
    
On Redhat and Centos, Singularity is available through EPEL.  The Singularity web site has detailed [installation instructions](https://www.sylabs.io/guides/3.2/user-guide/installation.html#install-the-centos-rhel-package-using-yum).

On HPC systems, if Singularity is available, it is probably installed
as a module.  Run

    module av

and look for `singularity`.

Once Singularity is installed, the next step is to download and
convert the image

    singularity pull docker://wlandry/sdpb:2.1.2

Depending on your version of Singularity, this will create a file
named `sdpb-2.1.2.simg` or `sdpb_2.1.2.sif`.

Singularity should automatically mount your home directory.  If your
Singularity image is named `sdpb-2.1.2.simg`, you can invoke the SDPB
programs by prepending the command with

    singularity exec sdpb-2.1.2.simg

So to convert the XML input file at `/home/user/input.xml`, run the command

    singularity exec sdpb-2.1.2.simg mpirun -n 4 pvm2sdp 1024 /home/user/input.xml /home/user/input

This uses 4 cores when running pvm2sdp.  You can change that number to
match your own machine.

To find a primal-dual solution, 

    singularity exec sdpb-2.1.2.simg mpirun -n 4 sdpb --precision=1024 --procsPerNode=4 -s /home/user/input

In theory, Singularity can be used to run jobs across multiple nodes.
We have not been able to make that work yet.  So for large, multi-node
jobs, it is still recommended to build from source using the
instructions in [Install.md](../Install.md).
