# Singularity

[Singularity](https://www.sylabs.io/docs/) is the preferred method for
running container on HPC systems.  Singularity can use existing Docker
images.

The first step is to download and convert the image

    singularity pull docker://wlandry/sdpb:2.1.0

This will create a file named `sdpb-2.1.0.simg` or `sdpb_2.1.0.sif`
depending on your version of Singularity.

Singularity should automatically mount your home directory.  If your
Singularity image is named `sdpb-2.1.0.simg`, you can invoke the SDPB
programs by prepending the command with

    singularity exec sdpb-2.1.0.simg

So to convert the XML input file at `/home/user/input.xml`, run the command

    singularity exec sdpb-2.1.0.simg mpirun -n 4 pvm2sdp 1024 /home/user/input.xml /home/user/input

This uses 4 cores when running pvm2sdp.  You can change that number to
match your own machine.

To find a primal-dual solution, 

    singularity exec sdpb-2.1.0.simg mpirun -n 4 sdpb --precision=1024 --procsPerNode=4 -s /home/user/input

In theory, Singularity can be used to run jobs across multiple nodes.
We are still working through some issues with that.
