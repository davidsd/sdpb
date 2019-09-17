# Docker

Docker is a relatively easy way to run binaries without having to
build SDPB.  This can be a convenient way to run SDPB on your laptop
or desktop.  For HPC centers, most of them support Singularity rather
than Docker, so you should read the guide for
[Singularity](Singularity.md).  There are instructions at the
[Docker website](https://www.docker.com/get-started) on how to install
it on Windows, Mac, or Linux.

There is a docker image on [Docker Hub](https://hub.docker.com/) named
`wlandry/sdpb:2.2.0` containing binaries for `scalar_blocks`,
`pvm2sdp`, `sdp2input`, and `sdpb`.

Suppose you have an input file `/my/project/input.xml`. To use this
file, and be able to write output, we will make the `/my/project`
directory visible to the docker image in the location `/usr/local/share/sdpb`.

The structure of the command to run docker is

    docker run <options> <image> <command>
    
In this case, `<options>` will be used to mount your directory in a
place that docker will see it.

    -v /my/project/:/usr/local/share/sdpb/
    
`<image>` is the image name.

        wlandry/sdpb:2.2.0

`<command>` is the command that you would normally use to run the SDPB
commands (see [Usage.md](Usage.md)).  The directory containing the
input file is mounted as `/usr/local/share/sdpb`.  So we first run `pvm2sdp` to
convert from xml

    mpirun --allow-run-as-root -n 4 pvm2sdp 1024 /usr/local/share/sdpb/input.xml /usr/local/share/sdpb/input
    
`mpirun` runs as root inside the docker container.  Running `mpirun` as
root is normally dangerous, but it is safe to do so inside the
container.  To allow mpirun to run as root, we add the option
`--allow-run-as-root`.  This uses 4 cores when running pvm2sdp.  You
can change that number to match your own machine.  Putting it all
together on a single line

    docker run -v /my/project/:/usr/local/share/sdpb wlandry/sdpb:2.2.0 mpirun --allow-run-as-root -n 4 pvm2sdp 1024 /usr/local/share/sdpb/input.xml /usr/local/share/sdpb/input

Running this command will populate the directory `/my/project/input`.
To search for primal-dual solutions

    docker run -v /my/project/:/usr/local/share/sdpb wlandry/sdpb:2.2.0 mpirun --allow-run-as-root -n 4 sdpb --precision=1024 --procsPerNode=4 -s /usr/local/share/sdpb/input

The results will be in `/my/project/`.  Note that the files may be
owned by root.
