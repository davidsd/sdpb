# Docker

[Docker](https://www.docker.com/) (or its alternative, [Podman](https://podman.io))
is a relatively easy way to run binaries without having to build SDPB.
This can be a convenient way to run SDPB on your laptop or desktop.

For HPC centers, most of them support Singularity rather
than Docker, so you should read the guide for
[Singularity](Singularity.md).

There are instructions at the [Docker website](https://www.docker.com/get-started) on how to install it on Windows, Mac,
or Linux.

## Get Docker image

### Download from Docker Hub

Docker image built from the `master` branch (updated regularly)
can be downloaded from [Docker Hub](https://hub.docker.com/r/bootstrapcollaboration/sdpb/tags) as

    docker pull bootstrapcollaboration/sdpb:master

You can also download image for a specific release, e.g., 2.7.0:

    docker pull bootstrapcollaboration/sdpb:2.7.0

The list of all available tags can be found in
[bootstrapcollaboration/sdpb](https://hub.docker.com/r/bootstrapcollaboration/sdpb/tags) repo on Docker Hub.
Images for SDPB 2.5.1 and earlier can be downloaded from [wlandry/sdpb](https://hub.docker.com/r/wlandry/sdpb/tags)
repo.

### Build from sources

You can also [build](https://docs.docker.com/engine/reference/commandline/build/) Docker image from sources by yourself
using [Dockerfile](../Dockerfile):

    docker build . --tag sdpb:master
    docker run sdpb:master sdpb --help

Our [Dockerfile](../Dockerfile) also contains separate `test` target which allows
to [run built-in tests](../test/run_all_tests.sh) inside Docker containter:

    docker build . --tag sdpb-test --target test
    docker run sdpb-test ./test/run_all_tests.sh

## Run Docker image

Suppose you have an input file `/my/project/input.json`. To use this
file, and be able to write output, we will make the `/my/project`
directory visible to the docker image in the location `/usr/local/share/sdpb`.

The structure of the [docker run](https://docs.docker.com/engine/reference/commandline/run/) command is

    docker run <options> <image> <command>
    
In this case, `<options>` will be used to mount your directory in a
place that docker will see it.

    -v /my/project/:/usr/local/share/sdpb/

`<image>` is the image name, e.g. `bootstrapcollaboration/sdpb:master`

`<command>` is the command that you would normally use to run the SDPB
commands (see [Usage.md](Usage.md)).  The directory containing the
input file is mounted as `/usr/local/share/sdpb`. So we first run `pmp2sdp` to
convert from json

    mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/input.json -o /usr/local/share/sdpb/sdp
    
`mpirun` runs as root inside the docker container.  Running `mpirun` as
root is normally dangerous, but it is safe to do so inside the
container. To allow `mpirun` to run as root, we add the option
`--allow-run-as-root`. This uses 4 cores when running pmp2sdp. You
can change that number to match your own machine.  Putting it all
together on a single line

    docker run -v /my/project/:/usr/local/share/sdpb bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/input.json -o /usr/local/share/sdpb/sdp

Running this command will create directory `/my/project/sdp`.
To search for primal-dual solutions

    docker run -v /my/project/:/usr/local/share/sdpb bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp

The results will be in `/my/project/`.

Note that the newly created files may be owned by root.
If you cannot remove them outside the container, run `rm` from the container, e.g.:

    docker run -v /my/project/:/usr/local/share/sdpb bootstrapcollaboration/sdpb:master rm /usr/local/share/sdpb/sdp

# Podman

Instead of Docker, one can also use [Podman](https://podman.io), which is compatible with Docker images and uses the
same syntax (the only difference is the additional `docker.io/` prefix):

    podman pull docker.io/bootstrapcollaboration/sdpb:master
    podman run -v /my/project/:/usr/local/share/sdpb bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/input.json -o /usr/local/share/sdpb/sdp