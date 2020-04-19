# Genius TCAD Docker Usage

Genius-TCAD uses a number of old versions of libraries, this makes it increasingly difficult to compile and run on the latest Linux distributions.  To work around this a Docker image of an older Ubuntu 14.04 distribution is used.

## Installation

### Docker Installation

Search for instructions on how to install Docker on your distribution.

### Container Installation

To install the Docker container navigate to the folder containing this document and execute the following command:

    ./build.sh

The software will be automatically downloaded and installed - this will take some time.


## Usage

To interact with the installed container execute the following command:

    ./bash.sh

Which will allow you to navigate and execute software within the container.

To move files to and from the container use the Docker command (from a host terminal):

    docker cp INST_NAME:/container/file/location /host/file/location

Use the ``docker ps`` to file the instance name of running containers.

