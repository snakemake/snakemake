# Test Remote iRODS

These are the instructions for testing Snakemake remote iRODS support locally
on your computer.

## Prerequisites

This requires the latest Docker installation (a.k.a. not the one from the
Ubuntu 16.04 repository). Follow the instructions on
https://docs.docker.com/engine/installation/linux/docker-ce/ubuntu/

## Build and run the Docker container

```
make run
```

## Stop and delete the Docker container

```
make stop
```

## Run the test

```
snakemake
```

## Touch the input file (for a new test)

```
make touch
```

## Example

```
make run
snakemake
snakemake
# nothing to do here
make touch
snakemake
```
