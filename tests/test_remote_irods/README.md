# Test Remote iRODS

These are the instructions for testing Snakemake remote iRODS support locally
on your computer.

## Prerequisites

1. The latest Docker installation (a.k.a. not the one from the Ubuntu 16.04
repository). Follow the instructions on
https://docs.docker.com/engine/installation/linux/docker-ce/ubuntu/
2. The irods icommand tools as described in https://packages.irods.org/ to
setup the repository and `apt-get install irods-icommands` for installation.
The password file has to be generated with `iinit`. A valid environment file
is located in `setup-data/irods_environment.json` (in this case the
authentication file is expected in `~/.irods/.irodsA`). This is necessary,
because the obfuscation of the password uses the uid, so the `.irodsA` file
can't be shipped.
3. Snakemake

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
snakemake -s Snakefile.local
```

## Touch the input file (for a new test)

```
make touch
```

## Show iRODS content

```
make ls
```

## Example

```
make run
make ls
snakemake -s Snakefile.local
make ls
snakemake -s Snakefile.local
# nothing to do here
make touch
snakemake -s Snakefile.local
```
