# Test Remote iRODS

## Run Docker

To set up the iRODS server for the test, make sure that the `Dockerfile` is
available and run the following commands:

```
docker build -t irods-server .
docker run \
    -d \
    -p 1247:1247 \
    --name provider \
    irods-server \
        -i run_irods
docker exec -u irods provider iput /incoming/infile
```

To touch the input file for a new timestamp, issue:

```
docker exec -u irods provider iput /incoming/infile -f
```

