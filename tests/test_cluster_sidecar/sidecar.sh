#!/bin/bash

set -ex

echo "FIRST_LINE"
echo "sidecar started" > sidecar.txt
sleep infinity &
pid=$!

catch()
{
    set -x
    kill -TERM $pid || true
    echo "sidecar stopped" >> sidecar.txt
    exit 0
}

trap catch SIGTERM SIGINT

wait
