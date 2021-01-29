#!/bin/bash

PID=$(bash -c "setsid $1 > /dev/null 2>&1 & echo \"\$!\"")
echo "$PID"

bash -c "sleep 2; kill -SIGKILL \"$PID\"" > /dev/null 2>&1 &