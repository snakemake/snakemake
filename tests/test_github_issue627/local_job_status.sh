#!/bin/bash

STATUS=$(ps aux | awk '{print $2}' | grep "$1")
if [ -z "$STATUS" ]; then
  echo "failed"
else
  echo "running"
fi