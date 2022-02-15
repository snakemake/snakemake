#!/bin/bash
set -x
echo `date`
echo cancel $* >>scancel.txt

list_descendants ()
{
  local children=$(ps -o pid= --ppid "$1")

  for pid in $children
  do
    list_descendants "$pid"
  done

  echo "$children"
}

for x in $*; do
    kill $(list_descendants $x)
done
