#!/bin/bash
echo cancel >>scancel.txt
kill $* &>/dev/null
