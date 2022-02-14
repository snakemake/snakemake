#!/bin/bash
echo mcancel >>scancel.txt
echo mcancel $*
kill $* &>/dev/null
