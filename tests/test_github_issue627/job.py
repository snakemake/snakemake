#!/usr/bin/env python

import time

file = open("job.txt", "w")

for i in range(300):
    file.write(str(i) + "\n")
    file.flush()
    time.sleep(1)
