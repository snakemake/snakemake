import pandas as pd

with open("pandas_version.txt", "w") as f:
    f.write(pd.__version__)
