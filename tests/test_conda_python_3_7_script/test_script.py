import platform
import PIL

with open('version.txt', 'w') as f:
    f.write(platform.python_version())
