import platform
import pillow

with open('version.txt', 'w') as f:
    f.write(platform.python_version())
