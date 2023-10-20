import platform

with open('version.txt', 'w') as f:
    f.write(platform.python_version())
