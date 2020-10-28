#!/bin/bash
mkdir -p DistantDir/WithOUTRights
mkdir -p DistantDir/WithRights
touch DistantDir/WithOUTRights/File
touch DistantDir/WithRights/File
chmod 000 DistantDir/WithOUTRights/
mkdir Input
ln -s ../DistantDir/WithOUTRights/File Input/FileWithOUTRights
ln -s ../DistantDir/WithRights/File Input/FileWithRights
