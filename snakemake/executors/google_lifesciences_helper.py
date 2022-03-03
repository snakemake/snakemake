#!/usr/bin/env python

# This is a helper script for the Google Life Sciences instance to be able to:
# 1. download a blob from storage, which is required at the onset of the Snakemake
#     gls.py download <bucket> <source> <destination>
# workflow step to obtain the working directory.
# 2. Upload logs back to storage (or some specified directory of files)
#    gls.py save <bucket> <source-dir> <destination-dir>
#    gls.py save <bucket> /google/logs/output source/logs

import argparse
import datetime

from google.cloud import storage
from glob import glob
import sys
import os


def download_blob(bucket_name, source_blob_name, destination_file_name):
    """Downloads a blob from the bucket."""
    storage_client = storage.Client()
    bucket = storage_client.get_bucket(bucket_name)
    blob = bucket.blob(source_blob_name)

    blob.download_to_filename(destination_file_name)

    print("Blob {} downloaded to {}.".format(source_blob_name, destination_file_name))


def save_files(bucket_name, source_path, destination_path):
    """given a directory path, save all files recursively to storage"""
    storage_client = storage.Client()
    bucket = storage_client.get_bucket(bucket_name)

    # destination path should be stripped of path indicators too
    bucket_name = bucket_name.strip("/")
    destination_path = destination_path.strip("/")

    # These are fullpaths
    filenames = get_source_files(source_path)
    print("\nThe following files will be uploaded: %s" % "\n".join(filenames))

    if not filenames:
        print("Did not find any filenames under %s" % source_path)

    # Do the upload!
    for filename in filenames:

        # The relative path of the filename from the source path
        relative_path = filename.replace(source_path, "", 1).strip("/")

        # The path in storage includes relative path from destination_path
        storage_path = os.path.join(destination_path, relative_path)
        full_path = os.path.join(bucket_name, storage_path)
        print(
            "{filename} -> {full_path}".format(filename=filename, full_path=full_path)
        )

        # Get the blob
        blob = bucket.blob(storage_path)
        if not blob.exists():
            print("Uploading %s to %s" % (filename, full_path))
            blob.upload_from_filename(filename)


def get_source_files(source_path):
    """Given a directory, return a listing of files to upload"""
    filenames = []
    if not os.path.exists(source_path):
        print("%s does not exist!" % source_path)
        sys.exit(0)

    for x in os.walk(source_path):
        for name in glob(os.path.join(x[0], "*")):
            if not os.path.isdir(name):
                filenames.append(name)
    return filenames


def add_ending_slash(filename):
    """Since we want to replace based on having an ending slash, ensure it's there"""
    if not filename.endswith("/"):
        filename = "%s/" % filename
    return filename


def blob_commands(args):
    if args.command == "download":
        download_blob(
            args.bucket_name, args.source_blob_name, args.destination_file_name
        )
    elif args.command == "save":
        save_files(args.bucket_name, args.source_path, args.destination_path)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    subparsers = parser.add_subparsers(dest="command")

    # Download file from storage
    download_parser = subparsers.add_parser("download", help=download_blob.__doc__)
    download_parser.add_argument("bucket_name", help="Your cloud storage bucket.")
    download_parser.add_argument("source_blob_name")
    download_parser.add_argument("destination_file_name")

    # Save logs to storage
    save_parser = subparsers.add_parser("save", help=save_files.__doc__)
    save_parser.add_argument("bucket_name", help="Your cloud storage bucket.")
    save_parser.add_argument("source_path")
    save_parser.add_argument("destination_path")

    args = parser.parse_args()
    blob_commands(args)


if __name__ == "__main__":
    main()
