import os
import sys
import subprocess
import requests_mock
import json

sys.path.insert(0, os.path.dirname(__file__))

from common import *


TES_URL = "http://localhost:8000"

TEST_POST_RESPONSE = {"id": "id_1"}

TEST_TASK = {"id": "id_1", "state": "COMPLETE"}


def _validate_task(task):
    print("\n>>>> _validate_task", file=sys.stderr)
    required_keys = ["inputs", "outputs", "executors"]
    return all(i in task.keys() for i in required_keys)


def _post_task(request, context):
    print("\n>>>> _post_task", file=sys.stderr)
    task = json.loads(request.body)
    print(task, file=sys.stderr)
    if _validate_task(task):
        context.status_code = 200
        return TEST_POST_RESPONSE
    else:
        context.status_code = 400
        return None


def _get_task(request, context):
    print("\n>>>> _get_task", file=sys.stderr)
    context.status_code = 200
    return TEST_TASK


def test_tes(requests_mock):
    requests_mock.register_uri("POST", "{}/v1/tasks".format(TES_URL), json=_post_task)
    requests_mock.register_uri(
        "GET", "{}/v1/tasks/id_1".format(TES_URL), json=_get_task
    )
    workdir = dpath("test_tes")
    run(workdir, snakefile="Snakefile", tes=TES_URL, no_tmpdir=True, cleanup=False)
