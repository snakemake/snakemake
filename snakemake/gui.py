__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2019, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import json
import os
import threading

from flask import Flask, render_template, request

from snakemake.common import __version__

LOCK = threading.Lock()

app = Flask("snakemake", template_folder=os.path.dirname(__file__))
# app.debug=True
app.extensions = {
    "dag": None,
    "run_snakemake": None,
    "progress": "",
    "log": [],
    "status": {"running": False},
    "args": None,
    "targets": [],
    "rule_info": [],
    "resources": [],
}


def register(run_snakemake, args):
    app.extensions["run_snakemake"] = run_snakemake
    app.extensions["args"] = dict(
        targets=args.target,
        cluster=args.cluster,
        workdir=args.directory,
        touch=args.touch,
        forcetargets=args.force,
        forceall=args.forceall,
        forcerun=args.forcerun,
        prioritytargets=args.prioritize,
        stats=args.stats,
        keepgoing=args.keep_going,
        jobname=args.jobname,
        immediate_submit=args.immediate_submit,
        ignore_ambiguity=args.allow_ambiguity,
        lock=not args.nolock,
        force_incomplete=args.rerun_incomplete,
        ignore_incomplete=args.ignore_incomplete,
        jobscript=args.jobscript,
        notemp=args.notemp,
        latency_wait=args.latency_wait,
    )

    target_rules = []

    def log_handler(msg):
        if msg["level"] == "rule_info":
            target_rules.append(msg["name"])

    run_snakemake(list_target_rules=True, log_handler=log_handler)
    for target in args.target:
        target_rules.remove(target)
    app.extensions["targets"] = args.target + target_rules

    resources = []

    def log_handler(msg):
        if msg["level"] == "info":
            resources.append(msg["msg"])

    run_snakemake(list_resources=True, log_handler=log_handler)
    app.extensions["resources"] = resources
    app.extensions["snakefilepath"] = os.path.abspath(args.snakefile)


def run_snakemake(**kwargs):
    args = dict(app.extensions["args"])
    args.update(kwargs)
    app.extensions["run_snakemake"](**args)


@app.route("/")
def index():
    args = app.extensions["args"]
    return render_template(
        "gui.html",
        targets=app.extensions["targets"],
        cores_label="Nodes" if args["cluster"] else "Cores",
        resources=app.extensions["resources"],
        snakefilepath=app.extensions["snakefilepath"],
        version=__version__,
        node_width=15,
        node_padding=10,
    )


@app.route("/dag")
def dag():
    if app.extensions["dag"] is None:

        def record(msg):
            if msg["level"] == "d3dag":
                app.extensions["dag"] = msg
            elif msg["level"] in ("error", "info"):
                app.extensions["log"].append(msg)

        run_snakemake(printd3dag=True, log_handler=record)
    return json.dumps(app.extensions["dag"])


@app.route("/log/<int:id>")
def log(id):
    log = app.extensions["log"][id:]
    return json.dumps(log)


@app.route("/progress")
def progress():
    return json.dumps(app.extensions["progress"])


def _run(dryrun=False):
    def log_handler(msg):
        level = msg["level"]
        if level == "progress":
            app.extensions["progress"] = msg
        elif level in ("info", "error", "job_info", "job_finished"):
            app.extensions["log"].append(msg)

    with LOCK:
        app.extensions["status"]["running"] = True
    run_snakemake(log_handler=log_handler, dryrun=dryrun)
    with LOCK:
        app.extensions["status"]["running"] = False
    return ""


@app.route("/run")
def run():
    _run()


@app.route("/dryrun")
def dryrun():
    _run(dryrun=True)


@app.route("/status")
def status():
    with LOCK:
        return json.dumps(app.extensions["status"])


@app.route("/targets")
def targets():
    return json.dumps(app.extensions["targets"])


@app.route("/get_args")
def get_args():
    return json.dumps(app.extensions["args"])


@app.route("/set_args", methods=["POST"])
def set_args():
    app.extensions["args"].update(
        {name: value for name, value in request.form.items() if not name.endswith("[]")}
    )
    targets = request.form.getlist("targets[]")
    if targets != app.extensions["args"]["targets"]:
        app.extensions["dag"] = None
    app.extensions["args"]["targets"] = targets
    return ""
