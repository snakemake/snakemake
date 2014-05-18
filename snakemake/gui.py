import json
import os
import threading

from flask import Flask, render_template

LOCK = threading.Lock()

app = Flask("Snakemake", template_folder=os.path.dirname(__file__))
app.extensions = {
    "dag": None, "run_snakemake": None, "progress": "", "log": [],
    "status": {"running": False}, "args": None, "targets": [], "rule_info": [],
    "resources": []}


def register(run_snakemake, args):
    app.extensions["run_snakemake"] = run_snakemake
    app.extensions["args"] = {
        "targets": args.target,
        "cluster": args.cluster
    }
    
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
        node_width=15,
        node_padding=10)


@app.route("/dag")
def dag():
    if app.extensions["dag"] is None:
        def record(msg):
            if msg["level"] == "d3dag":
                app.extensions["dag"] = msg
        run_snakemake(printd3dag=True, log_handler=record)
    return json.dumps(app.extensions["dag"])


@app.route("/log/<int:id>")
def log(id):
    log = app.extensions["log"][id:]
    import json
    return json.dumps(log)


@app.route("/progress")
def progress():
    return json.dumps(app.extensions["progress"])


@app.route("/run")
def run():
    def log_handler(msg):
        level = msg["level"]
        if level == "progress":
            app.extensions["progress"] = msg
        elif level in ("info", "error", "job_info"):
            app.extensions["log"].append(msg)

    with LOCK:
        app.extensions["status"]["running"] = True
    run_snakemake(log_handler=log_handler)
    with LOCK:
        app.extensions["status"]["running"] = False
    return ""


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
    app.extensions["args"] = request.form
    return ""
