import json
import os
import threading

from flask import Flask, render_template

LOCK = threading.Lock()

app = Flask("Snakemake", template_folder=os.path.dirname(__file__))
app.extensions = {
    "dag": None, "run_snakemake": None, "progress": "", "log": [],
    "status": {"running": False}, "args": None, "targets": [], "rule_info": []}


def register(run_snakemake, args):
    app.extensions["run_snakemake"] = run_snakemake
    app.extensions["args"] = {
        "targets": args.target,
        "cluster": args.cluster
    }
    run_snakemake(list_target_rules=True)
    target_rules = [msg["name"] for msg in app.extensions["rule_info"]]
    for target in args.target:
        target_rules.remove(target)
    app.extensions["targets"] = args.target + target_rules


def run_snakemake(**kwargs):
    args = dict(app.extensions["args"])
    args.update(kwargs)
    app.extensions["run_snakemake"](**args)


def log_handler(msg):
    level = msg["level"]
    if level == "d3dag":
        app.extensions["dag"] = msg
    elif level == "progress":
        app.extensions["progress"] = msg
    elif level == "rule_info":
        app.extensions["rule_info"].append(msg)
    elif level in ("info", "error", "job_info"):
        app.extensions["log"].append(msg)


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
        run_snakemake(printd3dag=True)
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
    with LOCK:
        app.extensions["status"]["running"] = True
    run_snakemake()
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
