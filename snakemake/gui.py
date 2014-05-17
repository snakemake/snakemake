import json
import os

from flask import Flask, render_template

app = Flask("Snakemake", template_folder=os.path.dirname(__file__))
app.extensions = {"dag": None, "run_snakemake": None, "progress": "", "log": []}


def register(run_snakemake):
    app.extensions["run_snakemake"] = run_snakemake


def run_snakemake(**kwargs):
    app.extensions["run_snakemake"](**kwargs)


def log_handler(msg):
    level = msg["level"]
    if level == "d3dag":
        app.extensions["dag"] = msg
    elif level == "progress":
        app.extensions["progress"] = msg
    elif level in ("info", "error", "job_info"):
        app.extensions["log"].append(msg)


@app.route("/")
def index():
    return render_template("gui.html", dag_width=900, dag_height=500, node_width=15, node_padding=10)


@app.route("/dag")
def dag():
    if app.extensions["dag"] is None:
        run_snakemake(printd3dag=True)
    return json.dumps(app.extensions["dag"])


@app.route("/log/<int:id>")
def log(id):
    log = app.extensions["log"][id:]
    print(log)
    import json
    return json.dumps(log)


@app.route("/progress")
def progress():
    return json.dumps(app.extensions["progress"])


@app.route("/run")
def run():
    try:
        run_snakemake()
    except Exception as e:
        print(e)
    return "test"
