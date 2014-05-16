from flask import Flask, jsonify, render_template
import os


app = Flask("Snakemake", template_folder=os.path.dirname(__file__))
app.extensions = {"dag": None, "run_snakemake": None, "progress": "", "log": []}


def register(run_snakemake):
    app.extensions["run_snakemake"] = run_snakemake


def run_snakemake(**kwargs):
    app.extensions["run_snakemake"](**kwargs)


def log_handler(msg):
    print("log")
    if msg["level"] == "d3dag":
        print(msg)
        app.extensions["dag"] = msg
    elif msg["level"] == "progress":
        app.extensions["progress"] = msg
    elif msg["level"] != "debug":
        app.extensions["log"].append(msg)


@app.route("/")
def index():
    return render_template("gui.html", dag_width=900, dag_height=500, node_width=15, node_padding=10)


@app.route("/dag")
def dag():
    run_snakemake(printd3dag=True)
    return jsonify(app.extensions["dag"])


@app.route("/log/<int:id>")
def log(id):
    log = app.extensions["log"]
    print(log)
    if id < len(log):
        return jsonify(log[id])
    return ""


@app.route("/progress")
def progress():
    return jsonify(app.extensions["progress"])


@app.route("/run")
def run():
    run_snakemake()
