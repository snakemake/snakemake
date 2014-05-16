from flask import Flask, jsonify, render_template
import os


app = Flask("Snakemake", template_folder=os.path.dirname(__file__))
app.extensions = {}


def register(run_snakemake):
    app.extensions["run_snakemake"] = run_snakemake


def run_snakemake(**kwargs):
    app.extensions["run_snakemake"](**kwargs)


def log_handler(msg):
    print("log")
    if msg["level"] == "d3dag":
        print(msg)
        app.extensions["dag"] = msg


@app.route("/")
def index():
    return render_template("gui.html", dag_width=900, dag_height=500, node_width=15, node_padding=10)


@app.route("/dag.json")
def dag():
    run_snakemake(printd3dag=True)
    return jsonify(app.extensions["dag"])
