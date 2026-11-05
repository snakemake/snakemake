import pickle
greeting = "hello"


def add_newline(text):
    return text + "\n"


message = pickle.loads(pickle.dumps(add_newline))(greeting)
with open(snakemake.output[0], "w") as f:
    f.write(message)
