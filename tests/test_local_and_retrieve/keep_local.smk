import os

LOGO_URL = os.environ["DUMMY_SERVER_URL"]

storage:
    provider="http"

storage http_local:
    provider="http",
    keep_local=True

rule keep_local_default:
    input: storage.http(LOGO_URL)
    output: "keep_local_default.flag"
    shell: "touch {output}"

rule keep_local_true:
    input: storage.http(LOGO_URL, keep_local=True)
    output: "keep_local_true.flag"
    shell: "touch {output}"

rule keep_local_false:
    input: storage.http(LOGO_URL, keep_local=False)
    output: "keep_local_false.flag"
    shell: "touch {output}"

rule keep_local_true_directive:
    input: storage.http_local(LOGO_URL)
    output: "keep_local_true_directive.flag"
    shell: "touch {output}"
