import os

LOGO_URL = os.environ["DUMMY_SERVER_URL"]

storage:
    provider="http",
    keep_local=True

storage http_ret:
    provider="http",
    keep_local=True,
    retrieve=False,


rule retrieve_default:
    input: storage.http(LOGO_URL)
    output: "retrieve_default.flag"
    shell: "touch {output}"

rule retrieve_true:
    input: storage.http(LOGO_URL, retrieve=True)
    output: "retrieve_true.flag"
    shell: "touch {output}"

rule retrieve_false:
    input: storage.http(LOGO_URL, retrieve=False)
    output: "retrieve_false.flag"
    shell: "touch {output}"

rule retrieve_false_directive:
    input: storage.http_ret(LOGO_URL)
    output: "retrieve_false_directive.flag"
    shell: "touch {output}"
