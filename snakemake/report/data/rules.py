import json


def render_rules(rules):
    return json.dumps(
        {
            rulename: {
                "input": rule.input,
                "output": rule.output,
                "conda_env": rule.conda_env,
                "container_img_url": rule.container_img_url,
                "code": rule.code,
                "n_jobs": rule.n_jobs,
            }
            for rulename, rule in rules.items()
        }
    )
