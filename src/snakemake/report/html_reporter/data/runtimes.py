import json


def render_runtimes(runtimes):
    return json.dumps(
        {
            "$schema": "https://vega.github.io/schema/vega-lite/v3.json",
            "description": "Runtimes of jobs.",
            "data": {"values": runtimes},
            "mark": "point",
            "encoding": {
                "x": {
                    "field": "runtime",
                    "type": "quantitative",
                    "axis": {"title": "runtime [s]", "labelAngle": -90},
                    "scale": {"type": "log"},
                },
                "y": {"field": "rule", "type": "nominal"},
                "color": {"value": "#007bff"},
            },
        }
    )
