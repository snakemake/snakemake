import json


def render_timeline(timeline):
    return json.dumps(
        {
            "$schema": "https://vega.github.io/schema/vega-lite/v3.json",
            "description": "Timeline of jobs.",
            "data": {"values": timeline},
            "mark": "point",
            "encoding": {
                "x": {
                    "field": "endtime",
                    "type": "temporal",
                    "timeUnit": "yearmonthdatehoursminutes",
                    "axis": {"labelAngle": -90, "title": "end date"},
                },
                "y": {"field": "rule", "type": "nominal"},
                "color": {"value": "#007bff"},
            },
        }
    )
