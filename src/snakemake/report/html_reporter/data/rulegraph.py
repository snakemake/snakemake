import json


def render_rulegraph(nodes, links, links_direct):
    return json.dumps(
        {
            "$schema": "https://vega.github.io/schema/vega/v5.json",
            "padding": 0,
            "signals": [
                {"name": "cx", "update": "width / 2"},
                {"name": "cy", "update": "height / 2"},
            ],
            "data": [
                {"name": "node-data", "values": nodes},
                {"name": "link-data", "values": links},
                {"name": "link-data-direct", "values": links_direct},
            ],
            "scales": [
                {
                    "name": "color",
                    "type": "ordinal",
                    "range": {"scheme": "category20c"},
                },
                {"name": "x", "type": "linear"},
                {"name": "y", "type": "linear"},
            ],
            "marks": [
                {
                    "name": "nodes",
                    "type": "symbol",
                    "zindex": 1,
                    "from": {"data": "node-data"},
                    "encode": {
                        "enter": {
                            "fill": {"scale": "color", "field": "rule"},
                            "x": {"field": "fx", "scale": "x"},
                            "y": {"field": "fy", "scale": "y"},
                            "tooltip": {"value": "Click to show rule details."},
                        },
                        "update": {"size": {"value": 70}},
                        "hover": {"size": {"value": 140}},
                    },
                    "transform": [
                        {
                            "type": "force",
                            "iterations": 1,
                            "static": True,
                            "forces": [{"force": "link", "links": "link-data"}],
                        },
                        {
                            "type": "force",
                            "iterations": 1,
                            "static": True,
                            "forces": [{"force": "link", "links": "link-data-direct"}],
                        },
                    ],
                },
                {
                    "name": "labels",
                    "type": "text",
                    "zindex": 2,
                    "from": {"data": "node-data"},
                    "encode": {
                        "enter": {
                            "fill": {"value": "black"},
                            "fontWeight": {"value": "normal"},
                            "text": {"field": "rule"},
                            "x": {"field": "fx", "scale": "x"},
                            "y": {"field": "fy", "scale": "y"},
                            "dx": {"value": -5},
                            "dy": {"value": -5},
                            "align": {"value": "right"},
                        }
                    },
                },
                {
                    "type": "path",
                    "from": {"data": "link-data-direct"},
                    "interactive": True,
                    "encode": {
                        "update": {
                            "stroke": {"value": "#ccc"},
                            "strokeWidth": {"value": 1.0},
                        },
                        "hover": {"strokeWidth": {"value": 4.0}},
                    },
                    "transform": [
                        {
                            "type": "linkpath",
                            "shape": "diagonal",
                            "sourceX": "datum.source.x",
                            "sourceY": "datum.source.y",
                            "targetX": "datum.target.x",
                            "targetY": "datum.target.y",
                        }
                    ],
                },
                {
                    "type": "path",
                    "from": {"data": "link-data"},
                    "interactive": True,
                    "encode": {
                        "enter": {
                            "tooltip": {
                                "signal": "{\"from rule\": datum['sourcerule'], \"to rule\": datum['targetrule']}"
                            }
                        },
                        "update": {
                            "stroke": {"value": "#ccc"},
                            "strokeWidth": {"value": 1.0},
                        },
                        "hover": {"strokeWidth": {"value": 4.0}},
                    },
                    "transform": [
                        {
                            "type": "linkpath",
                            "shape": "curve",
                            "orient": "horizontal",
                            "sourceX": "datum.source.x",
                            "sourceY": "datum.source.y",
                            "targetX": "datum.target.x",
                            "targetY": "datum.target.y",
                        }
                    ],
                },
            ],
        }
    )
