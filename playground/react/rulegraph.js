var rulegraph_spec = {
    "$schema": "https://vega.github.io/schema/vega/v5.json",
    "padding": 0,

    "signals": [
        { "name": "cx", "update": "width / 2" },
        { "name": "cy", "update": "height / 2" }
    ],

    "data": [
        {
            "name": "node-data",
            "values": [{ 'rule': 'c', 'fx': 10, 'fy': 0 }, { 'rule': 'd', 'fx': 10, 'fy': 50 }, { 'rule': 'e', 'fx': 10, 'fy': 100 }, { 'rule': 'a', 'fx': 10, 'fy': 150 }, { 'rule': 'b', 'fx': 10, 'fy': 200 }, { 'rule': 'all', 'fx': 10, 'fy': 250 }]
        },
        {
            "name": "link-data",
            "values": [{ 'target': 3, 'source': 0, 'value': 1 }, { 'target': 3, 'source': 0, 'value': 1 }, { 'target': 3, 'source': 0, 'value': 1 }, { 'target': 3, 'source': 0, 'value': 1 }, { 'target': 3, 'source': 0, 'value': 1 }, { 'target': 3, 'source': 0, 'value': 1 }, { 'target': 3, 'source': 0, 'value': 1 }, { 'target': 3, 'source': 0, 'value': 1 }, { 'target': 3, 'source': 0, 'value': 1 }, { 'target': 3, 'source': 0, 'value': 1 }, { 'target': 4, 'source': 0, 'value': 1 }, { 'target': 4, 'source': 0, 'value': 1 }, { 'target': 4, 'source': 0, 'value': 1 }, { 'target': 4, 'source': 0, 'value': 1 }, { 'target': 4, 'source': 0, 'value': 1 }, { 'target': 4, 'source': 0, 'value': 1 }, { 'target': 4, 'source': 0, 'value': 1 }, { 'target': 4, 'source': 0, 'value': 1 }, { 'target': 4, 'source': 0, 'value': 1 }, { 'target': 4, 'source': 0, 'value': 1 }, { 'target': 5, 'source': 3, 'value': 1 }, { 'target': 5, 'source': 1, 'value': 1 }, { 'target': 5, 'source': 2, 'value': 1 }]
        },
        {
            "name": "link-data-direct",
            "values": [{ 'target': 5, 'source': 4, 'value': 1 }]
        }
    ],

    "scales": [
        {
            "name": "color",
            "type": "ordinal",
            "range": { "scheme": "category20c" }
        },
        {
            "name": "x",
            "type": "linear"
        },
        {
            "name": "y",
            "type": "linear"
        }
    ],

    "marks": [
        {
            "name": "nodes",
            "type": "symbol",
            "zindex": 1,
            "from": { "data": "node-data" },
            "encode": {
                "enter": {
                    "fill": { "scale": "color", "field": "rule" },
                    "x": { "field": "fx", "scale": "x" },
                    "y": { "field": "fy", "scale": "y" },
                    "tooltip": { "value": "Click to show rule details." }
                },
                "update": {
                    "size": { "value": 70 }
                },
                "hover": {
                    "size": { "value": 140 }
                }
            },

            "transform": [
                {
                    "type": "force",
                    "iterations": 1,
                    "static": true,
                    "forces": [
                        {
                            "force": "link",
                            "links": "link-data"
                        }
                    ]
                },
                {
                    "type": "force",
                    "iterations": 1,
                    "static": true,
                    "forces": [
                        {
                            "force": "link",
                            "links": "link-data-direct"
                        }
                    ]
                }
            ]
        },
        {
            "name": "labels",
            "type": "text",
            "zindex": 2,
            "from": { "data": "node-data" },
            "encode": {
                "enter": {
                    "fill": { "value": "black" },
                    "fontWeight": { "value": "normal" },
                    "text": { "field": "rule" },
                    "x": { "field": "fx", "scale": "x" },
                    "y": { "field": "fy", "scale": "y" },
                    "dx": { "value": -5 },
                    "dy": { "value": -5 },
                    "align": { "value": "right" }
                }
            }
        },
        {
            "type": "path",
            "from": { "data": "link-data-direct" },
            "interactive": false,
            "encode": {
                "update": {
                    "stroke": { "value": "#ccc" },
                    "strokeWidth": { "value": 1.0 }
                }
            },
            "transform": [
                {
                    "type": "linkpath", "shape": "diagonal",
                    "sourceX": "datum.source.x", "sourceY": "datum.source.y",
                    "targetX": "datum.target.x", "targetY": "datum.target.y"
                }
            ]
        },
        {
            "type": "path",
            "from": { "data": "link-data" },
            "interactive": false,
            "encode": {
                "update": {
                    "stroke": { "value": "#ccc" },
                    "strokeWidth": { "value": 1.0 }
                }
            },
            "transform": [
                {
                    "type": "linkpath", "shape": "curve", "orient": "horizontal",
                    "sourceX": "datum.source.x", "sourceY": "datum.source.y",
                    "targetX": "datum.target.x", "targetY": "datum.target.y"
                }
            ]
        }
    ]
};