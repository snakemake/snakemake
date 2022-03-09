var runtimes_spec = {
    "$schema": "https://vega.github.io/schema/vega-lite/v3.json",
    "description": "Runtimes of jobs.",
    "data": { "values": [{ 'rule': 'a', 'runtime': 2.0280115604400635 }, { 'rule': 'c', 'runtime': 3.028017282485962 }, { 'rule': 'c', 'runtime': 3.020017147064209 }, { 'rule': 'c', 'runtime': 2.0160114765167236 }, { 'rule': 'c', 'runtime': 3.016017198562622 }, { 'rule': 'c', 'runtime': 1.0240058898925781 }, { 'rule': 'c', 'runtime': 3.016016960144043 }, { 'rule': 'c', 'runtime': 1.0160057544708252 }, { 'rule': 'c', 'runtime': 2.0120115280151367 }, { 'rule': 'c', 'runtime': 3.024017095565796 }, { 'rule': 'c', 'runtime': 1.0160057544708252 }, { 'rule': 'd', 'runtime': 1.0520060062408447 }, { 'rule': 'e', 'runtime': 0.012000083923339844 }] },
    "mark": "point",
    "encoding": {
        "x": {
            "field": "runtime", "type": "quantitative",
            "axis": { "title": "runtime [s]", "labelAngle": -90 },
            "scale": { "type": "log" }
        },
        "y": { "field": "rule", "type": "nominal" },
        "color": { "value": "#007bff" }
    }
};