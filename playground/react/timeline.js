var timeline_spec = {
    "$schema": "https://vega.github.io/schema/vega-lite/v3.json",
    "description": "Timeline of jobs.",
    "data": {
        "values": [{ 'rule': 'a', 'starttime': '2020-03-27T14:58:42.460061', 'endtime': '2020-03-27T14:58:44.488072' }, { 'rule': 'c', 'starttime': '2020-03-27T14:58:39.404043', 'endtime': '2020-03-27T14:58:42.432060' }, { 'rule': 'c', 'starttime': '2020-03-27T14:58:38.368037', 'endtime': '2020-03-27T14:58:41.388055' }, { 'rule': 'c', 'starttime': '2020-03-27T14:58:38.384037', 'endtime': '2020-03-27T14:58:40.400049' }, { 'rule': 'c', 'starttime': '2020-03-27T14:58:38.352037', 'endtime': '2020-03-27T14:58:41.368054' }, { 'rule': 'c', 'starttime': '2020-03-27T14:58:39.356043', 'endtime': '2020-03-27T14:58:40.380049' }, { 'rule': 'c', 'starttime': '2020-03-27T14:58:38.344037', 'endtime': '2020-03-27T14:58:41.360054' }, { 'rule': 'c', 'starttime': '2020-03-27T14:58:38.332037', 'endtime': '2020-03-27T14:58:39.348043' }, { 'rule': 'c', 'starttime': '2020-03-27T14:58:38.340037', 'endtime': '2020-03-27T14:58:40.352049' }, { 'rule': 'c', 'starttime': '2020-03-27T14:58:39.424043', 'endtime': '2020-03-27T14:58:42.448061' }, { 'rule': 'c', 'starttime': '2020-03-27T14:58:38.376037', 'endtime': '2020-03-27T14:58:39.392043' }, { 'rule': 'd', 'starttime': '2020-03-27T14:58:38.360037', 'endtime': '2020-03-27T14:58:39.412043' }, { 'rule': 'e', 'starttime': '2020-03-27T14:58:38.356037', 'endtime': '2020-03-27T14:58:38.368037' }]
    },
    "mark": "point",
    "encoding": {
        "x": {
            "field": "endtime", "type": "temporal",
            "timeUnit": "yearmonthdatehoursminutes",
            "axis": {
                "labelAngle": -90,
                "title": "creation date"
            }
        },
        "y": { "field": "rule", "type": "nominal" },
        "color": { "value": "#007bff" }
    }
};