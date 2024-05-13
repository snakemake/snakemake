import json
import pytest
from snakemake.utils import validate
import pandas as pd

CONFIG_SCHEMA = """$schema: "https://json-schema.org/draft/2020-12/schema#"
description: Configuration schema
type: object
properties:
  param:
    type: object
    description: parameter settings
    # Necessary in case config is empty!
    default: {}
    properties:
      foo:
        type: string
        default: bar
"""

BAR_SCHEMA = """definitions:
  bar:
    type: string
    description: bar entry
    default: foo
"""

BAR_JSON_SCHEMA = {
    "definitions": {
        "jsonbar": {"type": "string", "description": "bar entry", "default": "foo"}
    }
}

DF_SCHEMA = """$schema: "https://json-schema.org/draft/2020-12/schema#"
description: an entry in the sample sheet
properties:
  sample:
    type: string
    description: sample name/identifier
  condition:
    type: string
    description: sample condition that will be compared during differential expression analysis (e.g. a treatment, a tissue time, a disease)
  case:
    type: boolean
    default: true
    description: boolean that indicates if sample is case or control
  date:
    type: string
    format: date-time
    description: collection date
    default: "2018-10-10"

required:
  - sample
  - condition
"""


@pytest.fixture
def schemadir(tmpdir):
    p = tmpdir.mkdir("schema")
    return p


@pytest.fixture
def bar_schema(schemadir):
    p = schemadir.join("bar.schema.yaml")
    p.write(BAR_SCHEMA)
    return p


@pytest.fixture
def json_bar_schema(schemadir):
    p = schemadir.join("bar.schema.json")
    p.write(json.dumps(BAR_JSON_SCHEMA))
    return p


@pytest.fixture
def df_schema(schemadir):
    p = schemadir.join("df.schema.yaml")
    p.write(DF_SCHEMA)
    return p


@pytest.fixture
def config_schema(schemadir):
    p = schemadir.join("config.schema.yaml")
    p.write(CONFIG_SCHEMA)
    return p


@pytest.fixture
def config_schema_ref(schemadir, bar_schema, json_bar_schema):
    p = schemadir.join("config.ref.schema.yaml")
    p.write(
        CONFIG_SCHEMA
        + "\n".join(
            [
                "      bar:",
                '        default: "yaml"',
                '        $ref: "{bar}"'.format(
                    bar=str(bar_schema) + "#/definitions/bar"
                ),
                "      jsonbar:",
                '        default: "json"',
                '        $ref: "{bar}"'.format(
                    bar=str(json_bar_schema) + "#/definitions/jsonbar"
                ),
                "",
            ]
        )
    )
    return p


def test_config(config_schema):
    config = {}
    validate(config, str(config_schema), False)
    assert config == {}
    validate(config, str(config_schema))
    assert dict(config) == {"param": {"foo": "bar"}}


def test_config_ref(config_schema_ref):
    config = {}
    validate(config, str(config_schema_ref))
    assert config["param"]["foo"] == "bar"
    assert config["param"]["bar"] == "yaml"
    assert config["param"]["jsonbar"] == "json"
    # Make sure regular validator works
    config["param"]["bar"] = 1
    config["param"]["jsonbar"] = 2
    from snakemake.exceptions import WorkflowError

    with pytest.raises(WorkflowError):
        validate(config, str(config_schema_ref), False)


def test_dataframe(df_schema):
    df = pd.DataFrame([{"sample": "foo", "condition": "bar"}])
    validate(df, str(df_schema), False)
    assert sorted(df.columns) == sorted(["sample", "condition"])
    validate(df, str(df_schema))
    assert sorted(df.columns) == sorted(["sample", "condition", "case", "date"])
    assert df.case.loc[0]
