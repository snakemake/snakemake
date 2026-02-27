import json
import pytest
from snakemake.utils import validate
from snakemake.exceptions import WorkflowError
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

BAR_SCHEMA = """$schema: "https://json-schema.org/draft/2020-12/schema#"
definitions:
  bar:
    type: string
    description: bar entry
    default: foo
"""

BAR_JSON_SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema#",
    "definitions": {
        "jsonbar": {"type": "string", "description": "bar entry", "default": "foo"}
    },
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

RELATIVE_CONFIG_SCHEMA = """$schema: "https://json-schema.org/draft/2020-12/schema#"
type: object
properties:
  bar:
    $ref: "bar.schema.yaml"
"""

NESTED_BAR_SCHEMA = """$schema: "https://json-schema.org/draft/2020-12/schema#"
type: object
properties:
  baz:
    $ref: "baz.schema.yaml"
required: ["baz"]
"""

NESTED_BAZ_SCHEMA = """$schema: "https://json-schema.org/draft/2020-12/schema#"
type: object
properties:
  qux:
    type: integer
    default: 42
required: ["qux"]
"""

REMOTE_ID_CONFIG_SCHEMA = """$schema: "https://json-schema.org/draft/2020-12/schema#"
$id: "https://example.com/config.schema.yaml"
type: object
properties:
  bar:
    $ref: "bar.schema.yaml"
"""

EXTERNAL_BAR_SCHEMA = """$schema: "https://json-schema.org/draft/2020-12/schema#"
type: object
properties:
  foo:
    type: string
required: ["foo"]
"""

FRAGMENT_BAR_SCHEMA = """$schema: "https://json-schema.org/draft/2020-12/schema#"
definitions:
  frag_obj:
    type: object
    properties:
      baz:
        type: string
    required: ["baz"]
"""

FRAGMENT_CONFIG_SCHEMA = """$schema: "https://json-schema.org/draft/2020-12/schema#"
type: object
properties:
  bar:
    $ref: "bar.schema.yaml#/definitions/frag_obj"
required: ["bar"]
"""

ALLOF_CONFIG_SCHEMA = """$schema: "https://json-schema.org/draft/2020-12/schema#"
type: object
properties:
  foo:
    allOf:
      - type: object
        properties:
          bar:
            type: integer
            default: 42
          required: ['bar']
      - type: object
        properties:
          foo:
            type: string
            default: "foo"
        required: ["foo"]
"""

DEFS_CONFIG_SCHEMA = """$schema: "https://json-schema.org/draft/2020-12/schema#"
$defs:
  bar:
    type: object
    properties:
      qux:
        type: string
        default: "qux"
  baz:
    type: object
    properties:
      qux:
        type: integer
        default: 42
    required: ["qux"]
type: object
properties:
  foo:
    anyOf:
      - $ref: "#/$defs/bar"
      - $ref: "#/$defs/baz"
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


@pytest.fixture
def config_schema_relative(schemadir):
    p = schemadir.join("config.relative.schema.yaml")
    p.write(RELATIVE_CONFIG_SCHEMA)
    bar = schemadir.join("bar.schema.yaml")
    bar.write(EXTERNAL_BAR_SCHEMA)
    return p


@pytest.fixture
def config_schema_relative_nested(schemadir):
    p = schemadir.join("config.relative.nested.schema.yaml")
    p.write(RELATIVE_CONFIG_SCHEMA)
    bar = schemadir.join("bar.schema.yaml")
    bar.write(NESTED_BAR_SCHEMA)
    baz = schemadir.join("baz.schema.yaml")
    baz.write(NESTED_BAZ_SCHEMA)
    return p


@pytest.fixture
def config_schema_relative_ref_with_remote_id(schemadir):
    p = schemadir.join("config.relative_ref.remote_id.schema.yaml")
    p.write(REMOTE_ID_CONFIG_SCHEMA)
    bar = schemadir.join("bar.schema.yaml")
    bar.write(EXTERNAL_BAR_SCHEMA)
    return p


@pytest.fixture
def config_schema_relative_ref_with_fragment(schemadir):
    p = schemadir.join("config.relative_ref.fragment.schema.yaml")
    p.write(FRAGMENT_CONFIG_SCHEMA)
    bar = schemadir.join("bar.schema.yaml")
    bar.write(FRAGMENT_BAR_SCHEMA)
    return p


@pytest.fixture
def config_schema_relative_nested_with_default(schemadir):
    p = schemadir.join("config.relative.nested.default.schema.yaml")
    p.write(RELATIVE_CONFIG_SCHEMA)
    bar = schemadir.join("bar.schema.yaml")
    bar.write(NESTED_BAR_SCHEMA)
    baz = schemadir.join("baz.schema.yaml")
    baz.write(NESTED_BAZ_SCHEMA)
    return p


@pytest.fixture
def config_schema_allof_default(schemadir):
    p = schemadir.join("config.allof.default.schema.yaml")
    p.write(ALLOF_CONFIG_SCHEMA)
    return p


@pytest.fixture
def config_schema_defs(schemadir):
    p = schemadir.join("config.defs.schema.yaml")
    p.write(DEFS_CONFIG_SCHEMA)
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

    with pytest.raises(WorkflowError):
        validate(config, str(config_schema_ref), False)


def test_dataframe(df_schema):
    df = pd.DataFrame([{"sample": "foo", "condition": "bar"}])
    validate(df, str(df_schema), False)
    assert sorted(df.columns) == sorted(["sample", "condition"])
    validate(df, str(df_schema))
    assert sorted(df.columns) == sorted(["sample", "condition", "case", "date"])
    assert df.case.loc[0]


def test_config_ref_relative(config_schema_relative):
    config = {"bar": {}}
    with pytest.raises(WorkflowError):
        validate(config, str(config_schema_relative), set_default=False)
    config = {"bar": {"foo": "baz"}}
    validate(config, str(config_schema_relative), set_default=False)


def test_config_ref_relative_nested(config_schema_relative_nested):
    config = {"bar": {"baz": {}}}
    with pytest.raises(WorkflowError):
        validate(config, str(config_schema_relative_nested), set_default=False)
    config = {"bar": {"baz": {"qux": 19}}}
    validate(config, str(config_schema_relative_nested), set_default=False)


def test_config_ref_relative_with_remote_id(config_schema_relative_ref_with_remote_id):
    config = {"bar": {}}
    with pytest.raises(WorkflowError):
        validate(
            config, str(config_schema_relative_ref_with_remote_id), set_default=False
        )
    config = {"bar": {"foo": "baz"}}
    validate(config, str(config_schema_relative_ref_with_remote_id), set_default=False)


def test_config_ref_relative_with_fragment(config_schema_relative_ref_with_fragment):
    config = {"bar": None}
    with pytest.raises(WorkflowError):
        validate(
            config, str(config_schema_relative_ref_with_fragment), set_default=False
        )
    config = {"bar": {"baz": "value"}}
    validate(config, str(config_schema_relative_ref_with_fragment), set_default=False)


def test_config_allof_default(config_schema_allof_default):
    config = {"foo": {}}
    validate(config, str(config_schema_allof_default), set_default=True)
    assert config["foo"]["bar"] == 42
    assert config["foo"]["foo"] == "foo"


def test_config_anyof_via_defs_default(config_schema_defs):
    config = {"foo": {}}
    validate(config, str(config_schema_defs), set_default=True)
    assert config["foo"]["qux"] == "qux"
