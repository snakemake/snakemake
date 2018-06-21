import pytest
from snakemake.utils import validate
import pandas as pd

CONFIG_SCHEMA = """$schema: "http://json-schema.org/draft-07/schema#"
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

DF_SCHEMA = """$schema: "http://json-schema.org/draft-07/schema#"
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
def config_schema_ref(schemadir, bar_schema):
    p = schemadir.join("config.ref.schema.yaml")
    p.write(CONFIG_SCHEMA + "\n".join(
        ["      bar:",
         "        $ref: \"{bar}\"".format(bar=str(bar_schema))]))
    return p
    

def test_config(config_schema):
    config = {}
    validate(config, str(config_schema), False)
    assert config == {}
    validate(config, str(config_schema))
    assert dict(config) == {'param': {'foo': 'bar'}}


def test_config_ref(config_schema_ref):
    pass


def test_dataframe(df_schema):
    df = pd.DataFrame([{'sample':'foo', 'condition':'bar'}])
    validate(df, str(df_schema))
    assert sorted(df.columns) == sorted(['sample', 'condition'])
    df = validate(df, str(df_schema), True)
    assert sorted(df.columns) == sorted(['sample', 'condition', 'case'])
    assert df.case.loc[0] == True
