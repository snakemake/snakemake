from snakemake.utils import validate
import pandas as pd

config_schema = "config.schema.yaml"
df_schema = "dataframe.schema.yaml"


def test_config():
    config = {}
    validate(config, config_schema)
    assert config == {}
    validate(config, config_schema, True)
    assert dict(config) == {'param': {'foo': 'bar'}}


def test_dataframe():
    df = pd.DataFrame([{'sample':'foo', 'condition':'bar'}])
    validate(df, df_schema)
    assert sorted(df.columns) == sorted(['sample', 'condition'])
    df = validate(df, df_schema, True)
    assert sorted(df.columns) == sorted(['sample', 'condition', 'case'])
    assert df.case.loc[0] == True
