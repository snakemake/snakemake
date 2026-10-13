import inspect

from snakemake.common.configfile import _load_configfile
from snakemake.logging import logger
from snakemake.exceptions import WorkflowError


def validate(data, schema, set_default=True):
    """Validate data with JSON schema at given path.

    Args:
        data (object): data to validate. Can be a config dict or a pandas data frame.
        schema (str): Path to JSON schema used for validation. The schema can also be
            in YAML format. If validating a pandas data frame, the schema has to
            describe a row record (i.e., a dict with column names as keys pointing
            to row values). See https://json-schema.org. The path is interpreted
            relative to the Snakefile when this function is called.
        set_default (bool): set default values defined in schema. See
            https://python-jsonschema.readthedocs.io/en/latest/faq/ for more
            information
    """
    frame = inspect.currentframe().f_back
    workflow = frame.f_globals.get("workflow")

    if workflow and workflow.modifier.skip_validation:
        # skip if a corresponding modifier has been defined
        return

    from snakemake.sourcecache import LocalSourceFile, infer_source_file
    import jsonschema
    from jsonschema import Draft202012Validator, validators
    from referencing import Registry, Resource
    from pathlib import Path
    from urllib.parse import urlparse
    from urllib.request import url2pathname

    schemafile = infer_source_file(schema)

    if isinstance(schemafile, LocalSourceFile) and not schemafile.isabs() and workflow:
        # if workflow object is not available this has not been started from a workflow
        schemafile = workflow.current_basedir.join(
            schemafile.get_path_or_uri(secret_free=True)
        )

    source = (
        workflow.sourcecache.open(schemafile)
        if workflow
        else schemafile.get_path_or_uri(secret_free=False)
    )

    schema = _load_configfile(source, filetype="Schema")

    # Inject $id to enforce local reference resolution.
    # Ensures all relative $ref's are resolved relative to this local file,
    # even if it defines an explicit $id pointing to a remote location.
    # This is required because the retrieve_uri function (see below) does not
    # handle remote files.
    # This allows to mostly restores pre-9.6.0 behavior fixing the regression
    # caused by https://github.com/snakemake/snakemake/pull/3420 and reported
    # in https://github.com/snakemake/snakemake/issues/3648.
    # Note that the old (RefResolver based) implementation did handle remote
    # URI fetching and could therefore respect remote URIs defined as ID
    # within the config file schema.
    # To fully restore that behaviour, retrieve_uri would have to be expanded
    # to fetch remote resources and this injection would have to be condition
    # on the schema not defining an ID.
    # However, in the context of config file validation, resolving a reference
    # defined in a local schema file through a remote request seems to solve no
    # purpose: Either the schemas are identical, in which relying on the local,
    # known-to-exist version is more efficient, or they differ for some reason
    # in which case not using the ref text from the local schema file might
    # cause very surprising and hard to track down behaviour.
    # Therefore, the new (resolving based) implementation purposefully breaks
    # backwards-compatibility with the former (RefResolver based) one.
    # This is also made explicit through the
    # test_config_ref_relative_with_remote_id test added to test_schema.py.
    schema["$id"] = (
        Path(schemafile.get_path_or_uri(secret_free=False)).resolve().as_uri()
    )

    def retrieve_uri(uri):
        # Note:
        # Relative $ref's are resolved against the (referencing) schema's $id
        # by the referencing library before calling retrieve.
        # Above, this was set to the local file's (absolute) file:// URI.
        # Since _load_configfile expects a file handle/path, and not a URI,
        # it must be parsed to strip off the (URI) schema.
        return Resource.from_contents(
            contents=_load_configfile(
                url2pathname(urlparse(uri).path), filetype="Schema"
            )
        )

    resource = Resource.from_contents(contents=schema)
    # pyrefly: ignore [unexpected-keyword]
    registry = Registry(retrieve=retrieve_uri).with_resource(
        uri=schemafile.get_path_or_uri(secret_free=False), resource=resource
    )
    Validator = Draft202012Validator(schema, registry=registry)

    # Taken from https://python-jsonschema.readthedocs.io/en/latest/faq/
    def extend_with_default(validator_class):
        validate_properties = validator_class.VALIDATORS["properties"]

        def set_defaults(validator, properties, instance, schema):
            for property, subschema in properties.items():
                if "default" in subschema:
                    instance.setdefault(property, subschema["default"])

            for error in validate_properties(
                validator,
                properties,
                instance,
                schema,
            ):
                yield error

        return validators.extend(
            validator_class,
            {"properties": set_defaults},
        )

    if Validator.META_SCHEMA["$schema"] != schema["$schema"]:
        logger.warning(
            f"No validator found for JSON Schema version identifier '{schema['$schema']}'"
        )
        logger.warning(
            f"Defaulting to validator for JSON Schema version '{Validator.META_SCHEMA['$schema']}'"
        )
        logger.warning("Note that schema file may not be validated correctly.")
    Defaultvalidator = extend_with_default(Validator)

    def _validate_record(record):
        if set_default:
            Defaultvalidator(schema, registry=registry).validate(record)
            return record
        else:
            Validator.validate(record)

    def _validate_pandas(data):
        try:
            import pandas as pd

            if isinstance(data, pd.DataFrame):
                logger.debug("Validating pandas DataFrame")

                recordlist = []
                for i, record in enumerate(data.to_dict("records")):
                    # Exclude NULL values
                    record = {k: v for k, v in record.items() if pd.notnull(v)}
                    try:
                        recordlist.append(_validate_record(record))
                    except jsonschema.exceptions.ValidationError as e:
                        raise WorkflowError(
                            f"Error validating row {i} of data frame.", e
                        )

                if set_default:
                    newdata = pd.DataFrame(recordlist)
                    # Add missing columns and fill None values using positional alignment
                    # to avoid ValueError with duplicate indices (data.update() aligns on index)
                    for col in newdata.columns:
                        new_vals = newdata[col].to_numpy()
                        if col not in data.columns:
                            data[col] = new_vals
                        else:
                            cur_vals = data[col].to_numpy(copy=True)
                            na_mask = pd.isna(cur_vals)
                            if na_mask.any():
                                cur_vals[na_mask] = new_vals[na_mask]
                                data[col] = cur_vals

            else:
                return False
        except ImportError:
            return False
        return True

    def _validate_polars(data):
        try:
            import polars as pl

            if isinstance(data, pl.DataFrame):
                logger.debug("Validating polars DataFrame")

                recordlist = []
                for i, record in enumerate(data.iter_rows(named=True)):
                    # Exclude NULL values
                    record = {
                        k: v
                        for k, v in record.items()
                        if pl.Series(k, [v]).is_not_null().all()
                    }
                    try:
                        recordlist.append(_validate_record(record))
                    except jsonschema.exceptions.ValidationError as e:
                        raise WorkflowError(
                            f"Error validating row {i} of data frame.", e
                        )

                if set_default:
                    newdata = pl.DataFrame(recordlist)
                    # Add missing columns
                    newcol = [col for col in newdata.columns if col not in data.columns]
                    [
                        data.insert_column(
                            len(data.columns),
                            pl.lit(None, newdata[col].dtype).alias(col),
                        )
                        for col in newcol
                    ]
                    # Fill in None values with values from newdata
                    for i in range(data.shape[0]):
                        for j in range(data.shape[1]):
                            if data[i, j] is None:
                                data[i, j] = newdata[i, j]

            elif isinstance(data, pl.LazyFrame):
                # If a LazyFrame is being used, probably it is a large dataframe (so check only first 1000 records)
                logger.debug("Validating first 1000 rows of polars LazyFrame")

                recordlist = []
                for i, record in enumerate(
                    data.head(1000).collect().iter_rows(named=True)
                ):
                    # Exclude NULL values
                    record = {
                        k: v
                        for k, v in record.items()
                        if pl.Series(k, [v]).is_not_null().all()
                    }
                    try:
                        recordlist.append(_validate_record(record))
                    except jsonschema.exceptions.ValidationError as e:
                        raise WorkflowError(
                            f"Error validating row {i} of data frame.", e
                        )

                if set_default:
                    logger.warning("LazyFrame does not support setting default values.")

            else:
                return False
        except ImportError:
            return False
        return True

    if isinstance(data, dict):
        logger.debug("Validating dictionary")
        try:
            _validate_record(data)
        except jsonschema.exceptions.ValidationError as e:
            raise WorkflowError("Error validating config file.", e)
        logger.debug("Dictionary validated!")
    else:
        if _validate_pandas(data):
            logger.debug("Pandas dataframe validated!")
        elif _validate_polars(data):
            logger.debug("Polars dataframe validated!")
