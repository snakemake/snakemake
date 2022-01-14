# Changelog

### [6.13.1](https://www.github.com/snakemake/snakemake/compare/v6.13.0...v6.13.1) (2022-01-11)


### Bug Fixes

* --conda-frontend value not passed on to cluster jobs ([#1317](https://www.github.com/snakemake/snakemake/issues/1317)) ([df46ddb](https://www.github.com/snakemake/snakemake/commit/df46ddb37022b291a4feca22fd0fbcf8773e7d03))
* atomic job error display ([#1326](https://www.github.com/snakemake/snakemake/issues/1326)) ([aa2c265](https://www.github.com/snakemake/snakemake/commit/aa2c2652608d3e95ad7fb568df09ef1ae09e1def))
* fix source cache handling for remote source files retrieved via github() or gitlab() tags. ([#1322](https://www.github.com/snakemake/snakemake/issues/1322)) ([6e2ecd2](https://www.github.com/snakemake/snakemake/commit/6e2ecd26e48eb64fa04c9c38dde591857e03c722))
* typos in code examples ([#1324](https://www.github.com/snakemake/snakemake/issues/1324)) ([60010e4](https://www.github.com/snakemake/snakemake/commit/60010e4ef07b7ba9b89aa5f48ee90ff3cec85b75))

## [6.13.0](https://www.github.com/snakemake/snakemake/compare/v6.12.3...v6.13.0) (2021-12-21)


### Features

* allow prefix definition in module statements ([#1310](https://www.github.com/snakemake/snakemake/issues/1310)) ([29e6540](https://www.github.com/snakemake/snakemake/commit/29e6540aac95b08b5e386a8478bd2013334e5954))

### [6.12.3](https://www.github.com/snakemake/snakemake/compare/v6.12.2...v6.12.3) (2021-12-09)


### Bug Fixes

* fixed display of any exceptions and errors from within a workflow definition ([23d40d9](https://www.github.com/snakemake/snakemake/commit/23d40d99614a88fd3c596d05e6915509ae43d4ce))

### [6.12.2](https://www.github.com/snakemake/snakemake/compare/v6.12.1...v6.12.2) (2021-12-07)


### Bug Fixes

* rule inheritance within modules (did previously lead to key errors) ([#1292](https://www.github.com/snakemake/snakemake/issues/1292)) ([603e0a8](https://www.github.com/snakemake/snakemake/commit/603e0a87d2c7af57a8f1d397605bc501c50934e0))


### Documentation

* Fix typo in rules.rst (—draft-notebook) ([#1290](https://www.github.com/snakemake/snakemake/issues/1290)) ([f5c42cf](https://www.github.com/snakemake/snakemake/commit/f5c42cfdc68f1516cec71b8ead8d78225ae915e5))

### [6.12.1](https://www.github.com/snakemake/snakemake/compare/v6.12.0...v6.12.1) (2021-11-29)


### Bug Fixes

* set default number of nodes to 1 in test cases ([#1288](https://www.github.com/snakemake/snakemake/issues/1288)) ([f6e12b4](https://www.github.com/snakemake/snakemake/commit/f6e12b4798485be3a1bb240b4af44d57dd5c84b2))

## [6.12.0](https://www.github.com/snakemake/snakemake/compare/v6.11.1...v6.12.0) (2021-11-29)


### Features

* add flag --draft-notebook for generating a skeleton notebook for manual editing (e.g. in VSCode). ([#1284](https://www.github.com/snakemake/snakemake/issues/1284)) ([d279322](https://www.github.com/snakemake/snakemake/commit/d2793223f914790c07b25363cb9b314ef166cb3e))


### Bug Fixes

* issue [#1257](https://www.github.com/snakemake/snakemake/issues/1257) (missing logfile failure when using shadow directory) ([#1258](https://www.github.com/snakemake/snakemake/issues/1258)) ([426d92f](https://www.github.com/snakemake/snakemake/commit/426d92fd9610b61b414b7f0152d777c463c939a2))
* keep empty output and input dirs of --draft-notebook job ([f1181bd](https://www.github.com/snakemake/snakemake/commit/f1181bd41ea8b20fafd3975c2733ca1d439381dc))
* SameFileError [#1153](https://www.github.com/snakemake/snakemake/issues/1153) ([#1220](https://www.github.com/snakemake/snakemake/issues/1220)) ([ede313d](https://www.github.com/snakemake/snakemake/commit/ede313dcd31ea5f136b3b8f743e2265331475342))
* snakemake API using only 1 job as default ([#1283](https://www.github.com/snakemake/snakemake/issues/1283)) ([e92ad48](https://www.github.com/snakemake/snakemake/commit/e92ad4867feb456ce8ef3dc57fd8528affa64ae9))


### Documentation

* short tutorial updates ([#1286](https://www.github.com/snakemake/snakemake/issues/1286)) ([b653a44](https://www.github.com/snakemake/snakemake/commit/b653a44d105e4b3799425a695d75a08239dc0d6b))

### [6.11.1](https://www.github.com/snakemake/snakemake/compare/v6.11.0...v6.11.1) (2021-11-26)


### Bug Fixes

* provide temporary IPYTHONDIR for notebook execution in order to avoid race conditions in https://github.com/ipython/ipython/blob/master/IPython/paths.py#L20 upon execution of multiple notebooks at the same time. ([#1280](https://www.github.com/snakemake/snakemake/issues/1280)) ([4d70da1](https://www.github.com/snakemake/snakemake/commit/4d70da11f810224ddce192ae1472a6380898865f))


### Documentation

* move psutil import into benchmark methods to avoid needing it as a dependency for doc building ([6ffe38d](https://www.github.com/snakemake/snakemake/commit/6ffe38d1740294a7170765ab875b363f4ae82cd4))
* require sphinx>=3 ([1773875](https://www.github.com/snakemake/snakemake/commit/1773875fc8f2fddb09362410afb7c49c4406bfa3))
* skip lazy property ([2883718](https://www.github.com/snakemake/snakemake/commit/28837183fa55a6764621580983b3d724f3881a6a))

## [6.11.0](https://www.github.com/snakemake/snakemake/compare/v6.10.0...v6.11.0) (2021-11-25)


### Features

* fail with an error if snakemake cannot write job metadata. ([#1273](https://www.github.com/snakemake/snakemake/issues/1273)) ([cd968cd](https://www.github.com/snakemake/snakemake/commit/cd968cd03437ad6db1d791f5d7ae5295b9754137))


### Bug Fixes

* Adds fixes for the first two MREs in [#823](https://www.github.com/snakemake/snakemake/issues/823) ([#1215](https://www.github.com/snakemake/snakemake/issues/1215)) ([cfd2f89](https://www.github.com/snakemake/snakemake/commit/cfd2f890a0af57628f7b9278d8d43f59b7006825))
* env file usage after changes to source file handling (inspired by [#1233](https://www.github.com/snakemake/snakemake/issues/1233) and [#1211](https://www.github.com/snakemake/snakemake/issues/1211)). ([#1236](https://www.github.com/snakemake/snakemake/issues/1236)) ([3ac8e85](https://www.github.com/snakemake/snakemake/commit/3ac8e858a7b908326922c8f68cae512b1250e906))
* fixed code change detection when using modules ([#1264](https://www.github.com/snakemake/snakemake/issues/1264)) ([b571e09](https://www.github.com/snakemake/snakemake/commit/b571e09ce452f6a1a95395e1c3c8b9e3f83867ad))
* handle config file extension/overwriting more explicitly ([#1251](https://www.github.com/snakemake/snakemake/issues/1251)) ([d0a7bf2](https://www.github.com/snakemake/snakemake/commit/d0a7bf243c5df204136fa1f14706aab793793c68))
* Issue [#1253](https://www.github.com/snakemake/snakemake/issues/1253) (problems editing Jupyter Notebooks) ([#1255](https://www.github.com/snakemake/snakemake/issues/1255)) ([3398ddf](https://www.github.com/snakemake/snakemake/commit/3398ddffd1f68182af768ef4ea519e9a9ad4efaf))
* more informative nothing to be done message ([#1234](https://www.github.com/snakemake/snakemake/issues/1234)) ([368d265](https://www.github.com/snakemake/snakemake/commit/368d265ff3da984bd3a53b319dcb882d6916975b))
* only consider context of shell command for technical switches if called from snakemake rules. ([#1213](https://www.github.com/snakemake/snakemake/issues/1213)) ([4816a58](https://www.github.com/snakemake/snakemake/commit/4816a58653e466ca94b1482a1d947a856f5381b3))
* R encoding of pathlib.Path objects ([#1201](https://www.github.com/snakemake/snakemake/issues/1201)) ([bd516e9](https://www.github.com/snakemake/snakemake/commit/bd516e958af22e57c18cacf0cb22552c2a237bd8))
* Use 'snakemake.utils.update_config' instead of 'dict.update' ([#1126](https://www.github.com/snakemake/snakemake/issues/1126)) ([2658027](https://www.github.com/snakemake/snakemake/commit/2658027458dde4c10b3d6e1af7671564d175f9cb))

## [6.10.0](https://www.github.com/snakemake/snakemake/compare/v6.9.1...v6.10.0) (2021-10-21)


### Features

* Add more informative errors when evaluation of `--default-resources` fails ([#1192](https://www.github.com/snakemake/snakemake/issues/1192)) ([b3c4e68](https://www.github.com/snakemake/snakemake/commit/b3c4e687c87c75075393cef842b129dcec70e7f6))


### Bug Fixes

* add quotes to each item of the wait_for_files list ([#1160](https://www.github.com/snakemake/snakemake/issues/1160)) ([72856ed](https://www.github.com/snakemake/snakemake/commit/72856edd12fbe29d723731c6f596f05cd2b59c0e))
* caching process ([#1225](https://www.github.com/snakemake/snakemake/issues/1225)) ([0825a29](https://www.github.com/snakemake/snakemake/commit/0825a29e46c08b200efe6bd0c66acf1e6828eed8))
* enable usage of job grouping in GLS ([#1054](https://www.github.com/snakemake/snakemake/issues/1054)) ([d243c22](https://www.github.com/snakemake/snakemake/commit/d243c22ff494b63bd5e07b7c5bf1f6ff32539cde))
* Only --bind Snakemake when we're working with a Python script ([#1206](https://www.github.com/snakemake/snakemake/issues/1206)) ([1d79f62](https://www.github.com/snakemake/snakemake/commit/1d79f625b7262d66def71c779f2a2c091bc418d8))
* run dependencies with non-existent ancient files before the consuming job ([#1202](https://www.github.com/snakemake/snakemake/issues/1202)) ([84d1f64](https://www.github.com/snakemake/snakemake/commit/84d1f6451b12352eba5a8bfefcfcce8b2d98c5aa)), closes [#946](https://www.github.com/snakemake/snakemake/issues/946)
* status cmd repeats until killed by 11 *different* signals ([#1207](https://www.github.com/snakemake/snakemake/issues/1207)) ([8b28b57](https://www.github.com/snakemake/snakemake/commit/8b28b5740c34149c9b5df56dbbfa034219eb1574))
* typo in sourcecache use ([#1229](https://www.github.com/snakemake/snakemake/issues/1229)) ([8b54bc5](https://www.github.com/snakemake/snakemake/commit/8b54bc5db9d8e5c0bcb8f2c2ff141dc075e3e659))
* wms monitor arg parsing now accepts any kind of value ([#1181](https://www.github.com/snakemake/snakemake/issues/1181)) ([313de93](https://www.github.com/snakemake/snakemake/commit/313de932e2e2a4f2c530df18c1abb15d37eb3217))


### Documentation

* Clarification of --cluster-stats docs  &  elaborating on the situation where job ids are not passed to the status script ([#1221](https://www.github.com/snakemake/snakemake/issues/1221)) ([ed0e4a2](https://www.github.com/snakemake/snakemake/commit/ed0e4a27a2167a69a4fe1bcdf237dd27bb3732ca))
* Combine CHANGELOG.rst with CHANGELOG.md ([#1228](https://www.github.com/snakemake/snakemake/issues/1228)) ([19f5a43](https://www.github.com/snakemake/snakemake/commit/19f5a43261bd6ba548d6f01080640f0d4119871e))
* Mention required openssl dep for rust-script ([#1216](https://www.github.com/snakemake/snakemake/issues/1216)) ([fc8c5f6](https://www.github.com/snakemake/snakemake/commit/fc8c5f62c397a0239ef213ab45a26a1def50f9eb))
* Unpin docutils version ([#1230](https://www.github.com/snakemake/snakemake/issues/1230)) ([15a82bf](https://www.github.com/snakemake/snakemake/commit/15a82bfe402b3577bf19e6d2eca3b2fb86109628))

### [6.9.1](https://www.github.com/snakemake/snakemake/compare/v6.9.0...v6.9.1) (2021-09-30)


### Bug Fixes

* fix function call when creating report and hashes for between workflow caching ([#1198](https://www.github.com/snakemake/snakemake/issues/1198)) ([a4f6836](https://www.github.com/snakemake/snakemake/commit/a4f68365125c357f30510d0e61036f98b9d3aa69))

## [6.9.0](https://www.github.com/snakemake/snakemake/compare/v6.8.2...v6.9.0) (2021-09-29)


### Features

* autoconvert Path objects to str when passing to R or Julia scripts ([80ec513](https://www.github.com/snakemake/snakemake/commit/80ec51322f8134180c52c20b0a9dc6980df6c1bc))


### Bug Fixes

* fix source retrieval during between workflow caching and report generation ([2394ca4](https://www.github.com/snakemake/snakemake/commit/2394ca4a23a6b2792397bc9efc09945f01d1963b))

### [6.8.2](https://www.github.com/snakemake/snakemake/compare/v6.8.1...v6.8.2) (2021-09-29)


### Bug Fixes

* fix path returned by get_source() ([ee05315](https://www.github.com/snakemake/snakemake/commit/ee053153d2f44156171c127307cb110791b7624a))

### [6.8.1](https://www.github.com/snakemake/snakemake/compare/v6.8.0...v6.8.1) (2021-09-24)


### Bug Fixes

* async_run to allow nested event loops. ([#1170](https://www.github.com/snakemake/snakemake/issues/1170)) ([5dc6bbd](https://www.github.com/snakemake/snakemake/commit/5dc6bbd440ac46e81a926b6749969b98b7e33a9f))
* merging of pipe groups when multiple rules are chained together via pipes ([#1173](https://www.github.com/snakemake/snakemake/issues/1173)) ([de91d2c](https://www.github.com/snakemake/snakemake/commit/de91d2ccf53bd844b4dbf4f64dd087f4ee935be5))
* potential memory corruption caused by Google storage objects accessed from different threads ([#1174](https://www.github.com/snakemake/snakemake/issues/1174)) ([41a5071](https://www.github.com/snakemake/snakemake/commit/41a5071b750dca5d7fceec324d81d9a93c86bdb6))


### Performance Improvements

* more extensive caching of source files, including wrappers. ([#1182](https://www.github.com/snakemake/snakemake/issues/1182)) ([bdb75f8](https://www.github.com/snakemake/snakemake/commit/bdb75f828a3ae27ba97ea6cd5e71a34ac7b27eea))


### Documentation

* move note ([75a544b](https://www.github.com/snakemake/snakemake/commit/75a544ba528b30b43b861abc0ad464db4d6ae16f))
* polish ([47a7b62](https://www.github.com/snakemake/snakemake/commit/47a7b628686258a28dd870f20bf1f121b3a881c3))
* tutorial formatting ([594f5fb](https://www.github.com/snakemake/snakemake/commit/594f5fbb342e0722318641dea07d7da4c5eb8116))

## [6.8.0](https://www.github.com/snakemake/snakemake/compare/v6.7.0...v6.8.0) (2021-09-06)


### Features

* Add `shadow: "copy-minimal"` directive ([#1155](https://www.github.com/snakemake/snakemake/issues/1155)) ([1803f0b](https://www.github.com/snakemake/snakemake/commit/1803f0b4090d812df0c164653b26502fd130d326))
* support XRootD as a default remote provider ([#1017](https://www.github.com/snakemake/snakemake/issues/1017)) ([fe03157](https://www.github.com/snakemake/snakemake/commit/fe03157c31210984fce53c35d5fb87b20d278fe7))


### Bug Fixes

* AmbiguousRuleException bug caused by weak ordering of rules ([#1124](https://www.github.com/snakemake/snakemake/issues/1124)) ([7f54c39](https://www.github.com/snakemake/snakemake/commit/7f54c391f2821655ed168bcdafad6d07b96fcec7))
* Bugfix tes add files ([#1133](https://www.github.com/snakemake/snakemake/issues/1133)) ([8892bf2](https://www.github.com/snakemake/snakemake/commit/8892bf25d9d981a4032d5a1b525960ba3bdd1aec))
* Disable Persistence cache for snakemake jobs ([#1159](https://www.github.com/snakemake/snakemake/issues/1159)) ([7110f9d](https://www.github.com/snakemake/snakemake/commit/7110f9d2e7ee3f350bd1da3c5b4aab98c06725a1))
* efficient job status checking when using DRMAA API (this should yield much better parallelization and performance when using --drmaa) ([#1156](https://www.github.com/snakemake/snakemake/issues/1156)) ([ac004cb](https://www.github.com/snakemake/snakemake/commit/ac004cb19cebd4efb5e38f6039861a2810c702ff))
* improved error handling for cluster status scripts and smarter job selector choice in case of cluster submission (use greedy for single jobs). ([#1142](https://www.github.com/snakemake/snakemake/issues/1142)) ([48d2dd9](https://www.github.com/snakemake/snakemake/commit/48d2dd99a745fd54b74b1435cbb7e41e0ee1b4ac))
* Initialize assignments dictionary when setting rule-based resources ([#1154](https://www.github.com/snakemake/snakemake/issues/1154)) ([68c13fd](https://www.github.com/snakemake/snakemake/commit/68c13fd6fb2ad458e79bafe146499b601bf4bd0e))
* key error when handling FileNotFoundError in input functions. ([#1138](https://www.github.com/snakemake/snakemake/issues/1138)) ([d25f04d](https://www.github.com/snakemake/snakemake/commit/d25f04db820c9651835b7323baef5931d4f8dc0a))
* linting of remote snakefiles ([#1131](https://www.github.com/snakemake/snakemake/issues/1131)) ([2104e10](https://www.github.com/snakemake/snakemake/commit/2104e10d1d2c5e0f368e9c0fe95cc50f9d4847f1))


### Performance Improvements

* improve job selection performance in case of potential ambiguity that is resolved by comprehensive ruleorder statements. ([#1147](https://www.github.com/snakemake/snakemake/issues/1147)) ([921f4f7](https://www.github.com/snakemake/snakemake/commit/921f4f715e3814fc2e22a4f6527ff62e066cc5da))

## [6.7.0](https://www.github.com/snakemake/snakemake/compare/v6.6.1...v6.7.0) (2021-08-12)


### Features

* Add support for rust scripts (enabling directly integrated ad-hoc robust high performance scripting) ([#1053](https://www.github.com/snakemake/snakemake/issues/1053)) ([f0e8fa2](https://www.github.com/snakemake/snakemake/commit/f0e8fa285437a02ca7edcf87334bf00cb347064a))


### Bug Fixes

* Ga4gh tes bugfixes ([#1127](https://www.github.com/snakemake/snakemake/issues/1127)) ([af21d6c](https://www.github.com/snakemake/snakemake/commit/af21d6c2b125c22ef3dbc36a0a6a67a1874549c7))
* improved display of percentage of done jobs ([1fee8c0](https://www.github.com/snakemake/snakemake/commit/1fee8c06d6ed229d7e3757de3c693e755d01d1bb))
* improved error message in case of target rule misspecification ([83b1f5b](https://www.github.com/snakemake/snakemake/commit/83b1f5bbde437e13641be2160f4855f54043c046))


### Documentation

* fix contributing executors link ([#1112](https://www.github.com/snakemake/snakemake/issues/1112)) ([4bb58d1](https://www.github.com/snakemake/snakemake/commit/4bb58d12a44f77f79d47c5443f927cb6061677f5))
* Fix typo in file path in remote files documentation ([#1110](https://www.github.com/snakemake/snakemake/issues/1110)) ([9ce294f](https://www.github.com/snakemake/snakemake/commit/9ce294f6d5bdf72055a824ab610488a7f832a4d3))

### [6.6.1](https://www.github.com/snakemake/snakemake/compare/v6.6.0...v6.6.1) (2021-07-19)


### Bug Fixes

* avoid superfluous calls of conda info that have slowed down Snakemake since 6.4.1. ([#1099](https://www.github.com/snakemake/snakemake/issues/1099)) ([e990927](https://www.github.com/snakemake/snakemake/commit/e9909273c22a316dbd7301a243498e3c2a372642))

## [6.6.0](https://www.github.com/snakemake/snakemake/compare/v6.5.5...v6.6.0) (2021-07-16)


### Features

* Allow to mark all output files as temp with --all-temp ([#1097](https://www.github.com/snakemake/snakemake/issues/1097)) ([0ac3b38](https://www.github.com/snakemake/snakemake/commit/0ac3b3806c065d0ec3a551a5992faf30ddcf0576))

### [6.5.5](https://www.github.com/snakemake/snakemake/compare/v6.5.4...v6.5.5) (2021-07-16)


### Bug Fixes

* dummy release ([e4dca50](https://www.github.com/snakemake/snakemake/commit/e4dca508f6cbd3427d8580ef61f274f909ec8bab))

### [6.5.4](https://www.github.com/snakemake/snakemake/compare/v6.5.3...v6.5.4) (2021-07-16)


### Fixes

* Fixed --touch in combination with temp files (issue #1028) (@johanneskoester, @iromeo).

### Documentation

* Fix syntax error in docs/conf.py and update sphinx.ext.napoleon import ([#1084](https://www.github.com/snakemake/snakemake/issues/1084)) ([3e3fac2](https://www.github.com/snakemake/snakemake/commit/3e3fac2dbd5a8abad67e252f6181ad14bcfcb711))
* Improved pepfile (pepschema) documentation (@stolarczyk).
### \[6.5.3\] - 2021-07-06

-   Fixed a bug occuring when using --resources in the command line
    interface (@johanneskoester).
-   Minor improvements in the docs (@johanneskoester).

### \[6.5.2\] - 2021-07-02

-   Create directory pointed to by tmpdir resource if it does not yet
    exist (@johanneskoester).
-   Use a single core again in dryrun if --cores is not specified
    (@johanneskoester).
-   Bugfix for FTP remote provider (@jmeppley).
-   Improved documentation (@corneliusroemer).

### \[6.5.1\] - 2021-06-24

-   Extended best practices document (@johanneskoester)
-   Restore `-j all` behavior for local execution as a (deprecated) way
    of running Snakemake on all cores. Recommended now: `--cores all`
    (@johanneskoester).
-   Improved handling and better error messages for checkpoints
    (@johanneskoester).

## \[6.5.0\] - 2021-06-22

-   Allow to set the default profile via the environment variable
    $SNAKEMAKE_PROFILE.
-   There is a new default resource tmpdir (by default reflects the
    system setting), which is automatically used for temporary files by
    shell commands and scripts which properly consider the usual
    environment variables like $TMP, $TEMP, $TMPDIR (@johanneskoester).
-   The CLI flags --jobs and --cores are now separated, with --cores
    being responsible for local cores and global cores in the cluster
    case, and --jobs being responsible for number of jobs. Still -j and
    --jobs works as a fallback for local execution (@johanneskoester).
-   Added the ability to overwrite resources via --set-resources
    (@johanneskoester).
-   Various fixes for Windows execution (@melund).
-   Fixed a bug with fractional resources (@johanneskoester).
-   Fixed timeouts and other issues in google life science backend
    (@johanneskoester).
-   Fixed a bug with missing conda frontend definitions in subworkflows
    (@johanneskoester).
-   Skip envvar checking during linting (@johanneskoester).
-   Fixed a bug causing container images in modules to be ignored
    (@johanneskoester).

### \[6.4.1\] - 2021-05-27

-   Fixed bug in `workflow.source_path()` that occurred with modules
    included from remote locations (@johanneskoester).
-   Inform cluster jobs about conda/mamba/activate path such that they
    don't need to determine this themselves (@johanneskoester).

## \[6.4.0\] - 2021-05-20

-   Improvements in the docs (resource usage, best practices, remote
    files) (@johanneskoester, @admorris).
-   functions given to params can now safely open input files generated
    by previous rules. If they are not present, TBD will be displayed
    and function will be reevaluated immediately before the job is
    executed (i.e. when files are present) (@ASLeonard).
-   Connection pool for SFTP and FTP remote files, increasing download
    performance (@jmeppley).
-   Require correct minimum version of smart_open
    (@Redmar-van-den-Berg).
-   Added workflow.source_path(path), allowing to get the correct path
    relative to the current Snakefile, even when Snakefile is included
    via URL (@johanneskoester).
-   Fixed bugs in module system (@johanneskoester, @dlaehnemann).
-   Fixed issue with checkpoints and ruleorder where phantom
    dependencies are not properly removed from the DAG (@jmeppley,
    @johanneskoester).
-   Disable tibanna behavior that opens a browser window for each job
    (@nigiord).
-   Allow `Paramspace(..., filename_params="*")`, meaning that all
    columns of the paramspace will be encoded into the filename (@kpj).
-   Avoid PATH modification in cluster jobs (@johanneskoester).
-   For large sets of input files, pass files to wait for (FS latency)
    as a file instead of command line args (@kpj, @epruesse).

## \[6.3.0\] - 2021-04-29

-   Changed behavior of `workflow.snakefile` to always point to the
    current file instead of the main Snakefile (also in case of includes
    and modules) (@johanneskoester).
-   Fixed a typo in an error message (@nikostr).

## \[6.2.0\] - 2021-04-22

-   Support for integration of foreign workflow management systems by
    introducing a `handover` directive that passes on all resources to a
    particular rule (which can then invoke another workflow management
    system). See the docs ("Integrating foreign workflow management
    systems") (@johanneskoester).
-   Behavior improvement for temp handling of checkpoint rules
    (@epruesse).
-   Several improvements in the docs (@johanneskoester).

### \[6.2.1\] - 2021-04-20

-   Fixed a minor bug in the linter.

## \[6.2.0\] - 2021-04-20

-   Fixed several glitches in paramspace implementation (handling of
    bools, returning scalar values) (@kpj).
-   Fixed bugs in module implementation (@dlaehnemann,
    @johanneskoester).
-   Fall back to greedy scheduling solver if ILP solver needs more than
    10 sec (@johanneskoester).

### \[6.1.1\] - 2021-04-07

-   Fixed several small bugs of the new module system (@johanneskoester,
    @dlaehnemann).
-   Fixed archive based conda deployment (@johanneskoester).
-   Better handling of download and target attributed in the interactive
    report (@johanneskoester).

## \[6.1.0\] - 2021-04-01

-   Snakemake now uses **mamba** as the default conda frontend (which
    can be overwritten by specifying to use conda via the
    --conda-frontend flag) (@johanneskoester).
-   Profiles using --cluster option can now handle relative submit
    script paths in combination with arguments (@kdm9).
-   New AutoRemoteProvider, which infers the type of remote file
    protocol from the given URL (@kpj).
-   When using global container directive, container usage can be
    deactivated on a per rule base (@bilke).
-   Bugfixes for checkpoint handling (@johanneskoester).
-   Bugfixes for the module system (@johanneskoester, @dlaehnemann).
-   Various improvements for the tutorial.

### \[6.0.5\] - 2021-03-11

-   Fix bug (introduced with 6.0) when handling of HTML directories in
    report (@johanneskoester).

### \[6.0.4\] - 2021-03-11

-   Various textual improvements in the tutorial (@dlaehnemann).

### \[6.0.3\] - 2021-03-08

-   No longer use a shortened hash for naming conda environments in
    .snakemake/conda (@johanneskoester).
-   Various little updates to the docs (@johanneskoester).

### \[6.0.2\] - 2021-03-03

-   Fix race condition in conda checking code (@johanneskoester).

### \[6.0.1\] - 2021-03-03

-   Restored Python 3.5 compatibility by removing f-strings (@mbhall88)
-   Fix rendering issue in the docs.
-   Add gitpod dev environment and gitpod environment for the tutorial.

## \[6.0.0\] - 2021-02-26

-   Introduced a new module system, see
    <https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules>
    (@johanneskoester).
-   Introduced a rule inheritance mechanism, see
    <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rule-inheritance>
    (@johanneskoester).
-   Automatically containerize a conda-based pipeline with
    `--containerize`, see
    <https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#containerization-of-conda-based-workflows>
    (@johanneskoester).
-   Use temporary files for long shell commands (@epruesse).
-   Various fixes in the documentation (@ctb, @SilasK, @EthanHolleman).
-   Fixed a bug in job grouping that led to non-deterministic behavior
    (@johanneskoester).

### \[5.32.2\] - 2021-02-11

### Changed

-   Fixed infinite loading of results in Snakemake
    reports (@FelixMoelder)

### \[5.32.1\] - 2021-02-08

### Changed

-   Improved warning on wildcard constraints (@jheuel)
-   Improved logging from the new scheduler
    implementation (@johanneskoester)
-   Restored Python 3.5 compatibility by removing f-strings (@mbhall88)
-   Snakemake now automatically adds a global wildcard constraint for
    {scatteritem}, when scatter/gather support is used.
-   The zip variant of Snakemake reports is now compressed
    (@FelixMoelder).
-   Improved docs (@ctb).
-   Make output file removal in cluster mode more robust (@sebschmi).

## \[5.32.0\] - 2021-01-15

#### Changed

-   Handle accidental use of GLS backend with singularity (@vsoch).
-   Improved and extended WMS-monitor implementation (@vsoch).
-   Display index and total count in `{scatteritem}` when using the
    scatter-gather helper (@johanneskoester).
-   Fixed problems with jobid display when handling checkpoint updates
    (@johanneskoester, @jmeppley).
-   Fixed bug when checking for directory containment of output files
    (@jmeppley).
-   Implement --no-subworkflows treatment in combination with --cluster
    (@goi42).

### \[5.31.1\] - 2020-12-21

#### Changed

-   added wget again to the container image

## \[5.31.0\] - 2020-12-21

#### Added

-   The `Paramspace` helper for automatically exploring parameter spaces
    given as Pandas dataframes.
-   A new directive `name:` for setting rule names from variables.

#### Changed

-   Various small bug fixes for scheduling and checkpoint handling.
-   Automatically block R_LIBS, PYTHONPATH, PERL5LIB, and PERLLIB when
    using conda with --use-conda. This behavior can be deactivated with
    --conda-not-block-envvars.
-   Update container image to latest singularity.

### \[5.30.2\] - 2020-12-16

#### Changed

-   Fix permission issues with jobscripts on some systems (@Phhere).
-   Added notes on WSL to the tutorial (@RomainFeron).
-   Scheduler fixes (@johanneskoester).
-   Fixed a bug in checkpoint handling that led to hanging workflow
    execution (@jmeppley).
-   Pass cluster nodes to subworkflows (@votti).
-   Fix start time recording in metadata (@lparsons).
-   Fix time retrieval in reports (@johanneskoester).
-   Fix error when returning a Path from an input function (@sappjw).
-   Extending monitoring docs with some notes about future api changes
    (@vsoch).

## \[5.30.0\] - 2020-11-23

#### Added

-   Benchmarks now also report CPU time (@natir).

#### Changed

-   Fixed a reauthentication bug in Kubernetes support (@haizi-zh).

## \[5.29.0\] - 2020-11-19

#### Changed

-   Fixed several bugs in reports and scheduler.
-   Remove automatic (but buggy) encoding of csv/tsv files into HTML
    tables in the report (we will soon have a better alternative).
-   Fixed bug in kubernetes executor occurring with large source files.

## \[5.28.0\] - 2020-11-12

#### Added

-   Execution backend for GA4GH TES (task execution scheduler) an
    abstraction layer for various cluster and cloud queuing systems
    (@svedziok, @uniqueg).
-   script, notebook, wrapper and cwl directives now permit to use
    wildcards and params for composing paths (@johanneskoester).

#### Changed

-   Restored compatibility with Python 3.5 and 3.6 (@cclienti).
-   Various usability bug fixes (@goi43, @johanneskoester, @dcroote).
-   Better and more secure parsing of values when using --config
    (@bingxiao).

### \[5.27.4\] - 2020-11-03

#### Changed

-   Further speed improvements for DAG computation.
-   Fixed metadata migration errors occuring with long output file
    paths.
-   Add WorkflowHub specifications to the docs.
-   Fix group assignments.

### \[5.27.3\] - 2020-10-30

#### Changed

-   Added missing files to source distribution.

### \[5.27.2\] - 2020-10-30

#### Changed

-   DAG computation runtime has been improved by orders of magnitude, it
    is linear in the number of jobs now (@mhulsmann, @johanneskoester).
-   Stat calls have been dramatically reduced and are now performed in
    parallel (@johanneskoester).
-   Scheduler fixes (@FelixMoelder).
-   Directory support and other fixes for Google Life Sciences backend
    (@vsoch, @millerdz).
-   Support for panoptes monitor server (@fgypas).
-   Extended pathlib support (@mbhall88).
-   Vim plugin improvements (@troycomi).
-   Prevent jobs being rerun when input files are marked as ancient and
    another job in the DAG creates them.
-   Fixed --list-code-changes for included rules (@jbloom).

#### Added

-   Syntax highlighting for nano (@baileythegreen).

### \[5.26.1\] - 2020-10-01

#### Changed

-   Use coin ILP solver for scheduling by default (GLPK has bugs that
    can cause it to fail in certain situations).
-   If coin is not available, fall back to greedy scheduler.

## \[5.26.0\] - 2020-09-30

#### Added

-   Flag --max-inventory-time for setting maximum time spend on creating
    file inventory.
-   Flag --scheduler-ilp-solver for defining which solver to use for the
    ILP scheduler.

#### Changed

-   Fixed various bugs with the new scheduler (@FelixMoelder).
-   Fixed bug causing certain parameters not to be passed to the cluster
    (--set-scatter, --scheduler, --set-threads).
-   Updated docs and fixed of google backend (@vsoch).
-   Display jupyter notebook code in reports.
-   Improved scheduler behavior in order to directly remove temporary
    files if possible.

## \[5.25.0\] - 2020-09-18

#### Added

-   Simplified and more configurable support for scatter-gather
    processes (see docs).
-   Fully configurable DAG partitioning by grouping jobs at the command
    line. This should provide a vast additional improvement to
    scalability in cluster and cloud settings.

#### Changed

-   Depend on latest pulp, thereby enable Python >=3.8 compatibility
    again.
-   Fixes for snakefile handling in google life sciences backend
    (@vsoch).

### \[5.24.2\] - 2020-09-15

#### Changed

-   Fixed a bug in the linter that caused a false warning when using
    resources in shell commands.

### \[5.24.1\] - 2020-09-13

#### Changed

-   Depend on pulp \< 2.0, which includes the default coin cbc solver
    for all platforms.

## \[5.24.0\] - 2020-09-09

#### Added

-   Preemtion support for google cloud backend (@vsoch).

#### Changed

-   Fixed compatibility issues in new scheduler code (@dtrodrigues and
    @johanneskoester).
-   Improved error messages (@Sam-Tygier, @terrycojones)
-   Various small bug fixes.
-   Improved profile documentation (@johanneskoester).

## \[5.23.0\] - 2020-08-24

#### Added

-   Support for workflow configuration via portable encapsulated
    projects (PEPs, <https://pep.databio.org>).
-   A new ILP based default scheduler now ensures that temporary files
    are deleted as fast as possible (@FelixMoelder, @johanneskoester).

#### Changed

-   Fixed bug in modification date comparison for files in google
    storage (@vsoch).
-   Various small documentation improvements (@dcroote, @erjel,
    @dlaehnemann, @goi42).

### \[5.22.1\] - 2020-08-14

#### Changed

-   Fixed a missing dependency for google storage in cloud execution.

## \[5.22.0\] - 2020-08-13

#### Added

-   Added short option `-T` for CLI parameter `--restart-times`
    (@mbhall88).

#### Changed

-   Various small fixes for google storage and life sciences backends
    (@vsoch).

## \[5.21.0\] - 2020-08-11

#### Changed

-   Added default-remote-provider support for Azure storage
    (@andreas-wilm).
-   Various small bug fixes and documentation improvements.

### \[5.20.1\] - 2020-07-08

#### Changed

-   Fixed a bug that caused singularity args to be not passed on
    correctly when using script or conda.

## \[5.20.0\] - 2020-07-08

#### Changed

-   Exceptions in input functions are now handled in a smarter way, by
    choosing alternative paths in the DAG if available.
-   Debugging dag creation (--debug-dag) now gives more hints if
    alternative DAG paths are chosen.
-   Fixes for XRootD remote file implementation.
-   Improved CLI documentation.
-   Improved docs.
-   Various minor bug fixes.
-   Restored Python 3.5 compatibility.
-   Speed improvements for workdir cleanup.
-   Allow Path objects to be passed to expand.

### \[5.19.3\] - 2020-06-16

#### Changed

-   Performance improvements for DAG generation (up to 7x in the google
    cloud, anything from a little to massive in a cluster, depending on
    the overall filesystem performance).
-   Made harcoded bucket in google cloud executor configurable.
-   Improved speed of --unlock command.

### \[5.19.2\] - 2020-06-04

#### Changed

-   Fixed a bug in script and wrapper directives. Tried to decode a str.

### \[5.19.1\] - 2020-06-03

#### Changed

-   Fixed an issue with the parameter linting code, that could cause an
    index out of bounds exception.

## \[5.19.0\] - 2020-06-02

#### Added

-   The multiext function now allows arbitrary file extensions (no
    longer required to start with a "." (thanks to @jafors)
-   The include directive can now also take a Pathlib Path object
    (thanks to @mbhall88).

#### Changed

-   Jupyter notebook integration no longer automatically starts a
    browser.
-   Empty directories are cleaned up after workflow execution.
-   Fixed directory handling: no longer fail if the same job writes both
    a dir and a contained file.
-   Linter now recommends using spaces only for indentation.
-   Persistence dir "aux" has been renamed to "auxilliary" in order to
    make windows happy.
-   Linter now distinguishes awk syntax from regular variable usage.
-   Various bug fixes for Windows (thanks to @melund).

## \[5.18.0\] - 2020-05-21

#### Added

-   Native Google Cloud support via the (despite the name generic)
    lifesciences API.
-   Ability to optionally exchange the conda frontend to mamba (faster
    and sometimes more correct) instead of conda.

#### Changed

- Improved notebook integration experience, with various removed bugs and
  pitfalls.
- Auto-retry google storage API calls on transient or checksum errors.

## \[5.17.0\] - 2020-05-07

#### Added

- --envvars flag for passing secrets to cloud executors

#### Changed

- Wider thumbnail dialogs in report.
- Updated installation instructions.
- Various small kubernetes bug fixes.
- Bug fix for iRods remote files.

## \[5.16.0\] - 2020-04-29

#### Added

- Interactive jupyter notebook editing. Notebooks defined by rules can
  be interactively drafted and updated using snakemake --edit-notebook
  (see docs).

#### Changed

- Fixed group resource usage to occupy one cluster/cloud node.
- Minor bug fixes.

## \[5.15.0\] - 2020-04-21

#### Changed

-   The resource directive can now take strings, e.g. for defining a GPU
    model (see docs). This will e.g. be used for upcoming updates to
    cloud executors.
-   More extensive conda cleanup with --conda-cleanup-packages, meant
    for CI usage.
-   Further polish for reports.

## \[5.14.0\] - 2020-04-08

#### Changed

-   Redesigned HTML reports, with improved interface and performance.
-   For big data, HTML reports can now be stored as ZIP, where files are
    not anymore embedded but rather are stored in an auxilliary folder,
    such that they don't have to be in memory during report rendering.
-   Added subcategories to report (see docs).
-   Fixed a bug linter, leading to only one rule or snakefile to be
    linted.
-   Breaking change in CLI: added flags --conda-cleanup-envs and
    --conda-cleanup-pkgs, removed flag --cleanup-conda.
-   Fixed scheduling of pipe jobs, they are now always scheduled, fixing
    a hangup.
-   Corrected quoting of shell command for cluster submission.

## \[5.13.0\] - 2020-03-27

#### Added

- Allow to flag directories for inclusion in the report.

#### Changed

- Fixed hash computation for --cache in case of positional params
  arguments.
- Automatically restrict thread usage of linear algebra libraries to whatever
  is specified in the rule/job.

### \[5.12.3\] - 2020-03-24

#### Changed

-   Various minor bug fixes.

### \[5.12.2\] - 2020-03-24

#### Changed

-   Further improved linter output.

### \[5.12.1\] - 2020-03-24

#### Changed

-   Linter fixes

## \[5.12.0\] - 2020-03-24

#### Changed

-   Fixed the ability to supply functions for the thread directive.
-   Improved error messages for caching.

#### Added

-   A new "cache: true" directive that allows to annotate between
    workflow caching eligibility for rules in the workflow.

### \[5.11.2\] - 2020-03-19

#### Changed

-   Fixed a spurious error message complaining about missing singularity
    image if --use-singularity is not activated.

### \[5.11.1\] - 2020-03-16

#### Changed

-   Fixed a KeyError bug when executing a workflow that defines
    containers without --use-singularity.

## \[5.11.0\] - 2020-03-16

#### Changed

-   Fixes for environment modules and tibanna-based AWS execution.
-   Fixes for --default-resources defaults.
-   --cores is now a mandatory argument!
-   Automatic checksum validation for google storage.

#### Added

-   Azure storage authentication via SAS
-   A generic container directive that will in the future allow for
    other backends than just singularity. This deprecates the
    singularity directive, which will however stay functional at least
    until the next major release.
-   envvars directive for asserting environment variable existence. See
    docs.
-   support for AWS spot instances via --tibanna-config spot=true.
-   Automatic code quality linting via --lint.

## \[5.10.0\] - 2020-01-20

#### Added

-   Jupyter notebook integration, see docs. This enables interactive
    development of certain data analysis parts (e.g. for plotting).
-   Ability to overwrite thread definitions at the command line
    (`--threads rulename=3`), thereby improving scalability.
-   Requester pays configuration for google storage remote files.
-   Add keyword `allow_missing` to expand function, thereby allowing
    partical expansion by skipping wildcards for which no keywords are
    defined.

#### Changed

-   Various bug fixes, e.g. for between workflow caching and script
    execution.

### \[5.9.1\] - 2019-12-20

#### Changed

-   Added a missing module.

## \[5.9.0\] - 2019-12-20

#### Added

-   Support for per-rule environment module definitions to enable HPC
    specific software deployment (see docs).
-   Allow custom log handler defitions via --log-handler-script (e.g.
    post errors and progress to a slack channel or send emails).

-   Allow setting threads as a function of the given cores (see docs).

#### Changed

- Various minor fixes.

### \[5.8.2\] - 2019-12-16

#### Added

- Implemented a `multiext` helper, allowing to define a set of output
files that just differ by extension.

#### Changed

- Fixed a failure when caching jobs with conda environments.
- Fixed various minor bugs.
- Caching now allows to cache the output of rules using `multiext`.

### \[5.8.1\] - 2019-11-15

#### Changed

-   Fixed a bug by adding a missing module.

## \[5.8.0\] - 2019-11-15

#### Added

-   Blockchain based caching between workflows (in collaboration with
    Sven Nahnsen from QBiC), see [the
    docs](https://snakemake.readthedocs.io/en/v5.8.0/executing/caching.html).
-   New flag --skip-cleanup-scripts, that leads to temporary scripts
    (coming from script or wrapper directive) are not deleted (by Vanessa
    Sochat).

#### Changed

- Various bug fixes.

### \[5.7.4\] - 2019-10-23

#### Changed

-   Various fixes and adaptations in the docker container image and the
    test suite.

### \[5.7.1\] - 2019-10-16

#### Added

- Ability to print log files of failed jobs with --show-failed-logs.

#### Changed

- Fixed bugs in tibanna executor.
- Fixed handling of symbolic links.
- Fixed typos in help texts.
- Fixed handling of default resources.
- Fixed bugs in azure storage backend.

## \[5.7.0\] - 2019-10-07

#### Changed

-   Fixed various corner case bugs. Many thanks to the community for
    pull requests and reporting!
-   Container execution adapted to latest singularity.

#### Added

-   First class support for Amazon cloud execution via a new
    [Tibanna backend](https://snakemake.readthedocs.io/en/v5.7.0/executable.html#executing-a-snakemake-workflow-via-tibanna-on-amazon-web-services).
    Thanks to Soo Lee from Harvard Biomedical Informatics!
-   Allow multiple config files to be passed via the command line.
-   A new, more detailed way to visualize the DAG (--filegraph). Thanks
    to Henning Timm!
-   Pathlib compatibility added. Input and output files can now also be
    Path objects. Thanks to Frederik Boulund!
-   New azure storage remote provider. Transparently access input and
    output files on Microsoft Azure. Thanks to Sebastian Kurscheid!

## \[5.6.0\] - 2019-09-06

#### Changed

-   Fix compatibility with latest singularity versions.
-   Various bug fixes (e.g. in cluster error handling, remote providers,
    kubernetes backend).

#### Added

- Add --default-resources flag, that
  allows to define default resources for jobs (e.g. mem_mb, disk_mb), see
  [docs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources).
- Accept `--dry-run` as a synonym of `--dryrun`. Other Snakemake options
  are similarly hyphenated, so other documentation now refers to
  `--dry-run` but both (and also `-n`) will always be accepted
  equivalently.

### \[5.5.4\] - 2019-07-21

#### Changed

-   Reports now automatically include workflow code and configuration
    for improved transparency.

### \[5.5.3\] - 2019-07-11

#### Changed

-   Various bug fixes.
-   Polished reports.

### \[5.5.2\] - 2019-06-25

#### Changed

-   Various minor bug fixes in reports.
-   Speed improvements when using checkpoints.

### \[5.5.1\] - 2019-06-18

#### Changed

-   Improved report interface. In particular for large files.
-   Small TSV tables are automatically rendered as HTML with datatables.
-   Be more permissive with Snakefile choices: allow "Snakefile",
    "snakefile", "workflow/Snakefile", "workflow/snakefile".

## \[5.5.0\] - 2019-05-31

#### Added

- Script directives now also support Julia.

#### Changed

- Various small bug fixes.

### \[5.4.5\] - 2019-04-12

#### Changed

-   Fixed a bug with pipe output.
-   Cleaned up error output.

### \[5.4.4\] - 2019-03-22

#### Changed

-   Vastly improved performance of HTML reports generated with --report,
    via a more efficient encoding of dara-uri based download links.
-   Tighter layout, plus thumbnails and a lightbox for graphical results
    in HTML reports.
-   Bug fix for pipe groups.
-   Updated docs.
-   Better error handling in DRMAA executor.

### \[5.4.3\] - 2019-03-11

#### Changed

-   More robust handling of conda environment activation that should
    work with all setups where the conda is available when starting
    snakemake.
-   Fixed bugs on windows.

### \[5.4.2\] - 2019-02-15

#### Changed

-   Fixed a bug where git module cannot be imported from wrapper.

### \[5.4.1\] - 2019-02-14

#### Added

-   Warning when R script is used in combination with conda and R_LIBS
    environment variable is set. This can cause unexpected results and
    should be avoided.

#### Changed

-   Improved quoting of paths in conda commands.
-   Fixed various issues with checkpoints.
-   Improved error messages when combining groups with cluster config.
-   Fixed bugs in group implementation.
-   Fixed singularity in combination with shadow.

## \[5.4.0\] - 2018-12-18

#### Added

-   Snakemake now allows for data-dependent conditional re-evaluation of
    the job DAG via checkpoints. This feature also deprecates the
    `dynamic` flag. See [the
    docs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution).

### \[5.3.1\] - 2018-12-06

#### Changed

-   Various fixed bugs and papercuts, e.g., in group handling,
    kubernetes execution, singularity support, wrapper and script usage,
    benchmarking, schema validation.

## \[5.3.0\] - 2018-09-18

#### Added

-   Snakemake workflows can now be exported to CWL via the flag
    --export-cwl, see [the
    docs](https://snakemake.readthedocs.io/en/stable/executing/interoperability.html).

#### Changed

-   Fixed bug in script and wrapper execution when using
    `--use-singularity --use-conda`.
-   Add host argument to S3RemoteProvider.
-   Various minor bug fixes.

### \[5.2.4\] - 2018-09-10

#### Added

-   New command line flag --shadow-prefix

#### Changed

-   Fixed permission issue when using the script directive. This is a
    breaking change for scripts referring to files relative to the
    script directory (see the
    [docs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#external-scripts)).
-   Fixed various minor bugs and papercuts.
-   Allow URL to local git repo with wrapper directive
    (`git+file:///path/to/your/repo/path_to_file@@version`)

### \[5.2.2\] - 2018-08-01

#### Changed

-   Always print timestamps, removed the --timestamps CLI option.
-   more robust detection of conda command
-   Fixed bug in RMarkdown script execution.
-   Fixed a bug in detection of group jobs.

## \[5.2.0\] - 2018-06-28

#### Changed

-   Directory outputs have to marked with `directory`. This ensures
    proper handling of timestamps and cleanup. This is a breaking
    change. Implemented by Rasmus Ågren.
-   Fixed kubernetes tests, fixed kubernetes volume handling.
    Implemented by Andrew Schriefer.
-   jinja2 and networkx are not optional dependencies when installing
    via pip.
-   When conda or singularity directives are used and the corresponding
    CLI flags are not specified, the user is notified at the beginning
    of the log output.
-   Fixed numerous small bugs and papercuts and extended documentation.

### \[5.1.5\] - 2018-06-24

#### Changed

-   fixed missing version info in docker image.
-   several minor fixes to EGA support.

### \[5.1.4\] - 2018-05-28

#### Added

-   Allow `category` to be set.

#### Changed

-   Various cosmetic changes to reports.
-   Fixed encoding issues in reports.

### \[5.1.3\] - 2018-05-22

#### Changed

-   Fixed various bugs in job groups, shadow directive, singularity
    directive, and more.

### \[5.1.2\] - 2018-05-18

#### Changed

-   Fixed a bug in the report stylesheet.

## \[5.1.0\] - 2018-05-17

#### Added

-   A new framework for self-contained HTML reports, including results,
    statistics and topology information. In future releases this will be
    further extended.
-   A new utility snakemake.utils.validate() which allows to validate
    config and pandas data frames using JSON schemas.
-   Two new flags --cleanup-shadow and --cleanup-conda to clean up old
    unused conda and shadow data.

#### Changed

-   Benchmark repeats are now specified inside the workflow via a new
    flag repeat().
-   Command line interface help has been refactored into groups for
    better readability.

## \[5.0.0\] - 2018-05-11

#### Added

-   Group jobs for reduced queuing and network overhead, in particular
    with short running jobs.
-   Output files can be marked as pipes, such that producing and
    consuming job are executed simultaneously and interfomation is
    transferred directly without using disk.
-   Command line flags to clean output files.
-   Command line flag to list files in working directory that are not
    tracked by Snakemake.

#### Changed

-   Fix of --default-remote-prefix in case of input functions returning
    lists or dicts.
-   Scheduler no longer prefers jobs with many downstream jobs.

### \[4.8.1\] - 2018-04-25

#### Added

-   Allow URLs for the conda directive. # Changed
-   Various minor updates in the docs.
-   Several bug fixes with remote file handling.
-   Fix ImportError occuring with script directive.
-   Use latest singularity.
-   Improved caching for file existence checks. We first check existence
    of parent directories and cache these results. By this, large parts
    of the generated FS tree can be pruned if files are not yet present.
    If files are present, the overhead is minimal, since the checks for
    the parents are cached.
-   Various minor bug fixes.

## \[4.8.0\] - 2018-03-13

#### Added

-   Integration with CWL: the `cwl` directive allows to use CWL tool
    definitions in addition to shell commands or Snakemake wrappers.
-   A global `singularity` directive allows to define a global
    singularity container to be used for all rules that don't specify
    their own.
-   Singularity and Conda can now be combined. This can be used to
    specify the operating system (via singularity), and the software
    stack (via conda), without the overhead of creating specialized
    container images for workflows or tasks.

## \[4.7.0\] - 2018-02-19

#### Changed

-   Speedups when calculating dry-runs.
-   Speedups for workflows with many rules when calculating the DAG.
-   Accept SIGTERM to gracefully finish all running jobs and exit.
-   Various minor bug fixes.

## \[4.6.0\] - 2018-02-06

#### Changed

-   Log files can now be used as input files for other rules.
-   Adapted to changes in Kubernetes client API.
-   Fixed minor issues in --archive option.
-   Search path order in scripts was changed to fix a bug with leaked
    packages from root env when using script directive together with
    conda.

### \[4.5.1\] - 2018-02-01

#### Added

-   Input and output files can now tag pathlib objects. # ## Changed
-   Various minor bug fixes.

## \[4.5.0\] - 2018-01-18

#### Added

-   iRODS remote provider # ## Changed
-   Bug fix in shell usage of scripts and wrappers.
-   Bug fixes for cluster execution, --immediate-submit and
    subworkflows.

## \[4.4.0\] - 2017-12-21

#### Added

-   A new shadow mode (minimal) that only symlinks input files has been
    added.

#### Changed

-   The default shell is now bash on linux and macOS. If bash is not
    installed, we fall back to sh. Previously, Snakemake used the
    default shell of the user, which defeats the purpose of portability.
    If the developer decides so, the shell can be always overwritten
    using shell.executable().
-   Snakemake now requires Singularity 2.4.1 at least (only when running
    with --use-singularity).
-   HTTP remote provider no longer automatically unpacks gzipped files.
-   Fixed various smaller bugs.

### \[4.3.1\] - 2017-11-16

#### Added

-   List all conda environments with their location on disk via
    --list-conda-envs.

#### Changed

-   Do not clean up shadow on dry-run.
-   Allow R wrappers.

## \[4.3.0\] - 2017-10-27

#### Added

-   GridFTP remote provider. This is a specialization of the GFAL remote
    provider that uses globus-url-copy to download or upload files. # ##
    Changed
-   Scheduling and execution mechanisms have undergone a major revision
    that removes several potential (but rare) deadlocks.
-   Several bugs and corner cases of the singularity support have been
    fixed.
-   Snakemake now requires singularity 2.4 at least.

## \[4.2.0\] - 2017-10-10

#### Added

-   Support for executing jobs in per-rule singularity images. This is
    meant as an alternative to the conda directive (see docs), providing
    even more guarantees for reproducibility.

#### Changed

-   In cluster mode, jobs that are still running after Snakemake has
    been killed are automatically resumed.
-   Various fixes to GFAL remote provider.
-   Fixed --summary and --list-code-changes.
-   Many other small bug fixes.

## \[4.1.0\] - 2017-09-26

#### Added

-   Support for configuration profiles. Profiles allow to specify
    default options, e.g., a cluster submission command. They can be
    used via 'snakemake --profile myprofile'. See the docs for details.
-   GFAL remote provider. This allows to use GridFTP, SRM and any other
    protocol supported by GFAL for remote input and output files.
-   Added --cluster-status flag that allows to specify a command that
    returns jobs status. # ## Changed
-   The scheduler now tries to get rid of the largest temp files first.
-   The Docker image used for kubernetes support can now be configured
    at the command line.
-   Rate-limiting for cluster interaction has been unified.
-   S3 remote provider uses boto3.
-   Resource functions can now use an additional `attempt` parameter,
    that contains the number of times this job has already been tried.
-   Various minor fixes.

## \[4.0.0\] - 2017-07-24

#### Added

-   Cloud computing support via Kubernetes. Snakemake workflows can be
    executed transparently in the cloud, while storing input and output
    files within the cloud storage (e.g. S3 or Google Storage). I.e.,
    this feature does not need a shared filesystem between the cloud
    notes, and thereby makes the setup really simple.
-   WebDAV remote file support: Snakemake can now read and write from
    WebDAV. Hence, it can now, e.g., interact with Nextcloud or
    Owncloud.
-   Support for default remote providers: define a remote provider to
    implicitly use for all input and output files.
-   Added an option to only create conda environments instead of
    executing the workflow. # ## Changed
-   The number of files used for the metadata tracking of Snakemake
    (e.g., code, params, input changes) in the .snakemake directory has
    been reduced by a factor of 10, which should help with NFS and IO
    bottlenecks. This is a breaking change in the sense that Snakemake
    4.x won't see the metadata of workflows executed with Snakemake 3.x.
    However, old metadata won't be overwritten, so that you can always
    go back and check things by installing an older version of Snakemake
    again.
-   The google storage (GS) remote provider has been changed to use the
    google SDK. This is a breaking change, since the remote provider
    invocation has been simplified (see docs).
-   Due to WebDAV support (which uses asyncio), Snakemake now requires
    Python 3.5 at least.
-   Various minor bug fixes (e.g. for dynamic output files).

### \[3.13.3\] - 2017-06-23

#### Changed

-   Fix a followup bug in Namedlist where a single item was not returned
    as string.

### \[3.13.2\] - 2017-06-20

#### Changed

-   The --wrapper-prefix flag now also affects where the corresponding
    environment definition is fetched from.
-   Fix bug where empty output file list was recognized as containing
    duplicates (issue #574).

### \[3.13.1\] - 2017-06-20

#### Changed

-   Fix --conda-prefix to be passed to all jobs.
-   Fix cleanup issue with scripts that fail to download.

## \[3.13.0\] - 2017-06-12

#### Added

-   An NCBI remote provider. By this, you can seamlessly integrate any
    NCBI resouce (reference genome, gene/protein sequences, ...) as
    input file. # ## Changed
-   Snakemake now detects if automatically generated conda environments
    have to be recreated because the workflow has been moved to a new
    path.
-   Remote functionality has been made more robust, in particular to
    avoid race conditions.
-   `--config` parameter evaluation has been fixed for non-string types.
-   The Snakemake docker container is now based on the official debian
    image.

## \[3.12.0\] - 2017-05-09

#### Added

-   Support for RMarkdown (.Rmd) in script directives.
-   New option --debug-dag that prints all decisions while building the
    DAG of jobs. This helps to debug problems like cycles or unexpected
    MissingInputExceptions.
-   New option --conda-prefix to specify the place where conda
    environments are stored.

#### Changed

-   Benchmark files now also include the maximal RSS and VMS size of the
    Snakemake process and all sub processes.
-   Speedup conda environment creation.
-   Allow specification of DRMAA log dir.
-   Pass cluster config to subworkflow.

### \[3.11.2\] - 2017-03-15

#### Changed

-   Fixed fix handling of local URIs with the wrapper directive.

### \[3.11.1\] - 2017-03-14

#### Changed

-   --touch ignores missing files
-   Fixed handling of local URIs with the wrapper directive.

## \[3.11.0\] - 2017-03-08

#### Added

-   Param functions can now also refer to threads. # ## Changed
-   Improved tutorial and docs.
-   Made conda integration more robust.
-   None is converted to NULL in R scripts.

### \[3.10.2\] - 2017-02-28

#### Changed

-   Improved config file handling and merging.
-   Output files can be referred in params functions (i.e. lambda
    wildcards, output: ...)
-   Improved conda-environment creation.
-   Jobs are cached, leading to reduced memory footprint.
-   Fixed subworkflow handling in input functions.

## \[3.10.0\] - 2017-01-18

#### Added

-   Workflows can now be archived to a tarball with
    `snakemake --archive my-workflow.tar.gz`. The archive contains all
    input files, source code versioned with git and all software
    packages that are defined via conda environments. Hence, the archive
    allows to fully reproduce a workflow on a different machine. Such an
    archive can be uploaded to Zenodo, such that your workflow is
    secured in a self-contained, executable way for the future. # ##
    Changed
-   Improved logging.
-   Reduced memory footprint.
-   Added a flag to automatically unpack the output of input functions.
-   Improved handling of HTTP redirects with remote files.
-   Improved exception handling with DRMAA.
-   Scripts referred by the script directive can now use locally defined
    external python modules.

### \[3.9.1\] - 2016-12-23

#### Added

-   Jobs can be restarted upon failure (--restart-times). # ## Changed
-   The docs have been restructured and improved. Now available under
    snakemake.readthedocs.org.
-   Changes in scripts show up with --list-code-changes.
-   Duplicate output files now cause an error.
-   Various bug fixes.

## \[3.9.0\] - 2016-11-15

#### Added

-   Ability to define isolated conda software environments (YAML) per
    rule. Environments will be deployed by Snakemake upon workflow
    execution.
-   Command line argument --wrapper-prefix in order to overwrite the
    default URL for looking up wrapper scripts. # ## Changed
-   --summary now displays the log files correspoding to each output
    file.
-   Fixed hangups when using run directive and a large number of jobs
-   Fixed pickling errors with anonymous rules and run directive.
-   Various small bug fixes

### \[3.8.2\] - 2016-09-23

#### Changed

-   Add missing import in rules.py.
-   Use threading only in cluster jobs.

### \[3.8.1\] - 2016-09-14

#### Changed

-   Snakemake now warns when using relative paths starting with "./".
-   The option -R now also accepts an empty list of arguments.
-   Bug fix when handling benchmark directive.
-   Jobscripts exit with code 1 in case of failure. This should improve
    the error messages of cluster system.
-   Fixed a bug in SFTP remote provider.

## \[3.8.0\] - 2016-08-26

#### Added

-   Wildcards can now be constrained by rule and globally via the new
    `wildcard_constraints` directive (see the
    [docs](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-wildcards)).
-   Subworkflows now allow to overwrite their config file via the
    configfile directive in the calling Snakefile.
-   A method `log_fmt_shell` in the snakemake proxy object that is
    available in scripts and wrappers allows to obtain a formatted
    string to redirect logging output from STDOUT or STDERR.
-   Functions given to resources can now optionally contain an
    additional argument `input` that refers to the input files.
-   Functions given to params can now optionally contain additional
    arguments `input` (see above) and `resources`. The latter refers to
    the resources.
-   It is now possible to let items in shell commands be automatically
    quoted (see the
    [docs](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-rules)).
    This is usefull when dealing with filenames that contain
    whitespaces.

#### Changed

-   Snakemake now deletes output files before job exection. Further, it
    touches output files after job execution. This solves various
    problems with slow NFS filesystems.
-   A bug was fixed that caused dynamic output rules to be executed
    multiple times when forcing their execution with -R.
-   A bug causing double uploads with remote files was fixed. Various
    additional bug fixes related to remote files.
-   Various minor bug fixes.

### \[3.7.1\] - 2016-05-16

#### Changed

-   Fixed a missing import of the multiprocessing module.

## \[3.7.0\] - 2016-05-05

#### Added

-   The entries in `resources` and the `threads` job attribute can now
    be callables that must return `int` values.
-   Multiple `--cluster-config` arguments can be given to the Snakemake
    command line. Later one override earlier ones.
-   In the API, multiple `cluster_config` paths can be given as a list,
    alternatively to the previous behaviour of expecting one string for
    this parameter.
-   When submitting cluster jobs (either through `--cluster` or
    `--drmaa`), you can now use `--max-jobs-per-second` to limit the
    number of jobs being submitted (also available through Snakemake
    API). Some cluster installations have problems with too many jobs
    per second.
-   Wildcard values are now printed upon job execution in addition to
    input and output files. # ## Changed
-   Fixed a bug with HTTP remote providers.

### \[3.6.1\] - 2016-04-08

#### Changed

-   Work around missing RecursionError in Python \< 3.5
-   Improved conversion of numpy and pandas data structures to R
    scripts.
-   Fixed locking of working directory.

## \[3.6.0\] - 2016-03-10

#### Added

-   onstart handler, that allows to add code that shall be only executed
    before the actual workflow execution (not on dryrun).
-   Parameters defined in the cluster config file are now accessible in
    the job properties under the key "cluster".
-   The wrapper directive can be considered stable. # ## Changed
-   Allow to use rule/job parameters with braces notation in cluster
    config.
-   Show a proper error message in case of recursion errors.
-   Remove non-empty temp dirs.
-   Don't set the process group of Snakemake in order to allow kill
    signals from parent processes to be propagated.
-   Fixed various corner case bugs.
-   The params directive no longer converts a list `l` implicitly to
    `" ".join(l)`.

### \[3.5.5\] - 2016-01-23

#### Added

-   New experimental wrapper directive, which allows to refer to
    re-usable [wrapper
    scripts](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-wrappers).
    Wrappers are provided in the [Snakemake Wrapper
    Repository](https://bitbucket.org/snakemake/snakemake-wrappers).
-   David Koppstein implemented two new command line options to
    constrain the execution of the DAG of job to sub-DAGs (--until and
    --omit-from). # ## Changed
-   Fixed various bugs, e.g. with shadow jobs and --latency-wait.

### \[3.5.4\] - 2015-12-04

#### Changed

-   The params directive now fully supports non-string parameters.
    Several bugs in the remote support were fixed.

### \[3.5.3\] - 2015-11-24

#### Changed

-   The missing remote module was added to the package.

### \[3.5.2\] - 2015-11-24

#### Added

-   Support for easy integration of external R and Python scripts via
    the new [script
    directive](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-external-scripts).
-   Chris Tomkins-Tinch has implemented support for remote files:
    Snakemake can now handle input and output files from Amazon S3,
    Google Storage, FTP, SFTP, HTTP and Dropbox.
-   Simon Ye has implemented support for sandboxing jobs with [shadow
    rules](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-shadow-rules).

#### Changed

-   Manuel Holtgrewe has fixed dynamic output files in combination with
    multiple wildcards.
-   It is now possible to add suffixes to all shell commands with
    shell.suffix("mysuffix").
-   Job execution has been refactored to spawn processes only when
    necessary, resolving several problems in combination with huge
    workflows consisting of thousands of jobs and reducing the memory
    footprint.
-   In order to reflect the new collaborative development model,
    Snakemake has moved from my personal bitbucket account to
    <http://snakemake.bitbucket.org>.

### \[3.4.2\] - 2015-09-12

#### Changed

-   Willem Ligtenberg has reduced the memory usage of Snakemake.
-   Per Unneberg has improved config file handling to provide a more
    intuitive overwrite behavior.
-   Simon Ye has improved the test suite of Snakemake and helped with
    setting up continuous integration via Codeship.
-   The cluster implementation has been rewritten to use only a single
    thread to wait for jobs. This avoids failures with large numbers of
    jobs.
-   Benchmarks are now writing tab-delimited text files instead of JSON.
-   Snakemake now always requires to set the number of jobs with -j when
    in cluster mode. Set this to a high value if your cluster does not
    have restrictions.
-   The Snakemake Conda package has been moved to the bioconda channel.
-   The handling of Symlinks was improved, which made a switch to Python
    3.3 as the minimum required Python version necessary.

### \[3.4.1\] - 2015-08-05

#### Changed

-   This release fixes a bug that caused named input or output files to
    always be returned as lists instead of single files.

## \[3.4\] - 2015-07-18

#### Added

-   This release adds support for executing jobs on clusters in
    synchronous mode (e.g. qsub -sync). Thanks to David Alexander for
    implementing this.
-   There is now vim syntax highlighting support (thanks to Jay
    Hesselberth).
-   Snakemake is now available as Conda package.

#### Changed

-   Lots of bugs have been fixed. Thanks go to e.g. David Koppstein,
    Marcel Martin, John Huddleston and Tao Wen for helping with useful
    reports and debugging.

See [here](https://bitbucket.org/snakemake/snakemake/wiki/News-Archive)
for older changes.
