# Changelog

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
