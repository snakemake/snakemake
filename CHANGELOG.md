# Changelog


## [8.25.3](https://github.com/snakemake/snakemake/compare/v8.25.2...v8.25.3) (2024-11-11)


### Bug Fixes

* correctly set params in bash scripts ([#3188](https://github.com/snakemake/snakemake/issues/3188)) ([07ddab5](https://github.com/snakemake/snakemake/commit/07ddab5c97ae7f7b1a758d827dc6545df92ea644))
* more robust inference of source path that properly respects hosting provider urls without losing release or commit information ([#3195](https://github.com/snakemake/snakemake/issues/3195)) ([bd8212b](https://github.com/snakemake/snakemake/commit/bd8212b8f1b78a704d596a9d64040a48bccc223c))
* When generating a DAG or HTML rulegraph, use consistent colours ([#3189](https://github.com/snakemake/snakemake/issues/3189)) ([5f651d2](https://github.com/snakemake/snakemake/commit/5f651d2cd6a5964ec8b490781aeee05b1cf691a3))

## [8.25.2](https://github.com/snakemake/snakemake/compare/v8.25.1...v8.25.2) (2024-11-05)


### Bug Fixes

* include conda pinnings, conda post deploy script, and env modules for detection of software stack changes and corresponding rerun triggers ([#3184](https://github.com/snakemake/snakemake/issues/3184)) ([2aeaa46](https://github.com/snakemake/snakemake/commit/2aeaa46a06a34155ee24615d1417da103b54d14f))

## [8.25.1](https://github.com/snakemake/snakemake/compare/v8.25.0...v8.25.1) (2024-11-01)


### Bug Fixes

* ensure correct topological order when touching group job outputs ([#3181](https://github.com/snakemake/snakemake/issues/3181)) ([5924a3e](https://github.com/snakemake/snakemake/commit/5924a3e3993cb65365536ca4d225a1eb03c2fcd9))
* ensure version agnostic robust pickling of pandas, polars and numpy data structures passed as params to Python scripts or notebooks ([#3175](https://github.com/snakemake/snakemake/issues/3175)) ([eb11137](https://github.com/snakemake/snakemake/commit/eb1113713cbc4e9232aed6d106bec3615fa48632))

## [8.25.0](https://github.com/snakemake/snakemake/compare/v8.24.1...v8.25.0) (2024-10-29)


### Features

* add first 5 rules to group name (used e.g. when naming cluster/cloud jobs or logfiles) ([#3168](https://github.com/snakemake/snakemake/issues/3168)) ([5657122](https://github.com/snakemake/snakemake/commit/56571220a96afda4edf4b0578c697c9e94f8f15c))
* allow to mark input files of rules as ancient via the API or command line interface (and thereby also via workflow specific profiles). Putting this into a workflow specific profile (or specifying as argument) allows to overrule rerun triggers caused by file modification dates where the user knows better. ([#3171](https://github.com/snakemake/snakemake/issues/3171)) ([6f3aed3](https://github.com/snakemake/snakemake/commit/6f3aed321b48293b632246d29dd1bedc98e3d3b3))


### Bug Fixes

* skip storage object when cloning flags for shadowed IO ([#3174](https://github.com/snakemake/snakemake/issues/3174)) ([d733fed](https://github.com/snakemake/snakemake/commit/d733fed19c0fd4501fcd66e8528a801c978b1b53))
* use permission safe copying when hidden conda files are already present in a workdir. This avoids problems in case multiple people use the same workdir and workflow. ([#3169](https://github.com/snakemake/snakemake/issues/3169)) ([c98b2e7](https://github.com/snakemake/snakemake/commit/c98b2e7f71a391b99bffc54770654c9d74538ddf))


### Documentation

* add tutorial references and small syntax fix ([#3172](https://github.com/snakemake/snakemake/issues/3172)) ([6bee12a](https://github.com/snakemake/snakemake/commit/6bee12afdeee1621b50c96ecca0f3b2d6c3dc140))

## [8.24.1](https://github.com/snakemake/snakemake/compare/v8.24.0...v8.24.1) (2024-10-23)


### Bug Fixes

* fix bug with --edit-notebook sessions causing output files marked as incomplete, fix bug leading to missing log file after edit notebook sessions ([#3162](https://github.com/snakemake/snakemake/issues/3162)) ([19c6c0a](https://github.com/snakemake/snakemake/commit/19c6c0ab36da88adc9598ac18c20961c311eba28))
* proper error message if conda info fails ([#3157](https://github.com/snakemake/snakemake/issues/3157)) ([4f99c20](https://github.com/snakemake/snakemake/commit/4f99c201b31ad17a9b4f3ddb5be4b80c6f6f9a1f))

## [8.24.0](https://github.com/snakemake/snakemake/compare/v8.23.2...v8.24.0) (2024-10-21)


### Features

* subsample jobs to speed-up scheduler ([#3112](https://github.com/snakemake/snakemake/issues/3112)) ([e10feef](https://github.com/snakemake/snakemake/commit/e10feef262324357f9219d89214e2fccb2b85752))


### Documentation

* addition of interction/visualization/reporting tutorial ([#3159](https://github.com/snakemake/snakemake/issues/3159)) ([1d94bd1](https://github.com/snakemake/snakemake/commit/1d94bd10dbbc799b70a7a1d81c659e25de5fa5e0))
* fix tutorial step numbering ([2d7b9e9](https://github.com/snakemake/snakemake/commit/2d7b9e9402662556424e9564213f2efc9fda33c6))

## [8.23.2](https://github.com/snakemake/snakemake/compare/v8.23.1...v8.23.2) (2024-10-18)


### Bug Fixes

* ignore derived parameters when storing job metadata and inferring rerun necessity of jobs (this gets rid of spurious triggers caused by e.g. changed resources, threads, remote storage configuration) ([#3154](https://github.com/snakemake/snakemake/issues/3154)) ([73ce212](https://github.com/snakemake/snakemake/commit/73ce212d89275a320f0488918029197ea4850163))
* more robust handling of input metadata that avoids storing user or type specific local paths, and properly handles pipe or service files ([#3150](https://github.com/snakemake/snakemake/issues/3150)) ([756dc70](https://github.com/snakemake/snakemake/commit/756dc7045f8c2159da981d7000271ff87e8b0d10))


### Documentation

* Fix missing spaces in CLI help text ([#3146](https://github.com/snakemake/snakemake/issues/3146)) ([6416c36](https://github.com/snakemake/snakemake/commit/6416c3616b3628a389093e1e023707dff2b992a6))

## [8.23.1](https://github.com/snakemake/snakemake/compare/v8.23.0...v8.23.1) (2024-10-16)


### Bug Fixes

* fix bug in code change detection leading to spurious code change reporting when relying on older snakemake metadata ([#3144](https://github.com/snakemake/snakemake/issues/3144)) ([922d6e8](https://github.com/snakemake/snakemake/commit/922d6e863ec6b7d92f53f02845c68f5374380edd))

## [8.23.0](https://github.com/snakemake/snakemake/compare/v8.22.0...v8.23.0) (2024-10-14)


### Features

* More robust parameter and code change detection and transparent reporting of detected changes ([#3140](https://github.com/snakemake/snakemake/issues/3140)) ([576f588](https://github.com/snakemake/snakemake/commit/576f58868975e6763648224dc17191273049f971)). For setups using remote storage, this will lead to missing metadata records because the key used for storing the metadata is now the remote storage query instead of the path to the local copy. The reason is that the latter can be user-specific and therefore previously could have let led to e.g. different rerun behavior for different users.


### Documentation

* update installation docs to reflect changes in conda ecosystem ([#3141](https://github.com/snakemake/snakemake/issues/3141)) ([e83b8aa](https://github.com/snakemake/snakemake/commit/e83b8aaeac65711c327a5c86f82e231af33cfbf8))

## [8.22.0](https://github.com/snakemake/snakemake/compare/v8.21.0...v8.22.0) (2024-10-13)


### Features

* switch from toposort to graphlib ([#3109](https://github.com/snakemake/snakemake/issues/3109)) ([91e875d](https://github.com/snakemake/snakemake/commit/91e875d43fcb1cff247f82c743f9e2216ad328d7))


### Bug Fixes

* configfile `group` and `group-components` were not being registered ([#3135](https://github.com/snakemake/snakemake/issues/3135)) ([4397c7d](https://github.com/snakemake/snakemake/commit/4397c7d343289da5c6b2902249e3d78db5ac301e))
* remove paramiko dependency as issue has been fixed ([#3110](https://github.com/snakemake/snakemake/issues/3110)) ([1b43250](https://github.com/snakemake/snakemake/commit/1b43250782aaf92910001e6c2db924969956b103))

## [8.21.0](https://github.com/snakemake/snakemake/compare/v8.20.7...v8.21.0) (2024-10-12)


### Features

* support for specifying conda envs as directories ([#3132](https://github.com/snakemake/snakemake/issues/3132)) ([c54c95d](https://github.com/snakemake/snakemake/commit/c54c95d62b5395d78ab26297e0b84fc3f24dd017))


### Bug Fixes

* better error handling in evaluate function ([#3129](https://github.com/snakemake/snakemake/issues/3129)) ([04fb97f](https://github.com/snakemake/snakemake/commit/04fb97ff8b10fafdaf668368cb29e468d5bc07a4))
* notebook execution for apptainer ([#3131](https://github.com/snakemake/snakemake/issues/3131)) ([2e382c4](https://github.com/snakemake/snakemake/commit/2e382c48122c53712fd8db21ff06ca0e6c877e60))

## [8.20.7](https://github.com/snakemake/snakemake/compare/v8.20.6...v8.20.7) (2024-10-09)


### Bug Fixes

* install conda in container image ([#3127](https://github.com/snakemake/snakemake/issues/3127)) ([afa7bad](https://github.com/snakemake/snakemake/commit/afa7bad1a72118dc609f1499560c0847edc36413))
* remote pre-command not installing all required storage plugins ([#3116](https://github.com/snakemake/snakemake/issues/3116)) ([d829bb7](https://github.com/snakemake/snakemake/commit/d829bb7a57380354b3b4785d0755268ea1bf55d0))

## [8.20.6](https://github.com/snakemake/snakemake/compare/v8.20.5...v8.20.6) (2024-10-07)


### Bug Fixes

* rely on conda &gt;24.7.1 for conda env deployment and deprecate mamba support ([#3121](https://github.com/snakemake/snakemake/issues/3121)) ([9ece2db](https://github.com/snakemake/snakemake/commit/9ece2dbc57cec5a77ce909839b4c786be8e76ae1))

## [8.20.5](https://github.com/snakemake/snakemake/compare/v8.20.4...v8.20.5) (2024-09-25)


### Bug Fixes

* fixed check for remote conda env pinning and post-deploy files; fixed conda env cleanup ([#3103](https://github.com/snakemake/snakemake/issues/3103)) ([4d0a7e9](https://github.com/snakemake/snakemake/commit/4d0a7e9b0e56592aeda16d4961fd6acb1aabcca6))
* omit storage downloads during dryrun in workflows with checkpoints ([#3100](https://github.com/snakemake/snakemake/issues/3100)) ([151216a](https://github.com/snakemake/snakemake/commit/151216a64d7b19ceabab400b8f345d3e02d442ab))

## [8.20.4](https://github.com/snakemake/snakemake/compare/v8.20.3...v8.20.4) (2024-09-20)


### Bug Fixes

* cache conda envs to fix performance regression introduced in https://github.com/snakemake/snakemake/pull/1300 ([#3093](https://github.com/snakemake/snakemake/issues/3093)) ([66600c4](https://github.com/snakemake/snakemake/commit/66600c4a293446f5fdfff25075f5844fb2d720f2))
* Flatten conda pip dependencies for report rule info ([#3085](https://github.com/snakemake/snakemake/issues/3085)) ([56a1f20](https://github.com/snakemake/snakemake/commit/56a1f207ecf8343deab2b1583709fc9effc0ffb1))
* improve runtime complexity of post-job checkpoint handling ([#3096](https://github.com/snakemake/snakemake/issues/3096)) ([ba30781](https://github.com/snakemake/snakemake/commit/ba307816220c4313f82c0362fd6f00e783972569))


### Documentation

* Clarify the lookup function docstring ([#3091](https://github.com/snakemake/snakemake/issues/3091)) ([94177d5](https://github.com/snakemake/snakemake/commit/94177d5043bb5b15bf507c470374cdef9e8fe78f))
* Update lookup signature ([#3090](https://github.com/snakemake/snakemake/issues/3090)) ([655d6a1](https://github.com/snakemake/snakemake/commit/655d6a121642d781ea44b93bd21d1f4a36f80c36))

## [8.20.3](https://github.com/snakemake/snakemake/compare/v8.20.2...v8.20.3) (2024-09-09)


### Bug Fixes

* Add --cache to general_args for jobscripts ([#3080](https://github.com/snakemake/snakemake/issues/3080)) ([884498b](https://github.com/snakemake/snakemake/commit/884498bf6d4e171c8b4cef1463f9490689f0528d))

## [8.20.2](https://github.com/snakemake/snakemake/compare/v8.20.1...v8.20.2) (2024-09-09)


### Bug Fixes

* ensure that CLI help is formatted deterministically and set defaults are displayed properly (based on preliminary work by [@keszybz](https://github.com/keszybz)) ([#3081](https://github.com/snakemake/snakemake/issues/3081)) ([cbc2e2c](https://github.com/snakemake/snakemake/commit/cbc2e2c4086426a52f71a15384d18973ec8c9714))

## [8.20.1](https://github.com/snakemake/snakemake/compare/v8.20.0...v8.20.1) (2024-09-07)


### Bug Fixes

* add testcase data to package files ([#3077](https://github.com/snakemake/snakemake/issues/3077)) ([e90b098](https://github.com/snakemake/snakemake/commit/e90b09835404c3112c7c455accd2999f53653723))

## [8.20.0](https://github.com/snakemake/snakemake/compare/v8.19.3...v8.20.0) (2024-09-07)


### Features

* allow skipping apptainer and conda tests if not available ([#3074](https://github.com/snakemake/snakemake/issues/3074)) ([c3cd4e0](https://github.com/snakemake/snakemake/commit/c3cd4e01ad182433cbd4406f9007c167f805c0f5))
* MANIFEST.in: include tests/ in sdist ([#3073](https://github.com/snakemake/snakemake/issues/3073)) ([e24f3fc](https://github.com/snakemake/snakemake/commit/e24f3fc3078e4cee56334e2101939dd146c8f1b3))


### Bug Fixes

* add subfolders of report template to package data ([14a8f22](https://github.com/snakemake/snakemake/commit/14a8f2244f8d602f5c2988937ed722ab24c77259))
* Inconsistent Linting Output Formatting ([#3064](https://github.com/snakemake/snakemake/issues/3064)) ([90e51ae](https://github.com/snakemake/snakemake/commit/90e51ae68985d539ee716aa619fca6a8b625c60e))
* retry when downloading assets ([58e41b0](https://github.com/snakemake/snakemake/commit/58e41b0f6d1e7f7a6c7e2d915ba83ea786b237ed))
* skip asset download if files are already present ([#3076](https://github.com/snakemake/snakemake/issues/3076)) ([0d3d1e1](https://github.com/snakemake/snakemake/commit/0d3d1e1fd254e28f956243f4e01b51eb81b9a717))

## [8.19.3](https://github.com/snakemake/snakemake/compare/v8.19.2...v8.19.3) (2024-09-05)


### Bug Fixes

* add subfolders of report template to package data ([14a8f22](https://github.com/snakemake/snakemake/commit/14a8f2244f8d602f5c2988937ed722ab24c77259))
* add template data to setup call ([099c8e2](https://github.com/snakemake/snakemake/commit/099c8e2a4a82975a8ae1352ce7204e3cd6a3e83c))

## [8.19.2](https://github.com/snakemake/snakemake/compare/v8.19.1...v8.19.2) (2024-09-05)


### Bug Fixes

* add template data to setup call ([099c8e2](https://github.com/snakemake/snakemake/commit/099c8e2a4a82975a8ae1352ce7204e3cd6a3e83c))

## [8.19.1](https://github.com/snakemake/snakemake/compare/v8.19.0...v8.19.1) (2024-09-04)


### Bug Fixes

* fix issues with misinterpretation of max-jobs-per-timespan and max-jobs-per-seconds ([#3067](https://github.com/snakemake/snakemake/issues/3067)) ([d82453b](https://github.com/snakemake/snakemake/commit/d82453b7e54321c817c6516fca9d9e16c03c3977))
* pip deployment path ([#3062](https://github.com/snakemake/snakemake/issues/3062)) ([bf9305b](https://github.com/snakemake/snakemake/commit/bf9305b643ba5267285f56503809f275997d5a2e))
* return empty set if rate limiter at max ([#3060](https://github.com/snakemake/snakemake/issues/3060)) ([4e59963](https://github.com/snakemake/snakemake/commit/4e599633d331aa4857a8994afabe0df15b4241be))
* use pulps internal timeout for large scheduling problems, instead of stopit  ([#2938](https://github.com/snakemake/snakemake/issues/2938)) ([3b64e41](https://github.com/snakemake/snakemake/commit/3b64e414250de209a3921e44827779b1911f964a))
* Wrong linenumbers reported when linting ([#2985](https://github.com/snakemake/snakemake/issues/2985)) ([3a8bd36](https://github.com/snakemake/snakemake/commit/3a8bd367f30020fb3cf7d5efe64ffc7b63548202))


### Documentation

* update `doc-environment.yml` file and Documentation Setup documentation ([#3058](https://github.com/snakemake/snakemake/issues/3058)) ([a540a2e](https://github.com/snakemake/snakemake/commit/a540a2ec083f403dc96246109f65f6fca91ec488))

## [8.19.0](https://github.com/snakemake/snakemake/compare/v8.18.2...v8.19.0) (2024-08-29)


### Features

* check consistency of output file mtimes (must be newer than input files) ([#3050](https://github.com/snakemake/snakemake/issues/3050)) ([666cf62](https://github.com/snakemake/snakemake/commit/666cf6271d8dc5525dbea238bb12e496e437f316))
* print host name when executing workflow ([#3048](https://github.com/snakemake/snakemake/issues/3048)) ([b0ff787](https://github.com/snakemake/snakemake/commit/b0ff787b10ccf809789229c65901f363cfbc7a44))


### Bug Fixes

* `mem` and `disk` inference fixes ([#3040](https://github.com/snakemake/snakemake/issues/3040)) ([7530794](https://github.com/snakemake/snakemake/commit/7530794e4a5c64309e225b80f52932b404690756))
* avoid error accessing superclass in ioutils ([#3056](https://github.com/snakemake/snakemake/issues/3056)) ([a66a5f5](https://github.com/snakemake/snakemake/commit/a66a5f5bed1ea6c920e1b4bec0b0690566375261))
* disable execute after print compilation ([#3041](https://github.com/snakemake/snakemake/issues/3041)) ([86ed3cd](https://github.com/snakemake/snakemake/commit/86ed3cd7144f8ac8c54d4c67fe1d6023506ba4c7))
* download report assets upon package build such that reports become possible offline (cont. of [#2904](https://github.com/snakemake/snakemake/issues/2904)) ([#3026](https://github.com/snakemake/snakemake/issues/3026)) ([e8dad4b](https://github.com/snakemake/snakemake/commit/e8dad4bf4033b85e548cf5f8e37cb0443d6a959e))


### Documentation

* Add 'Editor integrations' section to Installation page ([#3045](https://github.com/snakemake/snakemake/issues/3045)) ([9a4006d](https://github.com/snakemake/snakemake/commit/9a4006d2e02eb64244001990f685d488756a873a))
* Fix typo (seesee to see) ([#3037](https://github.com/snakemake/snakemake/issues/3037)) ([de201fb](https://github.com/snakemake/snakemake/commit/de201fb0740123c9d9d30b6b929600062de485ff))
* various documentation fixes ([#3052](https://github.com/snakemake/snakemake/issues/3052)) ([b11460c](https://github.com/snakemake/snakemake/commit/b11460c0110bec6176bc44250336166602ffb2ed))

## [8.18.2](https://github.com/snakemake/snakemake/compare/v8.18.1...v8.18.2) (2024-08-21)


### Documentation

* recommending raw strings to get rid of syntax warnings ([#3022](https://github.com/snakemake/snakemake/issues/3022)) ([877b3a3](https://github.com/snakemake/snakemake/commit/877b3a3172ba1b5e86ec6b10d4c5b47bf2a96407))
* tutorial polishing ([16b1657](https://github.com/snakemake/snakemake/commit/16b16578d8bbf3d927ff26e4ff7f7434964ea38c))

## [8.18.1](https://github.com/snakemake/snakemake/compare/v8.18.0...v8.18.1) (2024-08-19)


### Bug Fixes

* add assets and use local file links to allow offline reports ([#2904](https://github.com/snakemake/snakemake/issues/2904)) ([9cd94f7](https://github.com/snakemake/snakemake/commit/9cd94f71077f4f91d566daa67d54ca40e1c5b276))
* use query from storage object in order to be able to reflect possible modifications (via StorageProvider.postprocess_query()) ([#3031](https://github.com/snakemake/snakemake/issues/3031)) ([3ddae58](https://github.com/snakemake/snakemake/commit/3ddae58e53eaf7020823a2a37145674d9a432389))


### Documentation

* clarify config file location ([6bd67d7](https://github.com/snakemake/snakemake/commit/6bd67d7b3284f054b893d72a3fc096ac711f54b4))

## [8.18.0](https://github.com/snakemake/snakemake/compare/v8.17.0...v8.18.0) (2024-08-14)


### Features

* show info on missing metadata ([#3014](https://github.com/snakemake/snakemake/issues/3014)) ([e502312](https://github.com/snakemake/snakemake/commit/e502312e72ad8dd1f54424eb95a7a0951e9fac62))


### Performance Improvements

* cache mtime of scripts and notebooks ([#2965](https://github.com/snakemake/snakemake/issues/2965)) ([405a056](https://github.com/snakemake/snakemake/commit/405a056230ba20114eac6caf083953148bb6d1c1))

## [8.17.0](https://github.com/snakemake/snakemake/compare/v8.16.0...v8.17.0) (2024-08-13)


### Features

* fix job rate limiting with --max-jobs-per-second and introduce the more flexible --max-jobs-per-timespan ([#3010](https://github.com/snakemake/snakemake/issues/3010)) ([9c31257](https://github.com/snakemake/snakemake/commit/9c3125702a4f24261948178d25dce0a2fea27466))


### Bug Fixes

* Allow hyphens in config keys given on the command line. ([#2998](https://github.com/snakemake/snakemake/issues/2998)) ([b70c0db](https://github.com/snakemake/snakemake/commit/b70c0db30212d84d9086954c34226528830bb4d9))
* allowing trailing '+' in name patterns ([#3002](https://github.com/snakemake/snakemake/issues/3002)) ([59150d3](https://github.com/snakemake/snakemake/commit/59150d3d7c55be7527369bdb677d8844f2e7d979))
* print message if not yet enough resources for executing further jobs ([b8df036](https://github.com/snakemake/snakemake/commit/b8df0364658071855b3dd66193aeceb4e80d26bf))
* unawaited coroutine sanitize_local_storage_copies ([#2972](https://github.com/snakemake/snakemake/issues/2972)) ([715c572](https://github.com/snakemake/snakemake/commit/715c57260316d2ba9c9e60f7dd4f1ed1b2a6c9f6))


### Documentation

* Change sha256 checksum in docs to more realistic example ([#2987](https://github.com/snakemake/snakemake/issues/2987)) ([16a5cf2](https://github.com/snakemake/snakemake/commit/16a5cf272959deefbe8919e6d3c8569d4b325991))
* Make it more clear that the cluster commands now require a plugin ([#2976](https://github.com/snakemake/snakemake/issues/2976)) ([74134cf](https://github.com/snakemake/snakemake/commit/74134cff09ae2e6ff725659ee1fe1fc9322e5a70))
* Update installation.rst to recommend Miniforge instead of Mambaforge ([#2975](https://github.com/snakemake/snakemake/issues/2975)) ([0fc7619](https://github.com/snakemake/snakemake/commit/0fc761998f31ab148bc7b1be7f33e21bab452441))
* use plain monospace font instead of theme default that changes &gt;= into ≥ ([cc17fc1](https://github.com/snakemake/snakemake/commit/cc17fc11cf83b03f9c528bab279d9f394789fb55))

## [8.16.0](https://github.com/snakemake/snakemake/compare/v8.15.2...v8.16.0) (2024-07-09)


### Features

* added snakemake.script.snakemake for type hinting ([#2917](https://github.com/snakemake/snakemake/issues/2917)) ([c85fb4b](https://github.com/snakemake/snakemake/commit/c85fb4b34ce983cb916bcce361d958d466ee1d07))


### Documentation

* use note directive to align with the capabilities of the used theme ([#2950](https://github.com/snakemake/snakemake/issues/2950)) ([fe27405](https://github.com/snakemake/snakemake/commit/fe274055dbd4663d7185114bc4f77ac862f8c7b3))

## [8.15.2](https://github.com/snakemake/snakemake/compare/v8.15.1...v8.15.2) (2024-07-05)


### Bug Fixes

* ensure that envvars in local storage prefix are not prematurely expanded by the shell ([#2943](https://github.com/snakemake/snakemake/issues/2943)) ([da50f27](https://github.com/snakemake/snakemake/commit/da50f27102e0aa761aa6ff724efe785af7a58b3a))
* fix circular import ([9e7d56f](https://github.com/snakemake/snakemake/commit/9e7d56f48d6048d8138c718eb577f65e157b42d2))

## [8.15.1](https://github.com/snakemake/snakemake/compare/v8.15.0...v8.15.1) (2024-07-04)


### Bug Fixes

* implement support for --touch on remote storage (if the storage provider supports it) ([#2941](https://github.com/snakemake/snakemake/issues/2941)) ([567094d](https://github.com/snakemake/snakemake/commit/567094d63242d6a885328d40f66bcd98204411a1))

## [8.15.0](https://github.com/snakemake/snakemake/compare/v8.14.0...v8.15.0) (2024-07-04)


### Features

* add `can_transfer_local_files` to executor plugin interface ([#2921](https://github.com/snakemake/snakemake/issues/2921)) ([85a6774](https://github.com/snakemake/snakemake/commit/85a6774542b110e16b9e56407a99995046eccfe1))


### Bug Fixes

* duplicate wildcards ([#2937](https://github.com/snakemake/snakemake/issues/2937)) ([5b6cc02](https://github.com/snakemake/snakemake/commit/5b6cc0283646a413ba8b9f71dea610faa3128f96))
* handling of missing attributes in input/output/params lists that have been guarded against misuse (sort, index) ([#2928](https://github.com/snakemake/snakemake/issues/2928)) ([1b75087](https://github.com/snakemake/snakemake/commit/1b75087fa10cfcfb06e41a6286cb51d10df3ee3f))
* improve error message in case of invalid default storage prefix ([3c1065c](https://github.com/snakemake/snakemake/commit/3c1065c9b441928f2d390b568a7e9f854256ad41))
* include more context to syntax errors when rule parsing fails ([#2924](https://github.com/snakemake/snakemake/issues/2924)) ([#2926](https://github.com/snakemake/snakemake/issues/2926)) ([3ecffcc](https://github.com/snakemake/snakemake/commit/3ecffccc7ea6df042765b82fe9f02a45b4daf727))
* parse f-string as-is from source ([#2930](https://github.com/snakemake/snakemake/issues/2930)) ([39fc8f7](https://github.com/snakemake/snakemake/commit/39fc8f71ce196008e938d5c549f4077cffc7682a))
* sanitize old local storage copies before evaluating params or resource functions ([#2939](https://github.com/snakemake/snakemake/issues/2939)) ([e3075dd](https://github.com/snakemake/snakemake/commit/e3075ddd4d9829df4732c80b372a249809710f95))
* suppress printing of non-constant or uninteresting defaults in --help ([#2936](https://github.com/snakemake/snakemake/issues/2936)) ([69add30](https://github.com/snakemake/snakemake/commit/69add3046a929dd7a3b5dc67d67b886dfde757ed))


### Documentation

* fix code block formatting ([#2929](https://github.com/snakemake/snakemake/issues/2929)) ([54a6461](https://github.com/snakemake/snakemake/commit/54a646144521d7799b7f810d6f5ee138ac25ca68))
* fix signature display ([bdd732d](https://github.com/snakemake/snakemake/commit/bdd732d1f98d0af0362115f74baa728bfc976851))

## [8.14.0](https://github.com/snakemake/snakemake/compare/v8.13.0...v8.14.0) (2024-06-11)


### Features

* add flag to mark files where path should not be modified ([#2888](https://github.com/snakemake/snakemake/issues/2888)) ([d142b46](https://github.com/snakemake/snakemake/commit/d142b46200df0e1d03b4de3ec7aafd9cdf714588))
* support per rule shell exec setting via resources ([#2862](https://github.com/snakemake/snakemake/issues/2862)) ([ab8d2dd](https://github.com/snakemake/snakemake/commit/ab8d2ddc38b610b0bf5b2e5368cfaf349921470a))


### Documentation

* update FAQ - recommend ensure function for failing on empty output ([#2910](https://github.com/snakemake/snakemake/issues/2910)) ([b035071](https://github.com/snakemake/snakemake/commit/b03507119806f1ada71f950086256c4fe63eb754))

## [8.13.0](https://github.com/snakemake/snakemake/compare/v8.12.0...v8.13.0) (2024-06-05)


### Features

* support for default value specification when using lookup helper function ([#2907](https://github.com/snakemake/snakemake/issues/2907)) ([08e88e2](https://github.com/snakemake/snakemake/commit/08e88e27f9e5ddcd5d9ebce7f79f45436262fddc))


### Documentation

* add badge for bioconda version ([#2902](https://github.com/snakemake/snakemake/issues/2902)) ([3f01348](https://github.com/snakemake/snakemake/commit/3f01348cda9040adc7af6c8173f9108343edde10))

## [8.12.0](https://github.com/snakemake/snakemake/compare/v8.11.6...v8.12.0) (2024-05-27)


### Features

* Include parameters in extended benchmarks ([#2887](https://github.com/snakemake/snakemake/issues/2887)) ([31a9c9b](https://github.com/snakemake/snakemake/commit/31a9c9bdb3147fe9c1f9144a584af389ec49720e))


### Bug Fixes

* fix corner case bug in input function exception handling ([#2895](https://github.com/snakemake/snakemake/issues/2895)) ([fc24292](https://github.com/snakemake/snakemake/commit/fc24292f14fd5f8bb00412c32568645cfbf5d2fa))
* fix quoting issues when passing complex apptainer args to spawned jobs ([#2898](https://github.com/snakemake/snakemake/issues/2898)) ([b07e2e0](https://github.com/snakemake/snakemake/commit/b07e2e07aaf7b9bd6fa6368e274f08028cdcea4e))
* properly restrict scheduler if --jobs/-j is used (in contrast to --cores) in local execution ([#2897](https://github.com/snakemake/snakemake/issues/2897)) ([6a276bb](https://github.com/snakemake/snakemake/commit/6a276bb25cc36d32285a196e0876e2b15040a863))
* typo in workflow specific profile default location ([#2878](https://github.com/snakemake/snakemake/issues/2878)) ([74627d3](https://github.com/snakemake/snakemake/commit/74627d3f07f549600dc8ad226a545f44191141e0))

## [8.11.6](https://github.com/snakemake/snakemake/compare/v8.11.5...v8.11.6) (2024-05-17)


### Bug Fixes

* fix opening of multiple checkout output files in the same input function when using remote storage ([2f8e719](https://github.com/snakemake/snakemake/commit/2f8e71927edfe636492bc14dfd7b90d3c226294c))

## [8.11.5](https://github.com/snakemake/snakemake/compare/v8.11.4...v8.11.5) (2024-05-16)


### Bug Fixes

* avoid premature deletion of local copies of remote storage input files used by multiple jobs ([#2874](https://github.com/snakemake/snakemake/issues/2874)) ([21ec649](https://github.com/snakemake/snakemake/commit/21ec649b6d69ff9eeeec229e1eb6058f417bba56))
* fix opening of checkpoint output files from remote storage ([#2873](https://github.com/snakemake/snakemake/issues/2873)) ([e7cb7fb](https://github.com/snakemake/snakemake/commit/e7cb7fb3e469057afc1f2c74008a40a716e76111))


### Documentation

* add link to code of conduct ([889a3bc](https://github.com/snakemake/snakemake/commit/889a3bca881e913764e4a4dc0eea70ee7ace8598))

## [8.11.4](https://github.com/snakemake/snakemake/compare/v8.11.3...v8.11.4) (2024-05-11)


### Bug Fixes

* fix missing await when opening checkpoint output ([#2868](https://github.com/snakemake/snakemake/issues/2868)) ([25a361b](https://github.com/snakemake/snakemake/commit/25a361bb14fcd66a26535093e0c05b5cf84e69d8))
* make checkpoint updates synchronous ([#2871](https://github.com/snakemake/snakemake/issues/2871)) ([b0e7ebd](https://github.com/snakemake/snakemake/commit/b0e7ebde7d95ac13b517ee4f200c92ebb1b1b805))


### Documentation

* update code of conduct email address ([3047683](https://github.com/snakemake/snakemake/commit/3047683bed9e28f9e5abc0d56f981aff18fbc801))

## [8.11.3](https://github.com/snakemake/snakemake/compare/v8.11.2...v8.11.3) (2024-05-03)


### Bug Fixes

* ignore errors when cleaning up runtime cache ([#2859](https://github.com/snakemake/snakemake/issues/2859)) ([6df7046](https://github.com/snakemake/snakemake/commit/6df70468bb29a9a5b172caa5fe5a1a8bde3f2ebe))
* show queries of remote storage files instead of local paths in summary ([#2860](https://github.com/snakemake/snakemake/issues/2860)) ([ba1db8e](https://github.com/snakemake/snakemake/commit/ba1db8eaf34a37a4d6029dd1bcddea49a775b912))

## [8.11.2](https://github.com/snakemake/snakemake/compare/v8.11.1...v8.11.2) (2024-05-02)


### Bug Fixes

* bug when requesting extended benchmark with slurm ([#2855](https://github.com/snakemake/snakemake/issues/2855)) ([0e039ff](https://github.com/snakemake/snakemake/commit/0e039ff5d0a0880c35297a9fd5b983aa2813226b))

## [8.11.1](https://github.com/snakemake/snakemake/compare/v8.11.0...v8.11.1) (2024-04-29)


### Bug Fixes

* check template rendering output for leaked input file paths ([#2850](https://github.com/snakemake/snakemake/issues/2850)) ([433302e](https://github.com/snakemake/snakemake/commit/433302ee990787d1aee5ad6d7cdbcc9533646305))
* do not distinguish between local and remote rules in dryrun ([74b99ec](https://github.com/snakemake/snakemake/commit/74b99ecee63a9922c4dbb6951cc63865955d198a))
* omit norun jobs when determining remote storage input file retrieval ([#2854](https://github.com/snakemake/snakemake/issues/2854)) ([37a7c7f](https://github.com/snakemake/snakemake/commit/37a7c7f9dbd7295e558298097301827f714db994))
* Prevent binary log files to crash snakemake execution with `show-failed-logs` ([#2827](https://github.com/snakemake/snakemake/issues/2827)) ([8a80bda](https://github.com/snakemake/snakemake/commit/8a80bdaed0d8524dffbdce71dd07eb123f71726d))
* replace pkg_resources for python 3.12 ([#2831](https://github.com/snakemake/snakemake/issues/2831)) ([ac144fc](https://github.com/snakemake/snakemake/commit/ac144fc7c68083cc56d1960f46b5f5e7888dd38e))
* return set instead of list when just --quiet ([#2829](https://github.com/snakemake/snakemake/issues/2829)) ([eeb57e2](https://github.com/snakemake/snakemake/commit/eeb57e26c4d85b520d5ae69a0f00183d4ea80eb8))
* small typo in error ([#2853](https://github.com/snakemake/snakemake/issues/2853)) ([325a715](https://github.com/snakemake/snakemake/commit/325a715b11f1763d835c67f70e2d301e903c7ebb))
* use keyword arguments in `_IOFile.open` ([#2847](https://github.com/snakemake/snakemake/issues/2847)) ([50c84dc](https://github.com/snakemake/snakemake/commit/50c84dc8a3b39143ea53edb77081448762832678))


### Documentation

* fix typo and link for RO Crate ([#2851](https://github.com/snakemake/snakemake/issues/2851)) ([cec0041](https://github.com/snakemake/snakemake/commit/cec0041a71cbd4ebb4e1f8ec4ca296127a65b24b))

## [8.11.0](https://github.com/snakemake/snakemake/compare/v8.10.8...v8.11.0) (2024-04-25)


### Features

* allow for more extensive benchmark file in jsonl format ([#2691](https://github.com/snakemake/snakemake/issues/2691)) ([de12463](https://github.com/snakemake/snakemake/commit/de12463604bc5f24f81907e33057ff42c31c7fc2))


### Bug Fixes

* only download input for local jobs in the main process, not within remote groups ([#2842](https://github.com/snakemake/snakemake/issues/2842)) ([97f428b](https://github.com/snakemake/snakemake/commit/97f428ba2921cebc2a311a4933c1c55536d34335))
* remove non-empty local copies of remote storage dirs ([#2845](https://github.com/snakemake/snakemake/issues/2845)) ([71b2b87](https://github.com/snakemake/snakemake/commit/71b2b87644ac0470aba1a96093e88e9a64df1176))
* retrieve files from storage if necessary when calling their open method (e.g. when accessing output files from a checkpoint) ([#2839](https://github.com/snakemake/snakemake/issues/2839)) ([5448208](https://github.com/snakemake/snakemake/commit/54482089893d12a0cab20a7b4b3d0af5a87aa772))

## [8.10.8](https://github.com/snakemake/snakemake/compare/v8.10.7...v8.10.8) (2024-04-19)


### Bug Fixes

* extend workflow test suite by wildcards with slashes (in order to detect bugs that can occur in executors) ([#2810](https://github.com/snakemake/snakemake/issues/2810)) ([fc9971b](https://github.com/snakemake/snakemake/commit/fc9971b227e43601c8b511c1c37e4dce3211e913))
* invalid extensions with multixt and cache ([#1808](https://github.com/snakemake/snakemake/issues/1808)) ([#2823](https://github.com/snakemake/snakemake/issues/2823)) ([a845b1c](https://github.com/snakemake/snakemake/commit/a845b1c62b106fb4c0f1b97485a3044bfbcab65b))
* module prefix added twice on expand ([#2814](https://github.com/snakemake/snakemake/issues/2814)) ([27416f4](https://github.com/snakemake/snakemake/commit/27416f4e233662162232de33e9b82205f7efe397))
* retrieve inputs of local rules from storage even if they are intermediate results ([#2811](https://github.com/snakemake/snakemake/issues/2811)) ([d158aa8](https://github.com/snakemake/snakemake/commit/d158aa849fb16fb3931d05fb052684c5d09862cd))


### Documentation

* add snakemake and reportplugins to req.txt to show all cli options ([#2817](https://github.com/snakemake/snakemake/issues/2817)) ([aaa9065](https://github.com/snakemake/snakemake/commit/aaa906565d102774b99b1515be0747d59e1828a6))

## [8.10.7](https://github.com/snakemake/snakemake/compare/v8.10.6...v8.10.7) (2024-04-12)


### Bug Fixes

* expand error on modules with prefix ([#1614](https://github.com/snakemake/snakemake/issues/1614)) ([#2803](https://github.com/snakemake/snakemake/issues/2803)) ([c6bfcd7](https://github.com/snakemake/snakemake/commit/c6bfcd762c595c7d009fb3c2cfff6231a1f3ffe4))
* failure message in cleanup_metadata ([#2800](https://github.com/snakemake/snakemake/issues/2800)) ([4c7fe55](https://github.com/snakemake/snakemake/commit/4c7fe55285fd4641d86a945110e546ad19dd9c2e))
* handle async method in cleanup_storage_objects ([#2806](https://github.com/snakemake/snakemake/issues/2806)) ([272205b](https://github.com/snakemake/snakemake/commit/272205beb0765e65c86a7b4e4f740b84ceac44b5))


### Documentation

* add Morten Lund as another maintainer. ([293bc05](https://github.com/snakemake/snakemake/commit/293bc056b6c44f17340543b2678f8a59e808a5a6))

## [8.10.6](https://github.com/snakemake/snakemake/compare/v8.10.5...v8.10.6) (2024-04-04)


### Bug Fixes

* only constrain by --max-threads if threads of job are already known ([#2790](https://github.com/snakemake/snakemake/issues/2790)) ([5f28fcd](https://github.com/snakemake/snakemake/commit/5f28fcd32bd8fa0f539c4fc8aabf07702787132b))

## [8.10.5](https://github.com/snakemake/snakemake/compare/v8.10.4...v8.10.5) (2024-04-04)


### Bug Fixes

* properly delete local copies of storage files after remote jobs ([#2793](https://github.com/snakemake/snakemake/issues/2793)) ([e3362b0](https://github.com/snakemake/snakemake/commit/e3362b03e11d9770ed6b01e82017973f4a652d4a))
* respect APPTAINER_CACHEDIR and allow env variables in --apptainer-prefix and --conda-prefix ([#2795](https://github.com/snakemake/snakemake/issues/2795)) ([b1694cd](https://github.com/snakemake/snakemake/commit/b1694cdf10ce46ff80236a8bb1d9c72e2f0fae0e))

## [8.10.4](https://github.com/snakemake/snakemake/compare/v8.10.3...v8.10.4) (2024-03-27)


### Documentation

* fix links ([420ea9c](https://github.com/snakemake/snakemake/commit/420ea9ca52d3de76f93e3ad30d922a97950c75de))

## [8.10.3](https://github.com/snakemake/snakemake/compare/v8.10.2...v8.10.3) (2024-03-27)


### Documentation

* mention maintainers ([be95e25](https://github.com/snakemake/snakemake/commit/be95e25d7ee65abe993f08b9d6ca726ee79ec7d4))

## [8.10.2](https://github.com/snakemake/snakemake/compare/v8.10.1...v8.10.2) (2024-03-26)


### Bug Fixes

* remove default packages from conda envs ([#2749](https://github.com/snakemake/snakemake/issues/2749)) ([027906c](https://github.com/snakemake/snakemake/commit/027906c519b5e39a628dda072f314e7b9f6343dc))
* use base64 encoding for passing default resources args to jobs ([#2780](https://github.com/snakemake/snakemake/issues/2780)) ([4735bc3](https://github.com/snakemake/snakemake/commit/4735bc3b81d9db96da0fb949893285ec0f77d076))

## [8.10.1](https://github.com/snakemake/snakemake/compare/v8.10.0...v8.10.1) (2024-03-26)


### Bug Fixes

* passing of --set-threads values to remote jobs ([#2775](https://github.com/snakemake/snakemake/issues/2775)) ([4fd767a](https://github.com/snakemake/snakemake/commit/4fd767a4e57043551e0551542a93dbc9d34df777))
* use base64 encoding when passing resources and threads to remote jobs (this solves issues with complex quoted resources) ([#2778](https://github.com/snakemake/snakemake/issues/2778)) ([a8ee4d8](https://github.com/snakemake/snakemake/commit/a8ee4d8a8228aeac7bccefdc3b868059455e38c4))

## [8.10.0](https://github.com/snakemake/snakemake/compare/v8.9.0...v8.10.0) (2024-03-22)


### Features

* expose ResourceSettings in TestWorkflowsBase ([#2770](https://github.com/snakemake/snakemake/issues/2770)) ([e7c323b](https://github.com/snakemake/snakemake/commit/e7c323b707091d759592a6c0b75e8f772b2c72c5))

## [8.9.0](https://github.com/snakemake/snakemake/compare/v8.8.0...v8.9.0) (2024-03-18)


### Features

* add function 'exists' for checking the prior existence of files or dirs before workflow execution while considering any remote storage settings. In addition: some bug fixes for error handling and the update/before_update functionality. ([ee96393](https://github.com/snakemake/snakemake/commit/ee9639385cb2c9cec425800088faee4e4bf77c9a))

## [8.8.0](https://github.com/snakemake/snakemake/compare/v8.7.0...v8.8.0) (2024-03-15)


### Features

* Allow smart_open 7.x ([#2745](https://github.com/snakemake/snakemake/issues/2745)) ([d52c1b1](https://github.com/snakemake/snakemake/commit/d52c1b1805d2e30124da2ac01f614fe4f5cb0257))


### Bug Fixes

* various error handling improvements, fixed logging/error behavior (stdout from dryrun, stderr otherwise) ([#2759](https://github.com/snakemake/snakemake/issues/2759)) ([d0d1f48](https://github.com/snakemake/snakemake/commit/d0d1f48e9bc1f8f994c2f243255b84a1b73e5208))

## [8.7.0](https://github.com/snakemake/snakemake/compare/v8.6.0...v8.7.0) (2024-03-13)


### Features

* add flag for marking output as being updated instead of rewritten (update("test.txt")) ([#2754](https://github.com/snakemake/snakemake/issues/2754)) ([9ba5d95](https://github.com/snakemake/snakemake/commit/9ba5d95a657952f06ad22cc415670f79c013d3f8))
* allow default storage provider to be explicitly set to none ([#2746](https://github.com/snakemake/snakemake/issues/2746)) ([ce519d7](https://github.com/snakemake/snakemake/commit/ce519d7fd2b9dde1f369cb0da27a71d910b46734))

## [8.6.0](https://github.com/snakemake/snakemake/compare/v8.5.5...v8.6.0) (2024-03-11)


### Features

* add setting for defining separate local storage prefix for remote jobs; improved ergonomics for semantic helper functions ([#2743](https://github.com/snakemake/snakemake/issues/2743)) ([5007e5c](https://github.com/snakemake/snakemake/commit/5007e5c038e0fa5e90a8d44ec4a75ec77c2cd91c))
* allow passing of lists of functions or single functions to expand ([#2741](https://github.com/snakemake/snakemake/issues/2741)) ([32e65df](https://github.com/snakemake/snakemake/commit/32e65df607127a93cc1857c415154137e3fd438d))


### Bug Fixes

* fix error message for invalid storage provider queries ([977951e](https://github.com/snakemake/snakemake/commit/977951ea541bceb97b6a77709fde863f6c638352))
* fix premature deletion of temp files in combination with checkpoints ([#2737](https://github.com/snakemake/snakemake/issues/2737)) ([b22ba5f](https://github.com/snakemake/snakemake/commit/b22ba5f32b893e8eb50f6e3a6f0874f1c932943b))

## [8.5.5](https://github.com/snakemake/snakemake/compare/v8.5.4...v8.5.5) (2024-03-07)


### Bug Fixes

* less frightering message when telling about missing output files as reason for running a job ([7b6c9d4](https://github.com/snakemake/snakemake/commit/7b6c9d4f620e020911192e0f311f1d9098391e0d))

## [8.5.4](https://github.com/snakemake/snakemake/compare/v8.5.3...v8.5.4) (2024-03-06)


### Bug Fixes

* fix bugs in --summary and --list-input-changes. Removed outdated statement in tutorial. ([#2735](https://github.com/snakemake/snakemake/issues/2735)) ([55c06d8](https://github.com/snakemake/snakemake/commit/55c06d853116cf402b90832649c805643712c437))


### Documentation

* fix explanation on how to use all cores in tutorial ([#2733](https://github.com/snakemake/snakemake/issues/2733)) ([6420428](https://github.com/snakemake/snakemake/commit/64204286f3061c494389ea502f1e4393e357efcd))
* fix syntax highlighting in tutorial ([2887604](https://github.com/snakemake/snakemake/commit/2887604ebe8a1e4c9c226780129312f556150dca))

## [8.5.3](https://github.com/snakemake/snakemake/compare/v8.5.2...v8.5.3) (2024-02-26)


### Bug Fixes

* error when detecting mime type during report creation ([#2721](https://github.com/snakemake/snakemake/issues/2721)) ([42dad42](https://github.com/snakemake/snakemake/commit/42dad4258b0b0028733b94cbde6ae63c534ba645))

## [8.5.2](https://github.com/snakemake/snakemake/compare/v8.5.1...v8.5.2) (2024-02-24)


### Bug Fixes

* when using remote storage: only wait for files if job did not error ([8c7ee91](https://github.com/snakemake/snakemake/commit/8c7ee91ac5425279955e60071fd0269bf2bf197b))


### Documentation

* fix apidocs theme ([1fb5467](https://github.com/snakemake/snakemake/commit/1fb546772bc90d96af2a390ced2d02db100ad810))

## [8.5.1](https://github.com/snakemake/snakemake/compare/v8.5.0...v8.5.1) (2024-02-24)


### Documentation

* api docs polishing ([58a6a13](https://github.com/snakemake/snakemake/commit/58a6a139a6fe2719de7d2d7c98e21f4d33a67cba))
* fix [#2698](https://github.com/snakemake/snakemake/issues/2698) ([#2714](https://github.com/snakemake/snakemake/issues/2714)) ([508080b](https://github.com/snakemake/snakemake/commit/508080b8e1185c45e4c3d7b413350c68444db90c))

## [8.5.0](https://github.com/snakemake/snakemake/compare/v8.4.12...v8.5.0) (2024-02-24)


### Features

* add ability to return input functions from input functions. Such nesting is evaluated 10 times at most. Beyond that, an error is thrown.  ([#2717](https://github.com/snakemake/snakemake/issues/2717)) ([7a47924](https://github.com/snakemake/snakemake/commit/7a47924627b30dc774efef0faae08e13e44c350a))
* support for report plugins ([#2700](https://github.com/snakemake/snakemake/issues/2700)) ([2f7d4b5](https://github.com/snakemake/snakemake/commit/2f7d4b5196979d368f6e2eede25ca2a61220691b))


### Bug Fixes

* fix wait for files in case of using remote storage and remote execution ([#2718](https://github.com/snakemake/snakemake/issues/2718)) ([eec3a5f](https://github.com/snakemake/snakemake/commit/eec3a5fd7a02c093dc14cfd0912fbaf8d7dd10cc))
* proper interpretation of standard resources given as strings (e.g. runtime as '5m'). Avoid the need to set additional quotes around size or timespan resources. Improved error messages for resource handling.  ([#2716](https://github.com/snakemake/snakemake/issues/2716)) ([b6636e9](https://github.com/snakemake/snakemake/commit/b6636e95f6d9baddb4a3eb5bc79351207f048b4d))

## [8.4.12](https://github.com/snakemake/snakemake/compare/v8.4.11...v8.4.12) (2024-02-20)


### Bug Fixes

* constrain dependencies to match conda experience ([#2710](https://github.com/snakemake/snakemake/issues/2710)) ([d9a7a13](https://github.com/snakemake/snakemake/commit/d9a7a13d6e427c025902a78c2d21c3f21cac242b))
* various bug fixes for resource parsing ([#2711](https://github.com/snakemake/snakemake/issues/2711)) ([d1daf0b](https://github.com/snakemake/snakemake/commit/d1daf0b4b2cc64b92adb9f58ac71c601ab5d29c8))

## [8.4.11](https://github.com/snakemake/snakemake/compare/v8.4.10...v8.4.11) (2024-02-19)


### Documentation

* fix heading ([1eda36d](https://github.com/snakemake/snakemake/commit/1eda36d9d7e545b05428548f8b826e941ad078d9))

## [8.4.10](https://github.com/snakemake/snakemake/compare/v8.4.9...v8.4.10) (2024-02-19)


### Bug Fixes

* properly handle --touch when using a storage provider ([#2705](https://github.com/snakemake/snakemake/issues/2705)) ([fca138d](https://github.com/snakemake/snakemake/commit/fca138d680ac0abad14df8dc208a84f1ae0dc0e2))

## [8.4.9](https://github.com/snakemake/snakemake/compare/v8.4.8...v8.4.9) (2024-02-15)


### Bug Fixes

* binary mem ([#2695](https://github.com/snakemake/snakemake/issues/2695)) ([18689b4](https://github.com/snakemake/snakemake/commit/18689b42dcc690deea7f27dc909379f603aaa1a9))
* Don't attempt to parse TBDString resources ([#2690](https://github.com/snakemake/snakemake/issues/2690)) ([f884b50](https://github.com/snakemake/snakemake/commit/f884b50eddf30f829b0f7d2bbda77d0148cf8cd7))


### Documentation

* added Mastodon follow label with just 'Follow' similar to X ([#2692](https://github.com/snakemake/snakemake/issues/2692)) ([7e36496](https://github.com/snakemake/snakemake/commit/7e364968fb601b2fd87dc5cd703e251b708a0e19))

## [8.4.8](https://github.com/snakemake/snakemake/compare/v8.4.7...v8.4.8) (2024-02-09)


### Bug Fixes

* fix bug causing FileNotFoundError when accessing checkpoint output. ([c81954d](https://github.com/snakemake/snakemake/commit/c81954d0bc10011a777aef7d0d5f089ee5748342))
* Fix collect-lookup attribute error ([#2687](https://github.com/snakemake/snakemake/issues/2687)) ([e39c74c](https://github.com/snakemake/snakemake/commit/e39c74c18d59f117f1e630a516ba00ee8a53c9e8))
* Fixed plot axis label on report ([#2683](https://github.com/snakemake/snakemake/issues/2683)) ([a4c2a03](https://github.com/snakemake/snakemake/commit/a4c2a03182d5f1bd5e7a9e389dc816cf7a9bc124))

## [8.4.7](https://github.com/snakemake/snakemake/compare/v8.4.6...v8.4.7) (2024-02-07)


### Documentation

* improve branch function docs ([e9d1a11](https://github.com/snakemake/snakemake/commit/e9d1a11ca67456d7304b580ae6e1f1ae01e712cd))

## [8.4.6](https://github.com/snakemake/snakemake/compare/v8.4.5...v8.4.6) (2024-02-06)


### Bug Fixes

* fix missing storage information when handling already completed checkpoints. This solves a bug causing failure to retrieve storage files in workflows with checkpoints. ([5791c60](https://github.com/snakemake/snakemake/commit/5791c6017195dd00a5d98c1da1e81fbfdde253b9))

## [8.4.5](https://github.com/snakemake/snakemake/compare/v8.4.4...v8.4.5) (2024-02-06)


### Bug Fixes

* for local execution, always unrestrictedly assume shared FS ([#2679](https://github.com/snakemake/snakemake/issues/2679)) ([0bee50b](https://github.com/snakemake/snakemake/commit/0bee50b24538e9e17f8679fe4e7022f335e32d72))
* support list of queries for storage provider ([#2674](https://github.com/snakemake/snakemake/issues/2674)) ([d53ef92](https://github.com/snakemake/snakemake/commit/d53ef92e79f5e8cf0187a1aa25c321c1478914d8))
* use default container image if nothing is provided ([#2677](https://github.com/snakemake/snakemake/issues/2677)) ([109c991](https://github.com/snakemake/snakemake/commit/109c991cea1ca757086488b3ca07f59cb843b46f))

## [8.4.4](https://github.com/snakemake/snakemake/compare/v8.4.3...v8.4.4) (2024-02-05)


### Bug Fixes

* fixed bug in handling of resource overrides for remote job submission ([5c06dd6](https://github.com/snakemake/snakemake/commit/5c06dd6a3273057144a8d5e1305c7cbdfc329b66))
* output of rulegraph, closes [#2656](https://github.com/snakemake/snakemake/issues/2656) ([#2671](https://github.com/snakemake/snakemake/issues/2671)) ([f9b9110](https://github.com/snakemake/snakemake/commit/f9b9110b2b0e25db25f76ceeabcf548b70aff839))
* SyntaxWarning due to invalid escape sequences in non-raw regex pattern string ([#2670](https://github.com/snakemake/snakemake/issues/2670)) ([3748d9d](https://github.com/snakemake/snakemake/commit/3748d9dc937195a54bc052994f6d2fa05373da42))

## [8.4.3](https://github.com/snakemake/snakemake/compare/v8.4.2...v8.4.3) (2024-02-02)


### Bug Fixes

* Do not scheduler execution message if no jobs are ready ([b1c4f47](https://github.com/snakemake/snakemake/commit/b1c4f47f894e344bd48b615f397c9ab03dc05f73))
* fix string resource definition in CLI and profile ([#2627](https://github.com/snakemake/snakemake/issues/2627)) ([bbd76ae](https://github.com/snakemake/snakemake/commit/bbd76ae53101a28f86386d12b3bdb10a89df2251))
* if report files are within storage, retrieve them from storage before loading into report ([60041bd](https://github.com/snakemake/snakemake/commit/60041bdbc3c6e941585f1876f3e5d598bfaa9caf))

## [8.4.2](https://github.com/snakemake/snakemake/compare/v8.4.1...v8.4.2) (2024-01-30)


### Bug Fixes

* allow lookup dpath or query to be a callable ([33f1637](https://github.com/snakemake/snakemake/commit/33f16379a5fb7f210769b888ab34fe8a8df96daf))
* fix error when passing callable as dpath or query of lookup function ([0e5b878](https://github.com/snakemake/snakemake/commit/0e5b8789c763201d00e5040f439436405d04b706))

## [8.4.1](https://github.com/snakemake/snakemake/compare/v8.4.0...v8.4.1) (2024-01-30)


### Fixes

* fixed resource handling in profiles

## [8.4.0](https://github.com/snakemake/snakemake/compare/v8.3.2...v8.4.0) (2024-01-29)


### Features

* add cols argument to lookup function; fix various minor bugs on cluster systems ([#2651](https://github.com/snakemake/snakemake/issues/2651)) ([ca7a602](https://github.com/snakemake/snakemake/commit/ca7a6022bacae77be45068adbf4386c2c93fb481))


### Bug Fixes

* batch bug [#2643](https://github.com/snakemake/snakemake/issues/2643) ([#2650](https://github.com/snakemake/snakemake/issues/2650)) ([2ecb21b](https://github.com/snakemake/snakemake/commit/2ecb21ba04088b9e6850447760f713784cf8b775))
* f-string in a more robust style? ([#2649](https://github.com/snakemake/snakemake/issues/2649)) ([2a50dc0](https://github.com/snakemake/snakemake/commit/2a50dc02bb709161d62d6f7dc5d6f2733e534c09))
* set ignore_incomplete to False in create_conda_envs  ([#2653](https://github.com/snakemake/snakemake/issues/2653)) ([4834a42](https://github.com/snakemake/snakemake/commit/4834a42180fb513670b310fcbaadd07a34adf0b7))
* Setting the value of ignore_incomplete Fixes [#2556](https://github.com/snakemake/snakemake/issues/2556) ([#2654](https://github.com/snakemake/snakemake/issues/2654)) ([05dac64](https://github.com/snakemake/snakemake/commit/05dac64666cca196e232beed4b6f1167dfbfc3bd))

## [8.3.2](https://github.com/snakemake/snakemake/compare/v8.3.1...v8.3.2) (2024-01-25)


### Bug Fixes

* do not require cores to be set for non-executing modes ([#2646](https://github.com/snakemake/snakemake/issues/2646)) ([30cf026](https://github.com/snakemake/snakemake/commit/30cf0261004be7a4bdce5ade563271e52c97ad6a))

## [8.3.1](https://github.com/snakemake/snakemake/compare/v8.3.0...v8.3.1) (2024-01-23)


### Documentation

* fix headings ([d947f85](https://github.com/snakemake/snakemake/commit/d947f85172f777ade0dc97bbe853c456517f28f6))

## [8.3.0](https://github.com/snakemake/snakemake/compare/v8.2.4...v8.3.0) (2024-01-23)


### Features

* implement semantic helper functions for input and param function handling ([#2344](https://github.com/snakemake/snakemake/issues/2344)) ([b4b5e51](https://github.com/snakemake/snakemake/commit/b4b5e51ee8601b81a7ac900b2d175603c5af90a9))
* support for continuously updated input (using Python queues) ([#2594](https://github.com/snakemake/snakemake/issues/2594)) ([db1c0ed](https://github.com/snakemake/snakemake/commit/db1c0edcadca499d51223d5cc72cbdfca1ff5d21))

## [8.2.4](https://github.com/snakemake/snakemake/compare/v8.2.3...v8.2.4) (2024-01-23)


### Bug Fixes

* fix exception when handling syntax error during parsing ([d5a7a56](https://github.com/snakemake/snakemake/commit/d5a7a564beea9850d1dac9f91429d8434a3aac46))

## [8.2.3](https://github.com/snakemake/snakemake/compare/v8.2.2...v8.2.3) (2024-01-19)


### Documentation

* handle overflow of content div ([b23e277](https://github.com/snakemake/snakemake/commit/b23e2776ce810d64225e7675039a33ac0e920bb9))

## [8.2.2](https://github.com/snakemake/snakemake/compare/v8.2.1...v8.2.2) (2024-01-19)


### Documentation

* add missing doc dependency ([7ba9c21](https://github.com/snakemake/snakemake/commit/7ba9c21ac949ae6952044944af68874c516791ff))
* fix typo in `rules.py` ([#2636](https://github.com/snakemake/snakemake/issues/2636)) ([8bc2919](https://github.com/snakemake/snakemake/commit/8bc291933201e205717f3346b743153b57764a6a))
* use sphinxawesome-theme instead of lutra ([4401e9c](https://github.com/snakemake/snakemake/commit/4401e9c2c53adb46609930b872b5da27b5af58fe))

## [8.2.1](https://github.com/snakemake/snakemake/compare/v8.2.0...v8.2.1) (2024-01-17)


### Bug Fixes

* do not require cores to be set for rule-level methods of the workflow API or the corresponding CLI commands (e.g. --lint). ([#2629](https://github.com/snakemake/snakemake/issues/2629)) ([2040468](https://github.com/snakemake/snakemake/commit/20404688e91beadd8790e3c6cdcb727bc47d597e))
* fix false complaints about rules with multiple output files ([#2628](https://github.com/snakemake/snakemake/issues/2628)) ([b1b4f5b](https://github.com/snakemake/snakemake/commit/b1b4f5b0adcd8066b7a6376e9b56778014f9921b))
* migration guide typo and wrong link ([#2625](https://github.com/snakemake/snakemake/issues/2625)) ([645f3d1](https://github.com/snakemake/snakemake/commit/645f3d1426c1dcd41adf813bdb833060775ffda5))

## [8.2.0](https://github.com/snakemake/snakemake/compare/v8.1.3...v8.2.0) (2024-01-16)


### Features

* add method to obtain group args for spawned jobs ([bd1b450](https://github.com/snakemake/snakemake/commit/bd1b450950c5b6e510778ca09fe11beffd9a6ff3))


### Bug Fixes

* properly resolve wildcards in group components ([#2620](https://github.com/snakemake/snakemake/issues/2620)) ([c788a46](https://github.com/snakemake/snakemake/commit/c788a465f25dc700103df6475d566d81e2c88bef))
* return set of rules when obtaining allowed rules for remote job ([2c44cf6](https://github.com/snakemake/snakemake/commit/2c44cf641335812ca0877e3987726403c8919415))

## [8.1.3](https://github.com/snakemake/snakemake/compare/v8.1.2...v8.1.3) (2024-01-15)


### Bug Fixes

* bug with preemptible rules ([#2616](https://github.com/snakemake/snakemake/issues/2616)) ([c6d7141](https://github.com/snakemake/snakemake/commit/c6d71410a042f2c0075a75fc024bd1e05fe0af15))
* do not pass snakefile as metadata when wms monitor flag is used ([#2573](https://github.com/snakemake/snakemake/issues/2573)) ([13b3205](https://github.com/snakemake/snakemake/commit/13b3205a6e808dc48fd7a3be9652b9f0c8887648))
* use default group settings if not execution workflow (fixes attribute error occurring with --report) ([#2617](https://github.com/snakemake/snakemake/issues/2617)) ([21e9964](https://github.com/snakemake/snakemake/commit/21e996421abd72122a1403a10eb74c111bde28be))

## [8.1.2](https://github.com/snakemake/snakemake/compare/v8.1.1...v8.1.2) (2024-01-12)


### Bug Fixes

* local mtime handling in case of storage plugins and cleaner error message for parallel storage retrieval ([#2611](https://github.com/snakemake/snakemake/issues/2611)) ([880b264](https://github.com/snakemake/snakemake/commit/880b2645c961a98c0eb33b9cdab75d7804860c0a))
* Migrate away from deprecated pulp API ([#2610](https://github.com/snakemake/snakemake/issues/2610)) ([fb26640](https://github.com/snakemake/snakemake/commit/fb2664008d6cffd09c52ec9fd9c994bd198ed69c))

## [8.1.1](https://github.com/snakemake/snakemake/compare/v8.1.0...v8.1.1) (2024-01-11)


### Bug Fixes

* deduplicate input files before retrieval from storage ([#2600](https://github.com/snakemake/snakemake/issues/2600)) ([37cf475](https://github.com/snakemake/snakemake/commit/37cf475a4205d306a3969ab0e3249f0c9f7d4e19))

## [8.1.0](https://github.com/snakemake/snakemake/compare/v8.0.1...v8.1.0) (2024-01-08)


### Features

* add --sdm short opt for --deployment ([#2551](https://github.com/snakemake/snakemake/issues/2551)) ([fd8b4b0](https://github.com/snakemake/snakemake/commit/fd8b4b08da778097071ea8c32c95e21e610c6614))


### Bug Fixes

* add mamba to docker image ([eb0c884](https://github.com/snakemake/snakemake/commit/eb0c88495b506300fdf7c1afc4c02d6b91c6a582))
* correctly report lineno ([#2584](https://github.com/snakemake/snakemake/issues/2584)) ([967a0d7](https://github.com/snakemake/snakemake/commit/967a0d7cb90c63c7be00b49fd535cf1029f63b5b))
* move apptainer into separate env in docker image ([94e9e2c](https://github.com/snakemake/snakemake/commit/94e9e2c33ee63551aa3630c7106344ee2fa11f4d))
* single line f-string format error in py3.12 ([#2588](https://github.com/snakemake/snakemake/issues/2588)) ([87c06c0](https://github.com/snakemake/snakemake/commit/87c06c0f5745f577c12db39852c6f763a2d41954))


### Documentation

* add note on google executor backends ([ff8683c](https://github.com/snakemake/snakemake/commit/ff8683c80511a78103c14442a3d3fb81bcbef2cc))
* diff 7 and 8 ([#2561](https://github.com/snakemake/snakemake/issues/2561)) ([ba22e07](https://github.com/snakemake/snakemake/commit/ba22e07af2ab8248888cb57f4bcb8a9eb2623977))

## [8.0.1](https://github.com/snakemake/snakemake/compare/v8.0.0...v8.0.1) (2023-12-21)


### Bug Fixes

* remove bash completion entrypoint (no longer supported, was too slow to be usable anyway) ([922b53a](https://github.com/snakemake/snakemake/commit/922b53aa0cba05da067cc67fccc6852bbc161edb))


### Documentation

* fix cli options rendering ([264c1a9](https://github.com/snakemake/snakemake/commit/264c1a92e824323b3060192d7a22aeb0c07678e0))
* fixes in migration guide ([f8adefa](https://github.com/snakemake/snakemake/commit/f8adefaa605bd87df40440a211276a6825d138ed))

## [8.0.0](https://github.com/snakemake/snakemake/compare/v7.32.2...v8.0.0) (2023-12-20)



### ⚠ BREAKING CHANGES

Snakemake 8 marks the beginning of decomposing Snakemake into a framework of plugins. This enables the democratization of method development within the Snakemake ecosystem.
We start with plugins for storage and execution backends. In the future, there will be plugins for the scheduling, metadata, software deployment, reporting, and many more.
This way, it will be possible to easily launch and explore new developments in workflow management and reproducible data analysis without the need to get your work merged into the main codebase of Snakemake and also without the need to develop a new workflow management system as a proof of concept.

In detail, Snakemake 8 introduces the following changes. Unfortunately it was unavoidable to break some usages (we apologize).
Nevertheless, we tried to ensure that every removed or modified feature has been replaced with an equivalent reimplementation, as outlined in our [migration docs](https://snakemake.readthedocs.io/en/latest/getting_started/migration.html#migrating-to-snakemake-8).
While Snakemake 8 has an even more thorough testing framework than any release before, and while it has been quite heavily tested in practice by us, you might initially experience bugs and glitches for which we want to apologize beforehand.
We think that the massive codebase improvements are worth it in the long run, and hope that everything goes well.
As always, any pull requests with test cases and pointers to bugs are more than welcome.

#### Detailed breaking changes

* removed the long time ago deprecated support for dynamic, version, and subworkflow (see [the migration docs](https://snakemake.readthedocs.io/en/latest/getting_started/migration.html#migrating-to-snakemake-8))
* migrated old remote providers into storage plugins (see [the migration docs](https://snakemake.readthedocs.io/en/latest/getting_started/migration.html#migrating-to-snakemake-8))
* migrated execution backends into plugins, including a change in the respective command line interfaces (see [the migration docs](https://snakemake.readthedocs.io/en/latest/getting_started/migration.html#migrating-to-snakemake-8))
* deprecates `--use-conda` and `--use-singularity` in favor of `--software-deployment-method conda` or `--software-deployment-method apptainer` and `--software-deployment-method conda apptainer` (see [the migration docs](https://snakemake.readthedocs.io/en/latest/getting_started/migration.html#migrating-to-snakemake-8))
* profile support is now versioned, such that different profiles can be written for different minimum Snakemake versions (see [the migration docs](https://snakemake.readthedocs.io/en/latest/getting_started/migration.html#migrating-to-snakemake-8))
* redesigned Snakemake API. It now uses a modern, dataclass based approach (see [the migration docs](https://snakemake.readthedocs.io/en/latest/getting_started/migration.html#migrating-to-snakemake-8))

### Features

* add ability to inject conda environments into running Snakefile ([#2479](https://github.com/snakemake/snakemake/issues/2479)) ([6140e29](https://github.com/snakemake/snakemake/commit/6140e29864cf0fcdd83194f218408867c54b730d))
* add functionality for deploying sources if no shared FS is assumed ([#2486](https://github.com/snakemake/snakemake/issues/2486)) ([76eac3c](https://github.com/snakemake/snakemake/commit/76eac3c570624f818ff43f4759a27a3e284bf03c))
* add option to control software deployment mode (shared or non shared FS) ([#2525](https://github.com/snakemake/snakemake/issues/2525)) ([04ec2c0](https://github.com/snakemake/snakemake/commit/04ec2c0262b2cb96cbcd7edbbb2596979c1703ae))
* allow detailed configuration of shared FS usage ([#2528](https://github.com/snakemake/snakemake/issues/2528)) ([0d34be9](https://github.com/snakemake/snakemake/commit/0d34be90040f937021eb25becb4ef5d5aae66473))
* allow environment variables in string values of profile (e.g. paths may now contain elements like $USER). ([58dc70c](https://github.com/snakemake/snakemake/commit/58dc70c513d08cb217647ec591df93b40517b650))
* allow python expressions in --set-resources ([#2521](https://github.com/snakemake/snakemake/issues/2521)) ([022a31e](https://github.com/snakemake/snakemake/commit/022a31e0e75416904589c62b58e54c952d930c69))
* allow to set latency_wait in executor test suite ([c0bca0b](https://github.com/snakemake/snakemake/commit/c0bca0bb8bf66fb34def7459f757140ce3825a25))
* automatically upload workflow sources to default storage provider if no shared FS is used ([a450c49](https://github.com/snakemake/snakemake/commit/a450c4998de3ab7deb0fb2bc19dc59fdc484309d))
* Faster ci test setup  ([#2489](https://github.com/snakemake/snakemake/issues/2489)) ([4798e8a](https://github.com/snakemake/snakemake/commit/4798e8ac226bded585e9fe31d43ae9e93a598780))
* implement precommand ([#2482](https://github.com/snakemake/snakemake/issues/2482)) ([ff0f979](https://github.com/snakemake/snakemake/commit/ff0f979b68b1e12be8151c7f5547c6a13ad3ee9a))
* redesigned Snakemake API. It now uses a modern, dataclass based approach ([#2403](https://github.com/snakemake/snakemake/issues/2403)) ([2be3bfa](https://github.com/snakemake/snakemake/commit/2be3bfa4841967928069a2a024554b8a86b699f1))
* support for external executor plugins ([#2305](https://github.com/snakemake/snakemake/issues/2305)) ([c9eaa4e](https://github.com/snakemake/snakemake/commit/c9eaa4e12e4a348f93e5ea5793faaec1fd547fac))
* version specific profile config files (profile/config.v8+.yaml with profile/config.yaml as fallback that matches any version) ([#2498](https://github.com/snakemake/snakemake/issues/2498)) ([47e5811](https://github.com/snakemake/snakemake/commit/47e581181f952884577f0237a1aa9457ee9554dd))

### Bug Fixes

* adapt to changes in snakemake-interface-executor-plugins ([635c68a](https://github.com/snakemake/snakemake/commit/635c68abe3c6e01a70803f2423718a99ea056a00))
* add storage provider args to deploy sources command ([67178e3](https://github.com/snakemake/snakemake/commit/67178e31e68dc165c503d49eaf40340a9ad65e90))
* add testcase for script directive to work with Python 3.7 and corresponding fix. ([0b4ae2e](https://github.com/snakemake/snakemake/commit/0b4ae2e1155378a78e0786ad4dfd5d44fa41c9ae))
* allow pepfile and pepschema to take pathlib ([#2546](https://github.com/snakemake/snakemake/issues/2546)) ([ca91661](https://github.com/snakemake/snakemake/commit/ca91661bf2ff215f005e8d9351fa3320d5cf5498))
* also inherit rule proxies if there is no rulename modifier specified in a use rule statement ([#2440](https://github.com/snakemake/snakemake/issues/2440)) ([1570289](https://github.com/snakemake/snakemake/commit/15702891b31635a79f31e857d09f3e285f71717b))
* assume at most 8GB memory for default resources. This way, we avoid exploding memory requirements for large input files that are very unlikely to be put entirely into memory by any tool. ([11c2ecc](https://github.com/snakemake/snakemake/commit/11c2eccfa72a2732d20a66bbf57995ee9dfd14ae))
* comparison to float in scheduler ([ef44d84](https://github.com/snakemake/snakemake/commit/ef44d844b81d671604a1ce285a43ca1d5ea59d96))
* detect job paths that leave and then enter a group. Such paths are invalid because then the group depends on itself. ([#2527](https://github.com/snakemake/snakemake/issues/2527)) ([5383a4d](https://github.com/snakemake/snakemake/commit/5383a4d85f8b7a3fb61fd0b457fcd1d008c2255f))
* ensure that auto deployment of default storage provider works in containers with read-only root home. ([1a347ff](https://github.com/snakemake/snakemake/commit/1a347ffd8b2144aee5c35ea17fec2262c5cc9c40))
* ensure that log and benchmark files are uploaded to storage as well ([#2545](https://github.com/snakemake/snakemake/issues/2545)) ([6aabb5d](https://github.com/snakemake/snakemake/commit/6aabb5db16634c077ba808dbd000a3ed67d6a3c0))
* ensure that targetjob is always forced. This fixes a bug causing run-directive rules to not being executed even when enforced via e.g. -R. ([#2448](https://github.com/snakemake/snakemake/issues/2448)) ([b2a60d5](https://github.com/snakemake/snakemake/commit/b2a60d5674c84f58c6e48af5d675ce4690a1d625))
* fix cache handling and unlock handling ([2f4d5e1](https://github.com/snakemake/snakemake/commit/2f4d5e11b36f199d3335839e9b5dc5d2094d11f8))
* fix nargs definition for --deploy-sources ([fc252c8](https://github.com/snakemake/snakemake/commit/fc252c80227e75a4fcf869f828d3c7d5d066f794))
* fix path handling when detective profiles ([fe63881](https://github.com/snakemake/snakemake/commit/fe63881371036fc6d55e6a113cb3ed5029b31d2d))
* fix storage handling on windows by converting all paths to posix paths ([#2519](https://github.com/snakemake/snakemake/issues/2519)) ([7864a76](https://github.com/snakemake/snakemake/commit/7864a76c41ad7839fb3b9aeb1c4468215ba0cb21))
* handle different f-string tokens in py3.12 ([#2485](https://github.com/snakemake/snakemake/issues/2485)) ([f2c7613](https://github.com/snakemake/snakemake/commit/f2c761320a5a73d6027ae3649843e6bf6a24f324))
* handle storage for local jobs; add test case ([6d978ef](https://github.com/snakemake/snakemake/commit/6d978ef31060393bd3a69b7973170e3c79473705))
* handling of group jobs when obtaining temp input files ([71be1de](https://github.com/snakemake/snakemake/commit/71be1de734032ae8c9e8b07fad57e36b321be421))
* import ([#2402](https://github.com/snakemake/snakemake/issues/2402)) ([2c831f1](https://github.com/snakemake/snakemake/commit/2c831f1fa98813cf5f69ecb046aad1364f514238))
* improved error handling for storage upload; fixed bugs caused by outdated calls to IOFile.exists(). ([720bb84](https://github.com/snakemake/snakemake/commit/720bb8400ecf8e6c3caf3ca6b47d7dda4f1f7ba9))
* improved error messages in case of invalid storage queries ([9671fd0](https://github.com/snakemake/snakemake/commit/9671fd0770fbe2d79d3546ca5f19c8fc39ffc25a))
* in addition to localrules statement,  infer that job is local if it has any input or output file that is marked as local ([#2541](https://github.com/snakemake/snakemake/issues/2541)) ([e8b682b](https://github.com/snakemake/snakemake/commit/e8b682be0bd52a01ec3ebdbbcb1ec018950819b2))
* only deactivate conda inject envs upon workflow tear down ([#2503](https://github.com/snakemake/snakemake/issues/2503)) ([e6dfdd4](https://github.com/snakemake/snakemake/commit/e6dfdd49cbaad228827f62a48cfe6419e4c1715e))
* Panoptes --wms-monitor-arg ([#2444](https://github.com/snakemake/snakemake/issues/2444)) ([98d2bdf](https://github.com/snakemake/snakemake/commit/98d2bdfd1e3aeddf189b272d3e4042632248c10f))
* proper reuse of rule proxies when importing several times from the same module ([#2404](https://github.com/snakemake/snakemake/issues/2404)) ([e867dda](https://github.com/snakemake/snakemake/commit/e867dda24dff306f42939ad0d4d93d32ec94f6e5))
* Restore backward compatibility for Google Life Sciences executor ([#2461](https://github.com/snakemake/snakemake/issues/2461)) ([5e3a464](https://github.com/snakemake/snakemake/commit/5e3a46476d78d5d52340a9ffa327d18a5e7e9828))
* shadow "full" mode ignore symlinks ([#2516](https://github.com/snakemake/snakemake/issues/2516)) ([1d58120](https://github.com/snakemake/snakemake/commit/1d5812027e5d1a8d5f2e45175fa38d28ca32763e))
* show failed logs in executor testcases ([92f7bf4](https://github.com/snakemake/snakemake/commit/92f7bf4b86f491073898b2385c3ed45a557a3b4d))
* Slack log service ([#2537](https://github.com/snakemake/snakemake/issues/2537)) ([26eb4ba](https://github.com/snakemake/snakemake/commit/26eb4babc28464a24c33d75703e7bf4f15c0e33f))
* sort report (sub-)categories in lexicographical order ([#2449](https://github.com/snakemake/snakemake/issues/2449)) ([d0705ad](https://github.com/snakemake/snakemake/commit/d0705adb4702ead9db0c26809859c7b769d800e1))
* update minimum snakemake-interface-storage-plugins version ([0ef7226](https://github.com/snakemake/snakemake/commit/0ef72262e29a5b22cdf016a4ab6aba4b8dbc686d))
* use temporary directory (faster, more likely local, always writable) for persistence and source cache in case of remote execution without shared fs ([#2502](https://github.com/snakemake/snakemake/issues/2502)) ([c8fa7ba](https://github.com/snakemake/snakemake/commit/c8fa7ba3ee70c9c62011a3839758a5eb8fde16f8))
* wait for logs before showing them on error ([a4ff328](https://github.com/snakemake/snakemake/commit/a4ff3280db0beb4f1a077ee880433f767c4ad142))


### Documentation

* document name directive with example ([#2534](https://github.com/snakemake/snakemake/issues/2534)) ([cce5551](https://github.com/snakemake/snakemake/commit/cce555142814f5bd1d73e68b9a17b772454817d4))
* fix syntax in cluster example ([#2460](https://github.com/snakemake/snakemake/issues/2460)) ([64e9645](https://github.com/snakemake/snakemake/commit/64e964554748bdee93bad1c7e6cd2924595c414f))
* notes on arm based machines in tutorial docs ([0586f04](https://github.com/snakemake/snakemake/commit/0586f04d2b03443a25deddd017f303156acdcd9c))
* **rust:** Fix typo on rust-script version ([#2488](https://github.com/snakemake/snakemake/issues/2488)) ([a79dd94](https://github.com/snakemake/snakemake/commit/a79dd94f330f1d5c1fba7d16cb5b08bd780d950d))

## [7.32.4](https://github.com/snakemake/snakemake/compare/v7.32.3...v7.32.4) (2023-08-18)


### Bug Fixes

* always sort report (sub-)categories in lexicographical order
* also inherit rule proxies if there is no rulename modifier specified in a use rule statement
* ensure that targetjob is always forced. This fixes a bug causing run-directive rules to not being executed even when enforced via e.g. -R.


## [7.32.3](https://github.com/snakemake/snakemake/compare/v7.32.2...v7.32.3) (2023-08-07)


### Bug Fixes

* fix bug occurring when using multiple `use rule` statements in combination with the rules object for referring to output of already defined rules.


## [7.32.2](https://github.com/snakemake/snakemake/compare/v7.32.1...v7.32.2) (2023-08-07)


### Bug Fixes

* unnecessary set Snakefile in AzBatch executor ([#2397](https://github.com/snakemake/snakemake/issues/2397)) ([78e6d6e](https://github.com/snakemake/snakemake/commit/78e6d6ec6a1ad930e40d6edfe9f7210232a674f2))

## [7.32.1](https://github.com/snakemake/snakemake/compare/v7.32.0...v7.32.1) (2023-08-05)


### Bug Fixes

* add missing spaces between lines that get concatenated. ([#2268](https://github.com/snakemake/snakemake/issues/2268)) ([7238458](https://github.com/snakemake/snakemake/commit/7238458c6d56c5d94787b93668718358ad44e9ef))
* better message about profile usage upon execution ([#2391](https://github.com/snakemake/snakemake/issues/2391)) ([cf8aea5](https://github.com/snakemake/snakemake/commit/cf8aea5862767f104c7e03c09369d401f25d50e7))
* do not overwrite default resources setting in azure batch executor ([#2395](https://github.com/snakemake/snakemake/issues/2395)) ([4aef3b9](https://github.com/snakemake/snakemake/commit/4aef3b93ddf029c532fabe38e05c262dd4237b5f))
* updating of non-dict config values gives error ([#2364](https://github.com/snakemake/snakemake/issues/2364)) ([b33aeec](https://github.com/snakemake/snakemake/commit/b33aeecdf51afba6012007a2b125b9c87b7b98f2))
* wrong rule names when nesting module imports ([#1817](https://github.com/snakemake/snakemake/issues/1817)) ([65c79a4](https://github.com/snakemake/snakemake/commit/65c79a48f956077839bb5ab1ea8d60a5f0ddecab))


### Documentation

* basics.rst: suggest VS Code instead of deprecated Atom as IDE ([#2368](https://github.com/snakemake/snakemake/issues/2368)) ([1357316](https://github.com/snakemake/snakemake/commit/135731605b974915cd2bd78b88a981f974bc7b78))

## [7.32.0](https://github.com/snakemake/snakemake/compare/v7.31.1...v7.32.0) (2023-08-03)


### Features

* add support for Kubernetes service account name spec ([#2254](https://github.com/snakemake/snakemake/issues/2254)) ([3370426](https://github.com/snakemake/snakemake/commit/3370426da7ee78af5de54689f623e2b5afa45f1f))


### Bug Fixes

* Enable values with an = sign in default_resources ([#2340](https://github.com/snakemake/snakemake/issues/2340)) ([c1c9229](https://github.com/snakemake/snakemake/commit/c1c922904f09c133e39872346a541e7cd216d0d2))
* Escape workdir paths for potential spaces in paths ([#2196](https://github.com/snakemake/snakemake/issues/2196)) ([9261f7e](https://github.com/snakemake/snakemake/commit/9261f7ea50a8ae424c015faca73b7811fb51d093))
* ga4gh executor resources ([#2042](https://github.com/snakemake/snakemake/issues/2042)) ([ad6eaef](https://github.com/snakemake/snakemake/commit/ad6eaef6bac05d4de682f59d9d4a088f143b5798))
* print exceptions when job is not a shell job ([#2385](https://github.com/snakemake/snakemake/issues/2385)) ([8a37b85](https://github.com/snakemake/snakemake/commit/8a37b8584f216ada10caffcbb8b731efd675376a))
* remote-azblob-sasToken-Authorization ([#1800](https://github.com/snakemake/snakemake/issues/1800)) ([bc854a7](https://github.com/snakemake/snakemake/commit/bc854a7e012cac751b708df83378fd5791e6e6fc))
* wms-monitor now gets data in correct json format ([#2347](https://github.com/snakemake/snakemake/issues/2347)) ([7fafa7a](https://github.com/snakemake/snakemake/commit/7fafa7ace72f8a727457f4abe6db2f9ed2d74d64))


### Documentation

* fix a copy&paste (?) mistake ([#2386](https://github.com/snakemake/snakemake/issues/2386)) ([d878847](https://github.com/snakemake/snakemake/commit/d87884749fd9450062f6fde5b7727867396e7a78))

## [7.31.1](https://github.com/snakemake/snakemake/compare/v7.31.0...v7.31.1) (2023-08-02)


### Bug Fixes

* require python &gt;=3.7 again (the python 3.9 dependency was unnecessary) ([#2372](https://github.com/snakemake/snakemake/issues/2372)) ([0d0e9c4](https://github.com/snakemake/snakemake/commit/0d0e9c4cf48a97952464e6da476ed7661d629ce3))


### Documentation

* update CHANGELOG.md: add minimum Python version bump ([#2370](https://github.com/snakemake/snakemake/issues/2370)) ([48e934d](https://github.com/snakemake/snakemake/commit/48e934dcf96e4e8fd30c81cab3674583bf049a45))

## [7.31.0](https://github.com/snakemake/snakemake/compare/v7.30.2...v7.31.0) (2023-07-26)


### Features

* Add support for Google Service Accounts and GCE VM network configuration ([#2318](https://github.com/snakemake/snakemake/issues/2318)) ([2b754aa](https://github.com/snakemake/snakemake/commit/2b754aae535ef76bd2dd34bc31d5c9f5c69363de))

## [7.30.2](https://github.com/snakemake/snakemake/compare/v7.30.1...v7.30.2) (2023-07-20)

### Breaking changes

* Bump minimum Python version from 3.7 to 3.9 ([#2369](https://github.com/snakemake/snakemake/issues/2369)) ([4608163](https://github.com/snakemake/snakemake/pull/2341/commits/4608163727bb32e216f1a26adc61d4c15d4b6a47))

### Bug Fixes

* do not allow setting benchmark and between-workflow caching for the same rule. The reason is that when the result is taken from cache, there is no way to fill the benchmark file with any reasonable values. ([#2335](https://github.com/snakemake/snakemake/issues/2335)) ([e2d64fa](https://github.com/snakemake/snakemake/commit/e2d64fad76b8ca1805eeaa48c0bf8d1fb7bf4736))
* ensure lazy evaluation of resource functions/callables (this also entails, for now, a removal of the thread statistics in the yellow job stats table); further, added some clarifying sentences about resource function evaluation to the docs ([#2356](https://github.com/snakemake/snakemake/issues/2356)) ([4c591b7](https://github.com/snakemake/snakemake/commit/4c591b72b31d6c6c36b43f1d7773d8317352fbc9))
* handle non-PEP440 versions of apptainer/singulariy ([#2337](https://github.com/snakemake/snakemake/issues/2337)) ([dea6ba8](https://github.com/snakemake/snakemake/commit/dea6ba8808793b88c7553880bde48711abb037f8))
* remote GS builds too many inventories; io:collect_mtime always uses uncached mtime ([#2266](https://github.com/snakemake/snakemake/issues/2266)) ([bad9115](https://github.com/snakemake/snakemake/commit/bad91152eeb70693e1459324f738a8c481378801))
* Solve apptainer version issue ([#2333](https://github.com/snakemake/snakemake/issues/2333)) ([a876e0f](https://github.com/snakemake/snakemake/commit/a876e0f5e187168eb269b504918c6aeff1496f16))
* SyntaxWarnings due to non-raw regex pattern strings ([#2359](https://github.com/snakemake/snakemake/issues/2359)) ([a08c0b0](https://github.com/snakemake/snakemake/commit/a08c0b071b2f9a9212117bbcf560fa67f1a02178))


### Documentation

* clarify minimum Snakemake version for profiles ([86dc277](https://github.com/snakemake/snakemake/commit/86dc277d530a557c9bdd6784b863f63ab859a1c7))
* clarify the channel priority in environment definition deployment.rst ([#2352](https://github.com/snakemake/snakemake/issues/2352)) ([76aa964](https://github.com/snakemake/snakemake/commit/76aa964c38b4aa069d9cce6f8f43c91c7d496cfb))
* fix typo (stackoverflow issue) ([#2365](https://github.com/snakemake/snakemake/issues/2365)) ([f770984](https://github.com/snakemake/snakemake/commit/f7709844cd932465859a2095edafcf9baa8c2bf7))
* note on using checkpoint mechanism only for input function, not for params or resources. ([#2353](https://github.com/snakemake/snakemake/issues/2353)) ([4be2f9d](https://github.com/snakemake/snakemake/commit/4be2f9dd9fb41dc169bae068753ceed9552248e7))

## [7.30.1](https://github.com/snakemake/snakemake/compare/v7.30.0...v7.30.1) (2023-06-28)


### Bug Fixes

* conda env inside script ([#1812](https://github.com/snakemake/snakemake/issues/1812)) ([49cac6a](https://github.com/snakemake/snakemake/commit/49cac6ac67fba360f2f35be7ab1972b2d8cc1f8b))

## [7.30.0](https://github.com/snakemake/snakemake/compare/v7.29.0...v7.30.0) (2023-06-28)


### Features

* allow profiles to be YTE templates; adapt to eido 2.0 ([#2325](https://github.com/snakemake/snakemake/issues/2325)) ([67d9ff2](https://github.com/snakemake/snakemake/commit/67d9ff20ea61186d3818e7bc1d33e4414058fc1f))

## [7.29.0](https://github.com/snakemake/snakemake/compare/v7.28.3...v7.29.0) (2023-06-21)


### Features

* introduce --workflow-profile for additional workflow specific profiles that overwrite global profiles; add ability to define key-value CLI flags like --set-threads or --set-resources as multi-level dictionaries in profile config yaml files ([#2310](https://github.com/snakemake/snakemake/issues/2310)) ([9675c17](https://github.com/snakemake/snakemake/commit/9675c17d4d7cbb95e589767974faa9219dd4154d))


### Bug Fixes

* addressing [#2197](https://github.com/snakemake/snakemake/issues/2197) by allowing 256 character account names in slurm ([#2198](https://github.com/snakemake/snakemake/issues/2198)) ([ab58c65](https://github.com/snakemake/snakemake/commit/ab58c652847c03a9f1529d2d7632f2788a5fadc4))
* removed distutils from snakemake ([#2312](https://github.com/snakemake/snakemake/issues/2312)) ([9b8c362](https://github.com/snakemake/snakemake/commit/9b8c3620e8c14e322ba15b7d044b9deab1854b2a))
* Update __init__.py to move "file" param to "print" ([#2291](https://github.com/snakemake/snakemake/issues/2291)) ([92352b6](https://github.com/snakemake/snakemake/commit/92352b69d14ef196b0253561c78fa04ffa25d73e))

## [7.28.3](https://github.com/snakemake/snakemake/compare/v7.28.2...v7.28.3) (2023-06-16)


### Bug Fixes

* Detect pandas availability to select serializer ([#2300](https://github.com/snakemake/snakemake/issues/2300)) ([e08a771](https://github.com/snakemake/snakemake/commit/e08a771f90aef84f3075e07c8d4e4c0f7881047c))


### Performance Improvements

* avoid superfluous mtime checks when the same file is referred to by multiple jobs ([#2284](https://github.com/snakemake/snakemake/issues/2284)) ([eb6e2e1](https://github.com/snakemake/snakemake/commit/eb6e2e161f01c61b139d95bcf1ddfa862f8029ba))


### Documentation

* update docs for azbatch and dockerhub ref ([#2298](https://github.com/snakemake/snakemake/issues/2298)) ([908dbf1](https://github.com/snakemake/snakemake/commit/908dbf143d4b1625fa6ee80f2b4eb713a6411a3f))

## [7.28.2](https://github.com/snakemake/snakemake/compare/v7.28.1...v7.28.2) (2023-06-13)


### Bug Fixes

* fix pandas import handling in metadata persistence ([27f7b40](https://github.com/snakemake/snakemake/commit/27f7b4014eaea66aa4e599aa854dda75822d30a0))

## [7.28.1](https://github.com/snakemake/snakemake/compare/v7.28.0...v7.28.1) (2023-06-11)


### Bug Fixes

* Bump yte from &gt;=1.0,&lt;2.0 to >=1.5.1,<2.0 ([#2275](https://github.com/snakemake/snakemake/issues/2275)) ([8c0b34f](https://github.com/snakemake/snakemake/commit/8c0b34f869e4f65ff2e47cf5f1e2863bd104f8e7))
* remove superfluous dependency ([aad61a0](https://github.com/snakemake/snakemake/commit/aad61a0131d7ca0f7393af23b98b1db702cd976d))

## [7.28.0](https://github.com/snakemake/snakemake/compare/v7.27.0...v7.28.0) (2023-06-11)


### Features

* Added native support for execution via Azure Batch ([#1953](https://github.com/snakemake/snakemake/issues/1953)) ([#2246](https://github.com/snakemake/snakemake/issues/2246)) ([0f9c49f](https://github.com/snakemake/snakemake/commit/0f9c49fe8643cca0e42e3b091cf9706a7feb877d))

## [7.27.0](https://github.com/snakemake/snakemake/compare/v7.26.0...v7.27.0) (2023-06-07)


### Features

* Allow the environment variable SNAKEMAKE_CONDA_PREFIX to be present without --use-conda ([#2263](https://github.com/snakemake/snakemake/issues/2263)) ([e4eba8d](https://github.com/snakemake/snakemake/commit/e4eba8d72b84aaa460c7d1b1ac54b607e844d782))


### Bug Fixes

* adapt linting rule to Python 3.11 ([a3a5c58](https://github.com/snakemake/snakemake/commit/a3a5c58cbbbe9a84b7383ce046b5981271288979))

## [7.26.0](https://github.com/snakemake/snakemake/compare/v7.25.4...v7.26.0) (2023-05-22)


### Features

* allow config files to be processed with YTE ([#2269](https://github.com/snakemake/snakemake/issues/2269)) ([8e1c22f](https://github.com/snakemake/snakemake/commit/8e1c22ff54e85ee941c6e0ac74dd594fce80efbb))

## [7.25.4](https://github.com/snakemake/snakemake/compare/v7.25.3...v7.25.4) (2023-05-12)


### Bug Fixes

* fix scrolling behavior in landing page of report for large workflows ([63c0c31](https://github.com/snakemake/snakemake/commit/63c0c31c222a921a843d83e330a6f91e430f209a))
* report spacing ([f3954b3](https://github.com/snakemake/snakemake/commit/f3954b33536314e6a252c044eee5f424dd234065))


### Documentation

* fix statement about logging ([#2252](https://github.com/snakemake/snakemake/issues/2252)) ([56c24b6](https://github.com/snakemake/snakemake/commit/56c24b6436cee9a4962d006bb201708c7a37c474))

## [7.25.3](https://github.com/snakemake/snakemake/compare/v7.25.2...v7.25.3) (2023-05-03)


### Bug Fixes

* fix missed wildcard constraints when using local rule inheritance ([#2242](https://github.com/snakemake/snakemake/issues/2242)) ([8e94785](https://github.com/snakemake/snakemake/commit/8e947858dd510cb4c813f24093b5be843fa4cf6c))

## [7.25.2](https://github.com/snakemake/snakemake/compare/v7.25.1...v7.25.2) (2023-04-28)


### Bug Fixes

* Fix inconsistencies between detailed summary and normal summary ([#2218](https://github.com/snakemake/snakemake/issues/2218)) ([d903123](https://github.com/snakemake/snakemake/commit/d9031236563b7dd8e31ed27208c9ad39699f765e))
* Fix race condition when creating lock directory ([#2225](https://github.com/snakemake/snakemake/issues/2225)) ([66ea4d1](https://github.com/snakemake/snakemake/commit/66ea4d199e3d9266b1b5fdb8752772e8137ffdea))
* quote paths given to singularity in order to ensure that it does not fail when paths contain whitespace ([#2190](https://github.com/snakemake/snakemake/issues/2190)) ([a572fb7](https://github.com/snakemake/snakemake/commit/a572fb7b8f00e39723cd98d6936f63171b26c8d9))


### Documentation

* added changelog info for &gt;v7.19.1 parsing error of "hh:mm:ss" time format in runtime resource ([#2189](https://github.com/snakemake/snakemake/issues/2189)) ([2889f38](https://github.com/snakemake/snakemake/commit/2889f3851b64cee7885fbf73f64a453eed5e806a))
* update misc/vim/Readme with info for packer.nvim ([#2095](https://github.com/snakemake/snakemake/issues/2095)) ([32166a7](https://github.com/snakemake/snakemake/commit/32166a7fce95312bfa6b6d3ae76bf94accf6d5de))
* Update workflow syntax with priority directive ([#2188](https://github.com/snakemake/snakemake/issues/2188)) ([af10db5](https://github.com/snakemake/snakemake/commit/af10db56b11badfab2aa4f3aa9fa4bbe3c05fe7d))

## [7.25.1](https://github.com/snakemake/snakemake/compare/v7.25.0...v7.25.1) (2023-04-28)


### Bug Fixes

* allow log directive in default target rule ([#2191](https://github.com/snakemake/snakemake/issues/2191)) ([86e9624](https://github.com/snakemake/snakemake/commit/86e962488dcd346cd0a29a2ff1b2dcd1abafb841))
* only consider global wildcard_constraints from the same module ([#2235](https://github.com/snakemake/snakemake/issues/2235)) ([c412b71](https://github.com/snakemake/snakemake/commit/c412b714a9fbe5cad9ad30de4a0b78b3c13068f6))
* Use `job.rule.name` attribute to fill rule field in summary ([#2217](https://github.com/snakemake/snakemake/issues/2217)) ([837c3fd](https://github.com/snakemake/snakemake/commit/837c3fd97b5a16ddb4f6b74bd2c2b5479d77bd8a))


### Documentation

* fix formatting ([087fe63](https://github.com/snakemake/snakemake/commit/087fe6307a72a577daefcea0cb150f69092138c7))
* replace `nosetest` with `pytest` ([#2211](https://github.com/snakemake/snakemake/issues/2211)) ([f6b3c47](https://github.com/snakemake/snakemake/commit/f6b3c47983bfe436f8ec33ab5830ba577fc38f90))

## [7.25.0](https://github.com/snakemake/snakemake/compare/v7.24.2...v7.25.0) (2023-03-23)


### Features

* added localrule directive ([#2180](https://github.com/snakemake/snakemake/issues/2180)) ([9c990b0](https://github.com/snakemake/snakemake/commit/9c990b076b51e9e2123ced56bfab176e03424770))
* Tes auth ([#2169](https://github.com/snakemake/snakemake/issues/2169)) ([3326a6f](https://github.com/snakemake/snakemake/commit/3326a6f7bbaf9c277400625649c3da903e251e2f))


### Bug Fixes

* always make sure that the original path of source cached files is properly passed into metadata persistence records ([#2179](https://github.com/snakemake/snakemake/issues/2179)) ([8bacbd0](https://github.com/snakemake/snakemake/commit/8bacbd0b19b0e372f5840b1fab838b79a6bab557))
* slurm batch job status queries ([#2167](https://github.com/snakemake/snakemake/issues/2167)) ([0bb69e4](https://github.com/snakemake/snakemake/commit/0bb69e429e161bf68b1f3e0b8f6fe3cbd6ed4dae))


### Documentation

* Change snakemake-tutorial download link to always be the latest ([#2183](https://github.com/snakemake/snakemake/issues/2183)) ([ae8a8f4](https://github.com/snakemake/snakemake/commit/ae8a8f4abca233e0c162e3fbd9dfb29ef38f98e6))
* fix typos in --help ([#2182](https://github.com/snakemake/snakemake/issues/2182)) ([09f0cbe](https://github.com/snakemake/snakemake/commit/09f0cbe9452a22bde768c7049eeb17f1c93db7ea))
* Improve error message when rule contains multiple run/shell/script/notebook/wrapper/template_engine/cwl keywords ([#2186](https://github.com/snakemake/snakemake/issues/2186)) ([cd5a3c4](https://github.com/snakemake/snakemake/commit/cd5a3c44eace393f24d9338b21d5ef3cb6298126))

## [7.24.2](https://github.com/snakemake/snakemake/compare/v7.24.1...v7.24.2) (2023-03-14)


### Bug Fixes

* fix index out of bounds error raised by usage of workflow.source_path called from input or params functions (thanks @AKBrueggemann) ([#2170](https://github.com/snakemake/snakemake/issues/2170)) ([cf8e6e8](https://github.com/snakemake/snakemake/commit/cf8e6e8995ecb4371c179851f4ded1d01cd1b7f9))
* limit length of failed logs decorations ([#2125](https://github.com/snakemake/snakemake/issues/2125)) ([6fc9243](https://github.com/snakemake/snakemake/commit/6fc92434f5aaed60c3d1e62bf4d33c68eeb6ed53))
* raise error if callable is passed to expand. ([#2171](https://github.com/snakemake/snakemake/issues/2171)) ([1f28476](https://github.com/snakemake/snakemake/commit/1f28476c35dc15a04e1032139b9dab5779801235))
* rounding for batch calculation ([#2064](https://github.com/snakemake/snakemake/issues/2064)) ([cbdbf9b](https://github.com/snakemake/snakemake/commit/cbdbf9b648124422ffe16c61b776f92c33c72ef8))

## [7.24.1](https://github.com/snakemake/snakemake/compare/v7.24.0...v7.24.1) (2023-03-09)


### Bug Fixes

* better job status queries for slurm executor  ([#2136](https://github.com/snakemake/snakemake/issues/2136)) ([a4df38c](https://github.com/snakemake/snakemake/commit/a4df38c56e935dde9c2745bed6afc13e1fed671f))
* get python version for script environment in a backwards compatible way that works down to python 2.7 ([#2161](https://github.com/snakemake/snakemake/issues/2161)) ([44e59b9](https://github.com/snakemake/snakemake/commit/44e59b9baa0842e19c5e0f2e05cf2fe5c9f47790))
* prevents DeprecationWarning caused by using old  draft of json schema ([#2152](https://github.com/snakemake/snakemake/issues/2152)) ([9791ffb](https://github.com/snakemake/snakemake/commit/9791ffb978a6d09d90d311cc98363c4d4efc2042))


### Performance Improvements

* Gfal2 remote provider using gfal2-python instead of  gfal2-utils. ([#2128](https://github.com/snakemake/snakemake/issues/2128)) ([0b9bfe5](https://github.com/snakemake/snakemake/commit/0b9bfe500a06e669d2557d897ed26550aec526d6))


### Documentation

* fix minor typos in a linting rule ([#2162](https://github.com/snakemake/snakemake/issues/2162)) ([71e1171](https://github.com/snakemake/snakemake/commit/71e1171a0dab824baed77bfe235268af9e095c1f))

## [7.24.0](https://github.com/snakemake/snakemake/compare/v7.23.1...v7.24.0) (2023-03-01)


### Features

* limit the number of input/output files in job properties ([#2149](https://github.com/snakemake/snakemake/issues/2149)) ([d93f091](https://github.com/snakemake/snakemake/commit/d93f091acea63a662dcb350c3f86c15fa9bdf721))


### Bug Fixes

* [#2130](https://github.com/snakemake/snakemake/issues/2130) by patching the protect() method so the path of files in subdirectories is properly resolved during write-protection ([#2131](https://github.com/snakemake/snakemake/issues/2131)) ([1a754fd](https://github.com/snakemake/snakemake/commit/1a754fd094bd13bb4a201f1c80a077656c89f995))
* `sre_constants` import because of deprecation ([#2139](https://github.com/snakemake/snakemake/issues/2139)) ([3b326db](https://github.com/snakemake/snakemake/commit/3b326dba22ef5358092c281479eafafe3480eeae))
* ensure user and group rw permissions for metadata files and source cache ([#2132](https://github.com/snakemake/snakemake/issues/2132)) ([cc51faa](https://github.com/snakemake/snakemake/commit/cc51faaa7d4f20896fc46b9fd67d062936d641bb))
* is_run error with local, group jobs ([#2133](https://github.com/snakemake/snakemake/issues/2133)) ([31bfcd5](https://github.com/snakemake/snakemake/commit/31bfcd5399540fc6cf52e3b76144e9abea6d4eab))
* require toposort &gt;= 1.10 ([#2145](https://github.com/snakemake/snakemake/issues/2145)) ([3cb54b8](https://github.com/snakemake/snakemake/commit/3cb54b8c62743897f20feb3fcf269a7357878434))


### Documentation

* Update modularization.rst ([#2137](https://github.com/snakemake/snakemake/issues/2137)) ([16954c7](https://github.com/snakemake/snakemake/commit/16954c7b633049df6646275139251097d574fd35))

## [7.23.1](https://github.com/snakemake/snakemake/compare/v7.23.0...v7.23.1) (2023-02-18)


### Bug Fixes

* batch collect jobs for scancel ([#2114](https://github.com/snakemake/snakemake/issues/2114)) ([0b1fe31](https://github.com/snakemake/snakemake/commit/0b1fe312e8c98a814b1c419940f35253f58f958e))

## [7.23.0](https://github.com/snakemake/snakemake/compare/v7.22.0...v7.23.0) (2023-02-18)


### Features

* changed report layout to display menu always left of the results. For fullscreen, one can still hide the menu, which leads to automatic growth of the results ([#2116](https://github.com/snakemake/snakemake/issues/2116)) ([d771b1b](https://github.com/snakemake/snakemake/commit/d771b1b5fc8344aaffe1f30388d4e4d31d4fe937))
* Publish docker images for amd64 & arm64 ([#2105](https://github.com/snakemake/snakemake/issues/2105)) ([4c898f5](https://github.com/snakemake/snakemake/commit/4c898f5587d832c45f0b534681f9502abe1de6ce))


### Bug Fixes

* use text/markdown for long_description_content_type ([#2112](https://github.com/snakemake/snakemake/issues/2112)) ([0241075](https://github.com/snakemake/snakemake/commit/02410755c51df21833199db70406b2179248380e))


### Performance Improvements

* Improve execution speed of cleanup_workdir (in dag)  ([#2103](https://github.com/snakemake/snakemake/issues/2103)) ([1fbc5f5](https://github.com/snakemake/snakemake/commit/1fbc5f5aee65bc8dd776765644d07051dd857670))

## [7.22.0](https://github.com/snakemake/snakemake/compare/v7.21.0...v7.22.0) (2023-02-12)


### Features

* add cleanup containers option ([#2088](https://github.com/snakemake/snakemake/issues/2088)) ([053e3b3](https://github.com/snakemake/snakemake/commit/053e3b37cfcc6c67bae6ac3660b82879b75acb4c))


### Bug Fixes

* assume shared filesystem by default when running with --flux ([#2075](https://github.com/snakemake/snakemake/issues/2075)) ([4bec2fd](https://github.com/snakemake/snakemake/commit/4bec2fd3bb0c48a1f38506a966cb64dc8c2d1021))
* properly handle NA values for paramspaces ([#2098](https://github.com/snakemake/snakemake/issues/2098)) ([6b6a880](https://github.com/snakemake/snakemake/commit/6b6a88074eac6e3a9aa8c89501fc9481f07ecc1d))

## [7.21.0](https://github.com/snakemake/snakemake/compare/v7.20.0...v7.21.0) (2023-01-30)


### Features

* ability to encode paramspaces into a single wildcard, via the newly introduced `single_wildcard` argument of `Paramspace`. ([#2069](https://github.com/snakemake/snakemake/issues/2069)) ([728ab3c](https://github.com/snakemake/snakemake/commit/728ab3cf59fb188bf1872a5cab3aba4519340a06))
* allow input, output, and params to be used in functions passed to report mark arguments ([#2081](https://github.com/snakemake/snakemake/issues/2081)) ([93ff8b6](https://github.com/snakemake/snakemake/commit/93ff8b604a152f50ca31389ecffe1a7527c5b5a8))


### Bug Fixes

* more robust encoding of params in persistent metadata storage. This way, pandas parameters do not lead to spurious rerun triggers. ([#2080](https://github.com/snakemake/snakemake/issues/2080)) ([106a4c3](https://github.com/snakemake/snakemake/commit/106a4c3f4b181d787a6c309b8ea2e655780413e7))
* more robust parsing of sacct output in slurm executor ([#2036](https://github.com/snakemake/snakemake/issues/2036)) ([fe651f8](https://github.com/snakemake/snakemake/commit/fe651f8a10b9ead94b07ab31efe2d560525fc3b6))
* Postprocess job groups in toposorted order for correct touch times ([#2073](https://github.com/snakemake/snakemake/issues/2073)) ([10b5849](https://github.com/snakemake/snakemake/commit/10b584916a869fa1147a9f42d1de4fec4120b441))

## [7.20.0](https://github.com/snakemake/snakemake/compare/v7.19.1...v7.20.0) (2023-01-18)


### Features

* add tes token ([#1966](https://github.com/snakemake/snakemake/issues/1966)) ([59a8fa0](https://github.com/snakemake/snakemake/commit/59a8fa04c4c6b113775fe11228b82510ecd36cb8))
* Add token auth to GitLab/GitHub hosting providers ([#1761](https://github.com/snakemake/snakemake/issues/1761)) ([e03a3b4](https://github.com/snakemake/snakemake/commit/e03a3b42eea89d512290bf98ee7d77ce2e17447c)), closes [#1301](https://github.com/snakemake/snakemake/issues/1301)
* allow for human friendly resource definitions (e.g. mem="5GB", runtime="1d") which deprecates slurm constrained time format (e.g. runtime="hh:mm:ss") ([#1861](https://github.com/snakemake/snakemake/issues/1861)) ([24610ac](https://github.com/snakemake/snakemake/commit/24610ac75849d543fc38c83fb2454fa4f9b42075)) ([#2154](https://github.com/snakemake/snakemake/issues/2154))



### Bug Fixes

* :bug: - fix hyperlink ([#2046](https://github.com/snakemake/snakemake/issues/2046)) ([9519d31](https://github.com/snakemake/snakemake/commit/9519d31b4ed2390d5c14f0f6a754ca665374d15c))
* Catch missing error stream in Slurm executor ([#2063](https://github.com/snakemake/snakemake/issues/2063)) ([c21fc7e](https://github.com/snakemake/snakemake/commit/c21fc7e528327b13d762c5db90ee0e40506cf0bd))
* correctly parse empty values in config cli ([#2032](https://github.com/snakemake/snakemake/issues/2032)) ([1b0689d](https://github.com/snakemake/snakemake/commit/1b0689dddf5eecbd8afd307c6df3dc31a32c338f))
* Correctly parse UserDicts in executors ([#2016](https://github.com/snakemake/snakemake/issues/2016)) ([e3926fa](https://github.com/snakemake/snakemake/commit/e3926fa4b44bedf99745b949080397919a522aa1))
* Fix handling of --jobs in no-exec state ([#2029](https://github.com/snakemake/snakemake/issues/2029)) ([e8e8222](https://github.com/snakemake/snakemake/commit/e8e8222a6f2192423aa55766304f9d1616a0d6e3))
* make `--show-failed-logs` handle empty log files ([#2039](https://github.com/snakemake/snakemake/issues/2039)) ([683c6f2](https://github.com/snakemake/snakemake/commit/683c6f2284e867457d7ed25e838ed3018da8f2d4)), closes [#2023](https://github.com/snakemake/snakemake/issues/2023)
* make python version check more robust ([#2058](https://github.com/snakemake/snakemake/issues/2058)) ([e685621](https://github.com/snakemake/snakemake/commit/e685621f3b4c0e21aa6f640dd571406d1d39e588))
* parsing error when last line is comment ([#2054](https://github.com/snakemake/snakemake/issues/2054)) ([a928dd4](https://github.com/snakemake/snakemake/commit/a928dd4391d186b6ddb582331047f09fbee03f03))
* prevent overriding of retries when set to 0 ([#2053](https://github.com/snakemake/snakemake/issues/2053)) ([a328f3e](https://github.com/snakemake/snakemake/commit/a328f3e009a53002781539dfe3e3c8c5d1738189))
* propagate attempt count from group to subjobs ([#2052](https://github.com/snakemake/snakemake/issues/2052)) ([da3f1c0](https://github.com/snakemake/snakemake/commit/da3f1c0a80ffabbe9d02ce1361bfa2374c546007))
* remove overflow from rulegraph div in report ([9a0aaa7](https://github.com/snakemake/snakemake/commit/9a0aaa703d6531d890b7116638dc515425f6ed34))
* skip type checks of missing dir in touch mode ([#2051](https://github.com/snakemake/snakemake/issues/2051)) ([ae00c25](https://github.com/snakemake/snakemake/commit/ae00c25541600dba14ad68f1145bab3c4455de19))
* slurm default_resources quoting ([#2043](https://github.com/snakemake/snakemake/issues/2043)) ([47d3fc3](https://github.com/snakemake/snakemake/commit/47d3fc3eef8ebe54df8b77f100fc2ab3fa36c190))
* Update list of python versions in classifiers ([#2020](https://github.com/snakemake/snakemake/issues/2020)) ([7a98100](https://github.com/snakemake/snakemake/commit/7a98100ba92b0174c8ead3ae715042c1ab710c61))
* use short argument name for `--chdir` for compatibility with Slurm &lt;=v17 ([#2040](https://github.com/snakemake/snakemake/issues/2040)) ([a9ed3ec](https://github.com/snakemake/snakemake/commit/a9ed3ec3e823442810388fce8a17ba4950bbdaa2))
* human friendly resource definitions introduce inability to parse slurm specific time format (e.g. "hh:mm:ss"). New time format (e.g. "1d") adds portability among various job schedulers and clusters ([#2154](https://github.com/snakemake/snakemake/issues/2154))


### Documentation

* Fix typo in SLURM help text ([#2049](https://github.com/snakemake/snakemake/issues/2049)) ([79b7025](https://github.com/snakemake/snakemake/commit/79b702528b9801d0283766f2da377d4b9daeebd2))
* mention XDG_CACHE_HOME ([#2057](https://github.com/snakemake/snakemake/issues/2057)) ([ec2ef45](https://github.com/snakemake/snakemake/commit/ec2ef45c5e49f97a0b58535f8e39152d2789d428))

## [7.19.1](https://github.com/snakemake/snakemake/compare/v7.19.0...v7.19.1) (2022-12-13)


### Bug Fixes

* improved default resources parsing (also allowing to deactivate a default resource via setting it to None) ([#2006](https://github.com/snakemake/snakemake/issues/2006)) ([e6cdb32](https://github.com/snakemake/snakemake/commit/e6cdb32fd20a9899f441b3ec8aca3f36710f7a4a))


### Documentation

* fix link ([4889c93](https://github.com/snakemake/snakemake/commit/4889c93a031acc3ada22c9bdfd27301cce2c107e))
* fix typo ([e1c3cc6](https://github.com/snakemake/snakemake/commit/e1c3cc6c6ddbf67e9d3cfc9dcef9bcf35b0060f3))
* fix typos ([e45b9e6](https://github.com/snakemake/snakemake/commit/e45b9e608c31af5a5c3389520486d46935549eb4))
* fix typos ([151095d](https://github.com/snakemake/snakemake/commit/151095ddc839860d33a27cafeb83260dc17d736b))
* format table ([4180a1b](https://github.com/snakemake/snakemake/commit/4180a1b6a8b7104e27a748c638275b98c5998200))
* polished text and table display ([413356c](https://github.com/snakemake/snakemake/commit/413356cceb8d3a96b24531ef350168e681c9d383))

## [7.19.0](https://github.com/snakemake/snakemake/compare/v7.18.2...v7.19.0) (2022-12-13)


### Features

* add  keyword to gridftp remote provide to specify the number or disable usage of multiple data stream ([#1974](https://github.com/snakemake/snakemake/issues/1974)) ([3e6675d](https://github.com/snakemake/snakemake/commit/3e6675df26bf65fa27006ea57a4c3cf36b89d6da))
* provide information about temp, pipe, and service files in --summary ([#1977](https://github.com/snakemake/snakemake/issues/1977)) ([c7c7776](https://github.com/snakemake/snakemake/commit/c7c7776f8722adf94e6a176174cb0a7564f11d9f))
* native SLURM support (--slurm, see docs) ([#1015](https://github.com/snakemake/snakemake/issues/1015)) ([c7ea059](https://github.com/snakemake/snakemake/commit/c7ea0590c396c67fa5d56042e21f678c20784d3b))


### Bug Fixes

* avoid logfile writing in case of dryrun; better hints in case of incomplete checkpoints ([#1994](https://github.com/snakemake/snakemake/issues/1994)) ([a022705](https://github.com/snakemake/snakemake/commit/a022705db14c9409b1ceadf6d5ae6367833e2131))
* handle case where zenodo deposition does not return files ([#2004](https://github.com/snakemake/snakemake/issues/2004)) ([b63c4a7](https://github.com/snakemake/snakemake/commit/b63c4a7e496ca7a3353e8a59b8ba493d65156cb5))
* issue [#1882](https://github.com/snakemake/snakemake/issues/1882) WorkflowError: Metadata can't be created as it already exists (Windows) ([#1971](https://github.com/snakemake/snakemake/issues/1971)) ([d4484e6](https://github.com/snakemake/snakemake/commit/d4484e61ef49be23fc2bef8bf879185e521f5376))
* json validation error with markdown cells ([#1986](https://github.com/snakemake/snakemake/issues/1986)) ([6c26f75](https://github.com/snakemake/snakemake/commit/6c26f757785226d796c80689fabae771f316af9f))

## [7.18.2](https://github.com/snakemake/snakemake/compare/v7.18.1...v7.18.2) (2022-11-10)


### Bug Fixes

* Change ratelimiter dependency to throttler ([#1958](https://github.com/snakemake/snakemake/issues/1958)) ([50b8f16](https://github.com/snakemake/snakemake/commit/50b8f1609a597dc9f25d2fd86c9cdda531bdc041))
* fixed problem with leaked modifications when inheriting multiple times from the same rule ([#1957](https://github.com/snakemake/snakemake/issues/1957)) ([2475cbc](https://github.com/snakemake/snakemake/commit/2475cbcd74d9c9f62e07617751979bb00025850a))
* forwarding --keep-incomplete to cluster executor ([#1951](https://github.com/snakemake/snakemake/issues/1951)) ([2894c7d](https://github.com/snakemake/snakemake/commit/2894c7dfaa854ebe34b1248897c90b4110d3962b))
* show input files on job error ([#1949](https://github.com/snakemake/snakemake/issues/1949)) ([ad21631](https://github.com/snakemake/snakemake/commit/ad2163187f031317d889e0cbe368176e1d48d13f))

## [7.18.1](https://github.com/snakemake/snakemake/compare/v7.18.0...v7.18.1) (2022-11-03)


### Bug Fixes

* regression ValueError introduced with 7.17.2 ([#1947](https://github.com/snakemake/snakemake/issues/1947)) ([53a4fca](https://github.com/snakemake/snakemake/commit/53a4fca8c67a3b58d61b146c8cfff3982889d77d))

## [7.18.0](https://github.com/snakemake/snakemake/compare/v7.17.2...v7.18.0) (2022-10-31)


### Features

* first try to match output files against input files while persisting wildcard values from the consuming job. This can dramatically reduce ambiuity problems. Thanks to [@descostesn](https://github.com/descostesn)! ([#1939](https://github.com/snakemake/snakemake/issues/1939)) ([d093907](https://github.com/snakemake/snakemake/commit/d093907417778c7693a05ed1f38fc40b8d34d9ba))

## [7.17.2](https://github.com/snakemake/snakemake/compare/v7.17.1...v7.17.2) (2022-10-28)


### Bug Fixes

* Consider source cache when setting search path for python scripts. This allows to import from Python modules next to scripts while deploying the workflow as a snakemake module, even from remote locations. ([#1940](https://github.com/snakemake/snakemake/issues/1940)) ([27be1d4](https://github.com/snakemake/snakemake/commit/27be1d41c397a974f33dcf93ccce331a80ab0198))

## [7.17.1](https://github.com/snakemake/snakemake/compare/v7.17.0...v7.17.1) (2022-10-28)


### Bug Fixes

* change source cache entries to keep the original name and folder structure, such that imports from e.g. scripts also work with remote modules (if specified as additional input files with workflow.source_path) ([#1936](https://github.com/snakemake/snakemake/issues/1936)) ([c34f3f6](https://github.com/snakemake/snakemake/commit/c34f3f64ac19d2c2eaab361d73d3144430538bb6))

## [7.17.0](https://github.com/snakemake/snakemake/compare/v7.16.2...v7.17.0) (2022-10-27)


### Features

* allow to define the cache mode per rule (this enables to exclude software envs from the caching hash value, which can be handy e.g. for download rules where the software version does not affect the result) ([#1933](https://github.com/snakemake/snakemake/issues/1933)) ([715e618](https://github.com/snakemake/snakemake/commit/715e6187e9a132c2a61f5bef34a1e10491680b0a))


### Performance Improvements

* cached os.pathconf() call in _record_path() ([#1920](https://github.com/snakemake/snakemake/issues/1920)) ([551badb](https://github.com/snakemake/snakemake/commit/551badb4eb9de582278185449d8fa9298ad7ae8c))

## [7.16.2](https://github.com/snakemake/snakemake/compare/v7.16.1...v7.16.2) (2022-10-26)


### Bug Fixes

* fix false rerun triggering downstream of checkpoints due to spurious parameter, code or software env changes ([638ea86](https://github.com/snakemake/snakemake/commit/638ea86c51ea5746c8b5452cc8a0e43108de15ef))
* remove redundant dot in expand call in multiext documentation ([#1921](https://github.com/snakemake/snakemake/issues/1921)) ([278beaa](https://github.com/snakemake/snakemake/commit/278beaa7c81b9e418fa42ac5f944e3d7e2cdfbbd))

## [7.16.1](https://github.com/snakemake/snakemake/compare/v7.16.0...v7.16.1) (2022-10-18)


### Bug Fixes

* conda create --no-shortcuts absent on Linux/MacOS (regression from [#1046](https://github.com/snakemake/snakemake/issues/1046)) ([#1916](https://github.com/snakemake/snakemake/issues/1916)) ([8a86a1e](https://github.com/snakemake/snakemake/commit/8a86a1e5ec7438492e2f1403b1b0fd81030255ad))
* fix typo in line display of exceptions ([#1912](https://github.com/snakemake/snakemake/issues/1912)) ([55e38a6](https://github.com/snakemake/snakemake/commit/55e38a6f80c74b24d7763975b0ae2826d75f23d9))

## [7.16.0](https://github.com/snakemake/snakemake/compare/v7.15.2...v7.16.0) (2022-10-14)


### Features

* k8s: add --k8s-cpu-scalar ([#1857](https://github.com/snakemake/snakemake/issues/1857)) ([a067a1b](https://github.com/snakemake/snakemake/commit/a067a1b6eb8f6432348bc257782faba40f89e805))


### Bug Fixes

* allow report generation to handle pathlib objects ([#1904](https://github.com/snakemake/snakemake/issues/1904)) ([7c34656](https://github.com/snakemake/snakemake/commit/7c346569106a36bbbce990576324af472ada9efd))
* fix false reruns after checkpoints ([#1907](https://github.com/snakemake/snakemake/issues/1907)) ([dc5af12](https://github.com/snakemake/snakemake/commit/dc5af12f54e774612bb1f8ead45e4597080dc100))

## [7.15.2](https://github.com/snakemake/snakemake/compare/v7.15.1...v7.15.2) (2022-10-08)


### Bug Fixes

* Comparison of rules and non-rule instances ([#1894](https://github.com/snakemake/snakemake/issues/1894)) ([bf01ece](https://github.com/snakemake/snakemake/commit/bf01ece0e9c51442daba02ecf2ef37aa276283d6))
* delay evaluation of tmpdir to actual job execution, and not submission. This way, tmpdir can be dependent on the node context. ([#1860](https://github.com/snakemake/snakemake/issues/1860)) ([4203556](https://github.com/snakemake/snakemake/commit/420355662ebea0bc72c9bd6bca1eea5259f3b43e))
* ensure that rule name string instead of object is passed to tabulate package ([#1898](https://github.com/snakemake/snakemake/issues/1898)) ([f9ff157](https://github.com/snakemake/snakemake/commit/f9ff157c5c534ba035bdf51a02fbbba5ad94dd61))
* issue 1846 ([#1888](https://github.com/snakemake/snakemake/issues/1888)) ([da2dfbd](https://github.com/snakemake/snakemake/commit/da2dfbd765aa1e2a8da36d9eaa1ac9fcffa5e921))
* lexicographically sorted rule display with --list, and trimmed rule docstrings ([#1880](https://github.com/snakemake/snakemake/issues/1880)) ([32128ae](https://github.com/snakemake/snakemake/commit/32128ae118d15f250e4b438735e30151cd6f27c5))


### Performance Improvements

* Average NamedList __getitem__ performance improvement ([#1825](https://github.com/snakemake/snakemake/issues/1825)) ([10451b7](https://github.com/snakemake/snakemake/commit/10451b7198a4a39149ce5e8ec82c17df1f18813b))

## [7.15.1](https://github.com/snakemake/snakemake/compare/v7.15.0...v7.15.1) (2022-10-04)


### Bug Fixes

* fix `--immediate-submit` ([#1851](https://github.com/snakemake/snakemake/issues/1851)) ([e358372](https://github.com/snakemake/snakemake/commit/e3583721f2dc620ce96876ecec58846d1cbe7bfd))
* Handle temp files for all jobs in a group. ([#1779](https://github.com/snakemake/snakemake/issues/1779)) ([d28b893](https://github.com/snakemake/snakemake/commit/d28b89363f007d303d733b2b12f517502867035c))


### Documentation

* small tweaks to flux documentation ([#1886](https://github.com/snakemake/snakemake/issues/1886)) ([f29b371](https://github.com/snakemake/snakemake/commit/f29b37106727c470b691076189e92f35e4cecfb6))
* various little fixes ([#1875](https://github.com/snakemake/snakemake/issues/1875)) ([b93f8e3](https://github.com/snakemake/snakemake/commit/b93f8e316cb2f51e72933e7a28872bdf523aec11))

## [7.15.0](https://github.com/snakemake/snakemake/compare/v7.14.2...v7.15.0) (2022-10-04)


### Features

* adding flux executor ([#1810](https://github.com/snakemake/snakemake/issues/1810)) ([40d2bd0](https://github.com/snakemake/snakemake/commit/40d2bd071984914ac511e7858690dfd16cefaf69))


### Bug Fixes

* Add back logging of run directives ([#1883](https://github.com/snakemake/snakemake/issues/1883)) ([a65559c](https://github.com/snakemake/snakemake/commit/a65559c1e21d65e5b4509b9565b472e734ab9f02))


### Documentation

* fix grammar in the intro ([#1859](https://github.com/snakemake/snakemake/issues/1859)) ([774bc6a](https://github.com/snakemake/snakemake/commit/774bc6aaa3105c94a551691498d9a5efb13ac216))
* fix typo ([#1843](https://github.com/snakemake/snakemake/issues/1843)) ([6572ad9](https://github.com/snakemake/snakemake/commit/6572ad91bed14dece6b01a26134007a25ef0c4b2))

## [7.14.2](https://github.com/snakemake/snakemake/compare/v7.14.1...v7.14.2) (2022-09-26)


### Bug Fixes

* reduce resource requirements for kubernetes tests ([#1876](https://github.com/snakemake/snakemake/issues/1876)) ([cb4b78a](https://github.com/snakemake/snakemake/commit/cb4b78a05ee08f7dafb561ba33bbe460ec097eb5))

## [7.14.1](https://github.com/snakemake/snakemake/compare/v7.14.0...v7.14.1) (2022-09-23)


### Bug Fixes

* allocation of local ssds in k8s tests ([#1870](https://github.com/snakemake/snakemake/issues/1870)) ([d0de4dc](https://github.com/snakemake/snakemake/commit/d0de4dccf6f6c749da0d4a30ef27fe0f9274995d))
* allow script directive to take pathlib Path ([#1869](https://github.com/snakemake/snakemake/issues/1869)) ([12cdc96](https://github.com/snakemake/snakemake/commit/12cdc961c86541d5c50533a6b8ff5df0cb6fd7d1))
* catch errors in remote.AUTO provider list ([#1834](https://github.com/snakemake/snakemake/issues/1834)) ([c613ed2](https://github.com/snakemake/snakemake/commit/c613ed217f1dfb16fc63fa7b06af4ddf4f3dd0b8))
* consistently use text output in conda shell commands and various little fixes for failing test cases due to conda package changes ([#1864](https://github.com/snakemake/snakemake/issues/1864)) ([4234fe7](https://github.com/snakemake/snakemake/commit/4234fe765e996dbb3a1738567299ebb4d8c28af0))
* declare associative arrays ([#1844](https://github.com/snakemake/snakemake/issues/1844)) ([90ae449](https://github.com/snakemake/snakemake/commit/90ae44943a6b7712d46677364752bf0a7cf91806))
* fix falsely triggered reruns if input files are obtained via workflow.source_path() ([#1862](https://github.com/snakemake/snakemake/issues/1862)) ([2dc2e6a](https://github.com/snakemake/snakemake/commit/2dc2e6aa3255b8b73e06b7af1ea646d5081799ce))
* fixed typos ([#1847](https://github.com/snakemake/snakemake/issues/1847)) ([a1e49b6](https://github.com/snakemake/snakemake/commit/a1e49b6f290daeb2012575ab3a85ca72cfe42747))
* k8s container volume mounts as list ([#1868](https://github.com/snakemake/snakemake/issues/1868)) ([5c54df3](https://github.com/snakemake/snakemake/commit/5c54df39c6a111dcfa7adaea15d0b81a3fc16b90))
* None type error when invoking Workflow object manually ([#1731](https://github.com/snakemake/snakemake/issues/1731)) ([dc45ccb](https://github.com/snakemake/snakemake/commit/dc45ccb9ee94e8fcdf386c8598bbae57e319ba10))
* request disk_mb resource from k8s ([#1858](https://github.com/snakemake/snakemake/issues/1858)) ([f68f166](https://github.com/snakemake/snakemake/commit/f68f166582abaeb45a1a093306626bb8abb0e0bb))
* respect shebang lines in post-deploy scripts (see deployment docs) ([#1841](https://github.com/snakemake/snakemake/issues/1841)) ([c26c4b6](https://github.com/snakemake/snakemake/commit/c26c4b6ff43e38797168ef7983c40b0c8a4b2f8c))

## [7.14.0](https://github.com/snakemake/snakemake/compare/v7.13.0...v7.14.0) (2022-08-27)


### Features

* add support for bash scripts in the script directive (beyond small shell commands) ([#1821](https://github.com/snakemake/snakemake/issues/1821)) ([c4cf8fd](https://github.com/snakemake/snakemake/commit/c4cf8fdb119e678e51fd1932392699957c63b6c4))


### Documentation

* fix small typo in FAQ ([#1832](https://github.com/snakemake/snakemake/issues/1832)) ([914172b](https://github.com/snakemake/snakemake/commit/914172baf366570fa2a2818746d3ed1417790e3a))

## [7.13.0](https://github.com/snakemake/snakemake/compare/v7.12.1...v7.13.0) (2022-08-25)


### Features

* add gitfile option to make it possible to use local git repos when importing modules ([#1376](https://github.com/snakemake/snakemake/issues/1376)) ([1a3b91f](https://github.com/snakemake/snakemake/commit/1a3b91fa8461d44bc1583f1d2e33ff2bae1360b3))


### Bug Fixes

* allow to use {wildcards} for group jobs in cluster config ([#1555](https://github.com/snakemake/snakemake/issues/1555)) ([f0ec73d](https://github.com/snakemake/snakemake/commit/f0ec73d2470d3ecf1d8500b52fa86021463a96fe))
* avoid "Admin" prompt when using conda on windows ([#1046](https://github.com/snakemake/snakemake/issues/1046)) ([552fadf](https://github.com/snakemake/snakemake/commit/552fadf7d0ea0f0ff3990295db2618a31c692f61))
* handle benmark bug that arise with singularity ([#1671](https://github.com/snakemake/snakemake/issues/1671)) ([10ef7c4](https://github.com/snakemake/snakemake/commit/10ef7c4d1cbbd89fa1760ec9f643e09ae4eb8bd9))
* Open Snakefile for reading with explicit encoding specified ([#1146](https://github.com/snakemake/snakemake/issues/1146)) ([ec1d859](https://github.com/snakemake/snakemake/commit/ec1d859dd2422fa551bbf58843907fc5eed244ff))
* remove superfluous comma causing TypeError in conda-frontend error message ([#1804](https://github.com/snakemake/snakemake/issues/1804)) ([87b013c](https://github.com/snakemake/snakemake/commit/87b013cbba420b71892ab9f482718c1cd59e6bda))


### Documentation

* explain SNAKEMAKE_PROFILE environment variable ([2b32bba](https://github.com/snakemake/snakemake/commit/2b32bba6efa53a9359bc67d8eb8ce06f1106e8ef))
* update contribution docs ([09a5595](https://github.com/snakemake/snakemake/commit/09a559504d7a4346e40d8b3a99ece46e903ff7a0))

## [7.12.1](https://github.com/snakemake/snakemake/compare/v7.12.0...v7.12.1) (2022-08-09)


### Bug Fixes

* Fix case of multiple scattergather processes ([#1799](https://github.com/snakemake/snakemake/issues/1799)) ([417aad4](https://github.com/snakemake/snakemake/commit/417aad4c1a5fe83bdf24729d0d4e70df112fd293))
* more comprehensive error reporting for RuleExceptions ([#1802](https://github.com/snakemake/snakemake/issues/1802)) ([1cd9512](https://github.com/snakemake/snakemake/commit/1cd95129a92e522e11689620473f68b1fe69fb42))

## [7.12.0](https://github.com/snakemake/snakemake/compare/v7.11.0...v7.12.0) (2022-07-29)


### Features

* print reason summary in case of dryrun ([#1778](https://github.com/snakemake/snakemake/issues/1778)) ([bd2a68b](https://github.com/snakemake/snakemake/commit/bd2a68bb8a7ddb64ad91bb822aa7d0342e884704))


### Bug Fixes

* Fix technical bugs in resource-scope documentation ([#1784](https://github.com/snakemake/snakemake/issues/1784)) ([878420c](https://github.com/snakemake/snakemake/commit/878420c4496562df7fe5fed4e6f877e670dcd533))
* move max_status_checks_per_second attribute setting before the wait thread of cluster backends is started to avoid missing attribute errors ([#1775](https://github.com/snakemake/snakemake/issues/1775)) ([a48e9d0](https://github.com/snakemake/snakemake/commit/a48e9d035efec3a91b50bbda43ddf92196d5084b))

## [7.11.0](https://github.com/snakemake/snakemake/compare/v7.10.0...v7.11.0) (2022-07-27)


### Features

* improved resource handling in groups and ability to define resource scopes (global or per node), see docs and --help ([#1218](https://github.com/snakemake/snakemake/issues/1218)) ([a8014d0](https://github.com/snakemake/snakemake/commit/a8014d030a2a3ea04743e30bcf5164801291378f))


### Bug Fixes

* fixed conda frontend detection and checking to also work with latest mambaforge ([#1781](https://github.com/snakemake/snakemake/issues/1781)) ([225e68c](https://github.com/snakemake/snakemake/commit/225e68cb05ecbf8fe001a884a416da0239cfd4d1))

## [7.10.0](https://github.com/snakemake/snakemake/compare/v7.9.0...v7.10.0) (2022-07-26)


### Features

* Support conda environment definitions to be passed as function pointers, similar to input, params, and resources ([#1300](https://github.com/snakemake/snakemake/issues/1300)) ([6f582f1](https://github.com/snakemake/snakemake/commit/6f582f180f5df9292e005623f093a3a0b2e597a7))


### Bug Fixes

* fix regression in workflow source acquisition of google life science executor ([#1773](https://github.com/snakemake/snakemake/issues/1773)) ([c07732e](https://github.com/snakemake/snakemake/commit/c07732e2c7d6e95fabb3f2570bf0a1b0af2ac2cb))
* limit filename length of temporary files generated by the persistence backend (metadata, incomplete markers, etc.) ([#1780](https://github.com/snakemake/snakemake/issues/1780)) ([59053e7](https://github.com/snakemake/snakemake/commit/59053e7017584ea94860cecb5ebecce66aed14ce))

## [7.9.0](https://github.com/snakemake/snakemake/compare/v7.8.5...v7.9.0) (2022-07-19)


### Features

* make it possible to exclude rules that will be imported when using 'use rule' statement ([#1717](https://github.com/snakemake/snakemake/issues/1717)) ([d9e0611](https://github.com/snakemake/snakemake/commit/d9e061178bd22307cc710bea28a5994e866260d9))


### Bug Fixes

* add lock free mechanism for avoiding race conditions when writing persistence information; consider corrupt metadata records as non-existent ([#1745](https://github.com/snakemake/snakemake/issues/1745)) ([71fe952](https://github.com/snakemake/snakemake/commit/71fe9527bb7011ba01d25fdd21c102c135412c04))
* conda python interpreter path on Windows ([#1711](https://github.com/snakemake/snakemake/issues/1711)) ([155c9d6](https://github.com/snakemake/snakemake/commit/155c9d6688a99db4b49c5b25d5a8a65a3aca532d))
* ensures that REncoder also checks for numpy.bool_ in encode_value ([#1749](https://github.com/snakemake/snakemake/issues/1749)) ([10a6e1d](https://github.com/snakemake/snakemake/commit/10a6e1de50ea7957e8685ba55b2cc115101cc23f))
* Move quiet default after profile parsing ([#1764](https://github.com/snakemake/snakemake/issues/1764)) ([6ade76d](https://github.com/snakemake/snakemake/commit/6ade76d1287c8f62056853491a4e67b08a4739a6))

## [7.8.5](https://github.com/snakemake/snakemake/compare/v7.8.4...v7.8.5) (2022-06-30)


### Documentation

* fix long description type for pypi (set to markdown) ([d8d9b8f](https://github.com/snakemake/snakemake/commit/d8d9b8f284d2863f672cd96fd96d2034806d52de))

## [7.8.4](https://github.com/snakemake/snakemake/compare/v7.8.3...v7.8.4) (2022-06-30)


### Bug Fixes

* only display a warning in case of non-strict channel priorities ([#1752](https://github.com/snakemake/snakemake/issues/1752)) ([b84fa33](https://github.com/snakemake/snakemake/commit/b84fa337f5c417360fc744202a419da3f37594b5))
* pass triggers and resources to subworkflow ([#1733](https://github.com/snakemake/snakemake/issues/1733)) ([fa7fb75](https://github.com/snakemake/snakemake/commit/fa7fb75d2314cd2965a9f578506e27aae4b07349))
* add pyproject.toml to use setuptools features ([#1725](https://github.com/snakemake/snakemake/issues/1725)) ([454bfd1](https://github.com/snakemake/snakemake/commit/454bfd17aba8631ea34d6a8971b1db60b2955965))


### Documentation

* add workflows.community metadata ([#1736](https://github.com/snakemake/snakemake/issues/1736)) ([8a42afc](https://github.com/snakemake/snakemake/commit/8a42afc053e3f398f1a299fb3234ca0ee847faae))

### [7.8.3](https://www.github.com/snakemake/snakemake/compare/v7.8.2...v7.8.3) (2022-06-20)


### Bug Fixes

* allow apptainer as a successor to singularity. ([#1706](https://www.github.com/snakemake/snakemake/issues/1706)) ([bcbdb0b](https://www.github.com/snakemake/snakemake/commit/bcbdb0bc44746961838f9603e42d0e171407c740))
* improved provenance trigger info ([#1720](https://www.github.com/snakemake/snakemake/issues/1720)) ([29d959d](https://www.github.com/snakemake/snakemake/commit/29d959d86341aee66e945f216cae41e9c531a4d1))
* small changes to make docs checkpoint example functional ([#1714](https://www.github.com/snakemake/snakemake/issues/1714)) ([1d4909e](https://www.github.com/snakemake/snakemake/commit/1d4909ef838ac79f8e7be9f29d76b969d358ef1b))

### [7.8.2](https://www.github.com/snakemake/snakemake/compare/v7.8.1...v7.8.2) (2022-06-08)


### Bug Fixes

* fixed bug in needrun computation of jobs downstream of checkpoints ([#1704](https://www.github.com/snakemake/snakemake/issues/1704)) ([c634b78](https://www.github.com/snakemake/snakemake/commit/c634b78b4d7c4f6ef59e46c94162893e42de6f73))

### [7.8.1](https://www.github.com/snakemake/snakemake/compare/v7.8.0...v7.8.1) (2022-05-31)


### Bug Fixes

* handling of remaining jobs when using --keep-going ([#1693](https://www.github.com/snakemake/snakemake/issues/1693)) ([87e4303](https://www.github.com/snakemake/snakemake/commit/87e430354b45c702e4d782bad8461fde44477a48))
* more robust calculation of number of jobs until ready for execution ([#1691](https://www.github.com/snakemake/snakemake/issues/1691)) ([fdfc717](https://www.github.com/snakemake/snakemake/commit/fdfc717f8764ce535801d06805aaa1d59186ec84))
* propagate rerun trigger info to cluster jobs; fix a bug leading to software stack trigger generating false positives in case of conda environments; fixed display of info message in case of provenance triggered reruns ([#1686](https://www.github.com/snakemake/snakemake/issues/1686)) ([503c70c](https://www.github.com/snakemake/snakemake/commit/503c70c7727e154f8fadf6ead088887d22a87a65))
* set channel priority in container system wide ([#1690](https://www.github.com/snakemake/snakemake/issues/1690)) ([41175b3](https://www.github.com/snakemake/snakemake/commit/41175b3001ff59db5d7ea8bf91354decea04f74d))

## [7.8.0](https://www.github.com/snakemake/snakemake/compare/v7.7.0...v7.8.0) (2022-05-24)


### Features

* automatically rerun jobs if parameters, code, input file set, or software stack changed (thanks to [@cclienty](https://www.github.com/cclienty) and [@timtroendle](https://www.github.com/timtroendle)). This also increases performance of DAG building by handling job "needrun" updates level wise, while avoiding to perform a full check for those jobs that are already downstream of a job that has been determined to require a rerun. ([#1663](https://www.github.com/snakemake/snakemake/issues/1663)) ([4c11893](https://www.github.com/snakemake/snakemake/commit/4c11893d2fda5824adff44d16d7741484e63efea))
* enable the definition of conda pin files in order to freeze an environment. This can drastically increase the robustness because it allows to freeze an environment at a working state. ([#1667](https://www.github.com/snakemake/snakemake/issues/1667)) ([53972bf](https://www.github.com/snakemake/snakemake/commit/53972bfddcca836d5abb8cdd452cbea40ab2571f))


### Bug Fixes

* fail with error if conda installation is not set to strict channel priorities ([#1672](https://www.github.com/snakemake/snakemake/issues/1672)) ([f1ffbf2](https://www.github.com/snakemake/snakemake/commit/f1ffbf28f04150c6e66297f242b768a22f80bd94))
* fix errors occurring when referring to input func via rules.<rulename>.input ([#1669](https://www.github.com/snakemake/snakemake/issues/1669)) ([28a4795](https://www.github.com/snakemake/snakemake/commit/28a47959bca2d135f82f9d8901b2b0aa228f30cb))
* parsing error when combining single line directive with multi-line directive in use rule statements ([#1662](https://www.github.com/snakemake/snakemake/issues/1662)) ([26e57d6](https://www.github.com/snakemake/snakemake/commit/26e57d69fc320adc972967a8046c5163b455456c))

## [7.7.0](https://www.github.com/snakemake/snakemake/compare/v7.6.2...v7.7.0) (2022-05-16)


### Features

* add flag `ensure` that allows to annotate that certain output files should be non-empty or agree with a given checksum ([#1651](https://www.github.com/snakemake/snakemake/issues/1651)) ([76f69d9](https://www.github.com/snakemake/snakemake/commit/76f69d9e21b9c9a9a01198862b66284bc3942d20))
* for small files, compare checksums to determine if job needs to run if input file is newer than output file ([#1568](https://www.github.com/snakemake/snakemake/issues/1568)) ([1ae85c6](https://www.github.com/snakemake/snakemake/commit/1ae85c6b57b3c6b5860214a9c0e3ab28c7c8c5dc))
* LockException ([#1276](https://www.github.com/snakemake/snakemake/issues/1276)) ([f5e6fa6](https://www.github.com/snakemake/snakemake/commit/f5e6fa68a68640ed39df86c256d676f2e4efddfc))
* new directive "retries" for annotating the number of times a job shall be restarted after a failure ([#1649](https://www.github.com/snakemake/snakemake/issues/1649)) ([c8d81d0](https://www.github.com/snakemake/snakemake/commit/c8d81d03de2885d5d0473084141e9f6abc5de445))


### Bug Fixes

* iRODS functionality - issue [#1510](https://www.github.com/snakemake/snakemake/issues/1510) ([#1611](https://www.github.com/snakemake/snakemake/issues/1611)) ([9c3767d](https://www.github.com/snakemake/snakemake/commit/9c3767d6ee13a2c149a4ffe1c0547cabec0346dd))


### Documentation

* singularity sometimes uses system /tmp explanation ([#1588](https://www.github.com/snakemake/snakemake/issues/1588)) ([170c1d9](https://www.github.com/snakemake/snakemake/commit/170c1d9d92de4cafc0da9567a6970b173161c7da))

### [7.6.2](https://www.github.com/snakemake/snakemake/compare/v7.6.1...v7.6.2) (2022-05-06)


### Bug Fixes

* fixed permission issues when using zenodo remote provider to access restricted depositions ([#1634](https://www.github.com/snakemake/snakemake/issues/1634)) ([510f534](https://www.github.com/snakemake/snakemake/commit/510f534ff55635e5c3ca677e0ccd8c5b5dd7ca0f))

### [7.6.1](https://www.github.com/snakemake/snakemake/compare/v7.6.0...v7.6.1) (2022-05-04)


### Bug Fixes

* check for skipped rules in case of local rule inheritance ([#1631](https://www.github.com/snakemake/snakemake/issues/1631)) ([9083ac1](https://www.github.com/snakemake/snakemake/commit/9083ac1f40daf3d284ce9b1ac2d4addde9b5b258))

## [7.6.0](https://www.github.com/snakemake/snakemake/compare/v7.5.0...v7.6.0) (2022-05-03)


### Features

* enable restricted access support in zenodo remote provider ([#1623](https://www.github.com/snakemake/snakemake/issues/1623)) ([692caf9](https://www.github.com/snakemake/snakemake/commit/692caf963d90313d2cd8117fecde097b228633ce))


### Bug Fixes

* avoid erroneous too early deletion of parent directories in case of failed jobs (thanks to @SichongP). ([#1601](https://www.github.com/snakemake/snakemake/issues/1601)) ([b0917e6](https://www.github.com/snakemake/snakemake/commit/b0917e6f07e356764880632495ec3567ec8555b4))
* ensure that rule inheritance considers the same globals and other settings as parent module ([#1621](https://www.github.com/snakemake/snakemake/issues/1621)) ([104cab9](https://www.github.com/snakemake/snakemake/commit/104cab97d9e1a7dda4c9948efa5883d5478d2229))
* issue 1615 - Switch formatting condition for dictionary ([#1617](https://www.github.com/snakemake/snakemake/issues/1617)) ([0771062](https://www.github.com/snakemake/snakemake/commit/0771062a07f0e2cfe9ee45a2276aa61b096eb6e1))
* multiext prefix computation in case it is used within a module that defines an additional prefix ([#1609](https://www.github.com/snakemake/snakemake/issues/1609)) ([fc6dfc6](https://www.github.com/snakemake/snakemake/commit/fc6dfc6469137a82382a36b9469190d967593759))
* remove redundant print ([#1608](https://www.github.com/snakemake/snakemake/issues/1608)) ([cc7e0e3](https://www.github.com/snakemake/snakemake/commit/cc7e0e3605bd65ffcb2d055e69761ff7337588ae))

## [7.5.0](https://www.github.com/snakemake/snakemake/compare/v7.4.0...v7.5.0) (2022-04-26)


### Features

* vim syntax updates ([#1584](https://www.github.com/snakemake/snakemake/issues/1584)) ([b8c77f6](https://www.github.com/snakemake/snakemake/commit/b8c77f6a2a1372a5c3ad8077ad36facf393bfacf))


### Bug Fixes

* properly use configfiles specified via CLI also if configfile specified via configfile directive is not present ([1e0649a](https://www.github.com/snakemake/snakemake/commit/1e0649ac37176a68bb2d8f4d1508ac8bb02463ff))


### Documentation

* checkpoint documentation ([#1562](https://www.github.com/snakemake/snakemake/issues/1562)) ([4cbfb47](https://www.github.com/snakemake/snakemake/commit/4cbfb4786a729a0c899a0a3e0427c1c1f0796c15))

## [7.4.0](https://www.github.com/snakemake/snakemake/compare/v7.3.8...v7.4.0) (2022-04-22)


### Features

* Allow paramspace to separate filename params with custom separator ([#1299](https://www.github.com/snakemake/snakemake/issues/1299)) ([8236e80](https://www.github.com/snakemake/snakemake/commit/8236e80794d0f9c9670238ba168770c0947e8379))


### Bug Fixes

* preserve dtypes across paramspace ([#1578](https://www.github.com/snakemake/snakemake/issues/1578)) ([70ce6a0](https://www.github.com/snakemake/snakemake/commit/70ce6a0feb8572ddcf888c3d377d631ea4a24370))
* use mambaforge for snakemake container image ([#1595](https://www.github.com/snakemake/snakemake/issues/1595)) ([b7e6906](https://www.github.com/snakemake/snakemake/commit/b7e6906926cae5fef6987adcf7b0294266d5faec))

### [7.3.8](https://www.github.com/snakemake/snakemake/compare/v7.3.7...v7.3.8) (2022-04-06)


### Bug Fixes

* support multiple input files for template_engine rules ([#1571](https://www.github.com/snakemake/snakemake/issues/1571)) ([aee7cf2](https://www.github.com/snakemake/snakemake/commit/aee7cf236611e5201feda152f5b7357b49b9f15b))

### [7.3.7](https://www.github.com/snakemake/snakemake/compare/v7.3.6...v7.3.7) (2022-04-05)


### Bug Fixes

* allow labels function to return None ([#1565](https://www.github.com/snakemake/snakemake/issues/1565)) ([fef74d6](https://www.github.com/snakemake/snakemake/commit/fef74d6406a04e29c115a699e76ac96e4a37cf9e))
* do not wrap whitespace in result info headers of reports ([653d0d0](https://www.github.com/snakemake/snakemake/commit/653d0d0b92d2556e0fa04a8208f37fd982dcb829))
* fixed detection of norun rules inside of modules ([#1566](https://www.github.com/snakemake/snakemake/issues/1566)) ([d2223d4](https://www.github.com/snakemake/snakemake/commit/d2223d41dfba057ab735395eac8339c27866c2ae))
* properly use retry mechanism in source cache ([#1564](https://www.github.com/snakemake/snakemake/issues/1564)) ([624a83d](https://www.github.com/snakemake/snakemake/commit/624a83d1bfc592a2a1878d5191e09f6c3d7ee7c2))

### [7.3.6](https://www.github.com/snakemake/snakemake/compare/v7.3.5...v7.3.6) (2022-04-02)


### Bug Fixes

* always recalculate job resources before job is scheduled as input might have changed or not have been present initially ([#1552](https://www.github.com/snakemake/snakemake/issues/1552)) ([44aacdb](https://www.github.com/snakemake/snakemake/commit/44aacdbb35879e1d7914aa105401541465387955))
* fixed handling of input functions and unpack when using the prefix setting of module definitions ([#1553](https://www.github.com/snakemake/snakemake/issues/1553)) ([d561e04](https://www.github.com/snakemake/snakemake/commit/d561e041a8919717046be0d39f197b2c6b937cb7))
* fixed parsing of subsequent use rule statements directly beneath each other ([#1548](https://www.github.com/snakemake/snakemake/issues/1548)) ([77d5a08](https://www.github.com/snakemake/snakemake/commit/77d5a08ab49d67b8cbd8ea4b6b6b7792edb38e3b))
* fix spurious missing file errors when using google storage ([#1541](https://www.github.com/snakemake/snakemake/issues/1541)) ([1b3ede1](https://www.github.com/snakemake/snakemake/commit/1b3ede19159856a982de65e6293ab064c0987352))
* proper error message if resource types do not match ([#1556](https://www.github.com/snakemake/snakemake/issues/1556)) ([1112321](https://www.github.com/snakemake/snakemake/commit/11123213188672a3b6e5acfffed18d9e0ccc8819))
* quote workdir in job exec prefix to allow to spaces in the workdir ([#1547](https://www.github.com/snakemake/snakemake/issues/1547)) ([c3a593e](https://www.github.com/snakemake/snakemake/commit/c3a593e8f7fc8e0dccec4e025f4cd9743bd80bc3))
* report error and possible cause if metadata cleanup fails ([#1554](https://www.github.com/snakemake/snakemake/issues/1554)) ([6866134](https://www.github.com/snakemake/snakemake/commit/68661341efa0a3de4e03de3fb1b8f3117de66efe))

### [7.3.5](https://www.github.com/snakemake/snakemake/compare/v7.3.4...v7.3.5) (2022-03-31)


### Bug Fixes

* do not remove existing temp files in case of dryrun ([#1543](https://www.github.com/snakemake/snakemake/issues/1543)) ([e820f97](https://www.github.com/snakemake/snakemake/commit/e820f973ad8ca99822be69c927c7d7bf6a89f54e))
* fixed bug in missing input file handling for cluster jobs ([#1544](https://www.github.com/snakemake/snakemake/issues/1544)) ([40e2eb2](https://www.github.com/snakemake/snakemake/commit/40e2eb2e6c6e31c7cc590c8d643e0640e9377aa4))


### Documentation

* explain automatic decompression strategy for http remote provider ([e6826b6](https://www.github.com/snakemake/snakemake/commit/e6826b6a740ba5b8877f12732c4ad95194833e07))

### [7.3.4](https://www.github.com/snakemake/snakemake/compare/v7.3.3...v7.3.4) (2022-03-30)


### Bug Fixes

* better error messages in case of missing files after latency period ([#1528](https://www.github.com/snakemake/snakemake/issues/1528)) ([5b394c0](https://www.github.com/snakemake/snakemake/commit/5b394c0319cfc5f8000b616d0f3c911a9091d05b))
* correct handling of exceptions in input functions that are generators ([#1536](https://www.github.com/snakemake/snakemake/issues/1536)) ([d9a56aa](https://www.github.com/snakemake/snakemake/commit/d9a56aaf75c5f70ba0217d9d461d839fa3013f2e))
* obtaining conda prefix when using in combination with singularity ([#1535](https://www.github.com/snakemake/snakemake/issues/1535)) ([99b22d3](https://www.github.com/snakemake/snakemake/commit/99b22d33f2100a7e4cf2d080a2272959a712d055))
* proper error message in case of missing git when checking for source files ([#1534](https://www.github.com/snakemake/snakemake/issues/1534)) ([92887a3](https://www.github.com/snakemake/snakemake/commit/92887a33dc3674f942bd0355edfbba53b810f18f))
* throw error message in case of target rule that depends on a pipe. ([#1532](https://www.github.com/snakemake/snakemake/issues/1532)) ([b9e9a7e](https://www.github.com/snakemake/snakemake/commit/b9e9a7eff4b6e3349dde6b90eec9f5a37ef69ce7))


### Documentation

* display rust-script env. ([950d8ba](https://www.github.com/snakemake/snakemake/commit/950d8ba785a384fa47fcda3d6fb948799a259e0e))
* zenodo example ([76159ae](https://www.github.com/snakemake/snakemake/commit/76159ae22539e38923712e487371a5f32d7cb3cf))

### [7.3.3](https://www.github.com/snakemake/snakemake/compare/v7.3.2...v7.3.3) (2022-03-28)


### Bug Fixes

* better error message in case of failing to create conda env ([#1526](https://www.github.com/snakemake/snakemake/issues/1526)) ([e7a461c](https://www.github.com/snakemake/snakemake/commit/e7a461ce7b5626e603784da27aa4c87649f5edec))
* fix singularity logging messages causing conda fail ([#1523](https://www.github.com/snakemake/snakemake/issues/1523)) ([7797595](https://www.github.com/snakemake/snakemake/commit/77975952ce7df346729f76a767ad4f475b385306))
* more robust handling of incompletely evaluated parameters (any interaction with them will result in a string <TBD> now). ([#1525](https://www.github.com/snakemake/snakemake/issues/1525)) ([3d4c768](https://www.github.com/snakemake/snakemake/commit/3d4c768aafbdca67a9032ad9e3b73449a1fadb0d))


### Documentation

* details on benchmarked results ([64fea09](https://www.github.com/snakemake/snakemake/commit/64fea0921f6a35dbea96435debb114012603ffc2))

### [7.3.2](https://www.github.com/snakemake/snakemake/compare/v7.3.1...v7.3.2) (2022-03-25)


### Bug Fixes

* fixed code change detection ([#1513](https://www.github.com/snakemake/snakemake/issues/1513)) ([67298c6](https://www.github.com/snakemake/snakemake/commit/67298c6167ccaef5f9fbd03ec6b4fe65d86e9ca3))
* modify dag and workflow display in report to also work for big DAGs ([#1517](https://www.github.com/snakemake/snakemake/issues/1517)) ([1364dfb](https://www.github.com/snakemake/snakemake/commit/1364dfbc9db58541aaf5600bca61230f9eb4ecbc))


### Documentation

* Clarify the use of conda with notebook directive ([#1515](https://www.github.com/snakemake/snakemake/issues/1515)) ([aefb1eb](https://www.github.com/snakemake/snakemake/commit/aefb1eb0a2d62faa6108670f3a11d58a1d797c41))

### [7.3.1](https://www.github.com/snakemake/snakemake/compare/v7.3.0...v7.3.1) (2022-03-23)


### Bug Fixes

* add about page to report, including embedded packages and licenses ([#1511](https://www.github.com/snakemake/snakemake/issues/1511)) ([142a452](https://www.github.com/snakemake/snakemake/commit/142a45256f1b192246dd8e9843abedb24badecc6))
* in google live science backend, save multiple logs per rule name and overwrite existing logs ([#1504](https://www.github.com/snakemake/snakemake/issues/1504)) ([9e92d63](https://www.github.com/snakemake/snakemake/commit/9e92d63b9e68b29ccd680c34171994b0a2041efb))
* in rules from imported modules, exclude modified paths from module prefixing ([#1494](https://www.github.com/snakemake/snakemake/issues/1494)) ([1e73db0](https://www.github.com/snakemake/snakemake/commit/1e73db0325f407529108acc689a915ff23611b5a))
* Replaced pathlib relative_to with os.relpath ([#1505](https://www.github.com/snakemake/snakemake/issues/1505)) ([dc65e29](https://www.github.com/snakemake/snakemake/commit/dc65e2921163e9b069c13f79dc0488be21452905))
* update for minimum of Python 3.7 ([#1509](https://www.github.com/snakemake/snakemake/issues/1509)) ([62024e2](https://www.github.com/snakemake/snakemake/commit/62024e2bfd6d5735763f37a0f4bf43a16f229443))

## [7.3.0](https://www.github.com/snakemake/snakemake/compare/v7.2.1...v7.3.0) (2022-03-21)


### Features

* Support for machine_type for kubernetes executor ([#1291](https://www.github.com/snakemake/snakemake/issues/1291)) ([12d6f67](https://www.github.com/snakemake/snakemake/commit/12d6f67a19a55ead092c80b2e5b49e1836fb4e5f))


### Bug Fixes

* always wait for input files before starting jobs, also upon local execution and within group jobs. This should add further robustness against NFS latency issues. ([#1486](https://www.github.com/snakemake/snakemake/issues/1486)) ([cab2adb](https://www.github.com/snakemake/snakemake/commit/cab2adbc2278a2c1689414d2a3f172bb1d5c84d1))
* cleaned up and rewritten execution backend structure, (fixing [#1475](https://www.github.com/snakemake/snakemake/issues/1475), [#860](https://www.github.com/snakemake/snakemake/issues/860), [#1007](https://www.github.com/snakemake/snakemake/issues/1007), [#1008](https://www.github.com/snakemake/snakemake/issues/1008)) (PR [#1491](https://www.github.com/snakemake/snakemake/issues/1491)) ([e87cc97](https://www.github.com/snakemake/snakemake/commit/e87cc979bea0567e1cd97722d385f472857df83c))
* do not skip local conda env creation per se when having no shared FS, because it is still needed for local jobs. Instead, decide for each env whether it is needed locally or not. ([#1490](https://www.github.com/snakemake/snakemake/issues/1490)) ([3f03c5d](https://www.github.com/snakemake/snakemake/commit/3f03c5d303fbdd9e05aa13a4d93bce08cade32b2))
* fixed temp file deletion for group jobs ([#1487](https://www.github.com/snakemake/snakemake/issues/1487)) ([d030443](https://www.github.com/snakemake/snakemake/commit/d030443548a9851a82bcce618b24a9e24a8b545d))
* improve robustness when retrieving remote source files, fixed usage of local git repos as wrapper prefixes (in collaboration with [@cokelaer](https://www.github.com/cokelaer) and @Smeds) ([#1495](https://www.github.com/snakemake/snakemake/issues/1495)) ([e16531d](https://www.github.com/snakemake/snakemake/commit/e16531d6d35b5eb3f7d19008e3e9c4432c4b2e69))
* mtime inventory for google storage was accidentally setting a float instead of a proper mtime object ([#1484](https://www.github.com/snakemake/snakemake/issues/1484)) ([7c762c7](https://www.github.com/snakemake/snakemake/commit/7c762c7e5204f95ca85157ba5fe5ab061b8abdfa))
* render empty caption if nothing defined in report flag ([013a6e8](https://www.github.com/snakemake/snakemake/commit/013a6e8459d0659e05546c849f84151860686004))


### Documentation

* clarify namespacing when using modules. ([dbed4a3](https://www.github.com/snakemake/snakemake/commit/dbed4a3f160106feb15a51d2e8cfcafae531ea57))
* separate api docs ([ded7da9](https://www.github.com/snakemake/snakemake/commit/ded7da90258284f06d4e9263e667cd632cdc12ae))
* separate api docs ([#1499](https://www.github.com/snakemake/snakemake/issues/1499)) ([5cf275a](https://www.github.com/snakemake/snakemake/commit/5cf275ab9c556dd1828a0618799bcdba0c561e70))

### [7.2.1](https://www.github.com/snakemake/snakemake/compare/v7.2.0...v7.2.1) (2022-03-14)


### Bug Fixes

* add missing report.templates.components module to setup.py ([cb4e3fe](https://www.github.com/snakemake/snakemake/commit/cb4e3feaa192ed9c53ddb7c965cb5b71297710c9))


### Documentation

* add install info of development (git) version to docs ([#1477](https://www.github.com/snakemake/snakemake/issues/1477)) ([2a2d6cd](https://www.github.com/snakemake/snakemake/commit/2a2d6cd15cf48278ff17470c7c1323e7ebc40bbd))

## [7.2.0](https://www.github.com/snakemake/snakemake/compare/v7.1.1...v7.2.0) (2022-03-13)


### Features

* improved reports: more interactive and modern interface, ability to define a label based representation of files ([#1470](https://www.github.com/snakemake/snakemake/issues/1470)) ([d09df0c](https://www.github.com/snakemake/snakemake/commit/d09df0c6b02494829345f0af0fa8811007afa28b))


### Bug Fixes

* always deploy conda envs in main process when assuming a shared file system (fixes issue [#1463](https://www.github.com/snakemake/snakemake/issues/1463)) ([#1472](https://www.github.com/snakemake/snakemake/issues/1472)) ([79788eb](https://www.github.com/snakemake/snakemake/commit/79788eb5e8bd404e507ef7e54a0caa6103d90c4e))
* do not wait for named or containerized conda envs ([#1473](https://www.github.com/snakemake/snakemake/issues/1473)) ([6b1d09c](https://www.github.com/snakemake/snakemake/commit/6b1d09c1e270348e8ef77d6ad8c24e1ca540215c))
* implement lock-free source file caching. This avoids hangs on network file systems like NFS. ([#1464](https://www.github.com/snakemake/snakemake/issues/1464)) ([9520e98](https://www.github.com/snakemake/snakemake/commit/9520e988a32f0c5369b4f2c68fdb741f21daa1a4))

### [7.1.1](https://www.github.com/snakemake/snakemake/compare/v7.1.0...v7.1.1) (2022-03-07)


### Bug Fixes

* quote jobid passed to status script to support multi-cluster Slurm setup ([#1459](https://www.github.com/snakemake/snakemake/issues/1459)) ([0232201](https://www.github.com/snakemake/snakemake/commit/023220160c6146810e3da2b277439441e8af9827))

## [7.1.0](https://www.github.com/snakemake/snakemake/compare/v7.0.4...v7.1.0) (2022-03-04)


### Features

* Zenodo remote provider for transparent storage on and retrieval from Zenodo ([#1455](https://www.github.com/snakemake/snakemake/issues/1455)) ([4586ef7](https://www.github.com/snakemake/snakemake/commit/4586ef7c9e5945568e9994a013235574c24d582f))


### Bug Fixes

* disable mtime retrieval from github api for now. This quickly exceeds rate limits. ([1858bb9](https://www.github.com/snakemake/snakemake/commit/1858bb912823da2021f88a9c0cdabe1ee1083575))
* display change warnings only for jobs that won't be executed otherwise ([086f60f](https://www.github.com/snakemake/snakemake/commit/086f60f142721a6085b105bc4bbe12cccc9cee02))
* work around segfault with >100 jobs in google life sciences backend ([#1451](https://www.github.com/snakemake/snakemake/issues/1451)) ([2c0fee2](https://www.github.com/snakemake/snakemake/commit/2c0fee2faec33185ca7fcd2276901977857e2c64))

### [7.0.4](https://www.github.com/snakemake/snakemake/compare/v7.0.3...v7.0.4) (2022-03-03)


### Bug Fixes

* more details on input and output exceptions (missing input, protected output, etc.) ([#1453](https://www.github.com/snakemake/snakemake/issues/1453)) ([8d64af2](https://www.github.com/snakemake/snakemake/commit/8d64af2cb905fef95585055c7b69fd1c45d44108))

### [7.0.3](https://www.github.com/snakemake/snakemake/compare/v7.0.2...v7.0.3) (2022-03-02)


### Bug Fixes

* fix a bug leading to duplicate conda env initializations; fix display of jobs and output files with changes ([994b151](https://www.github.com/snakemake/snakemake/commit/994b1510766083df7f22d10c0e6e4bb65ffdd710))
* preserve empty names input or output file lists in params or resource functions ([0d19ab0](https://www.github.com/snakemake/snakemake/commit/0d19ab0e6fcabe61b49d5ef9f2b293b9bcc06534))
* remove accidental pdb statement ([9c935f1](https://www.github.com/snakemake/snakemake/commit/9c935f1566b976392393aeb00acf0e39eb159e19))
* remove deprecated and add missing arguments to internal functions ([93a7e39](https://www.github.com/snakemake/snakemake/commit/93a7e39d9f225fac5ff5cb8cbe14500a09986ab3))

### [7.0.2](https://www.github.com/snakemake/snakemake/compare/v7.0.1...v7.0.2) (2022-03-01)


### Bug Fixes

* add local marker for input files in cufflinks example. fixes issue [#1362](https://www.github.com/snakemake/snakemake/issues/1362) ([90bc88b](https://www.github.com/snakemake/snakemake/commit/90bc88b84282b477f481e368ad657056e131cbdc))
* failure to properly apply default remote prefix in combination with the unpack marker ([#1448](https://www.github.com/snakemake/snakemake/issues/1448)) ([82666f1](https://www.github.com/snakemake/snakemake/commit/82666f1b2b043f0a8de739d7027aba66eccdaee3))
* set mtime for cached source files [WIP] ([#1443](https://www.github.com/snakemake/snakemake/issues/1443)) ([dd27209](https://www.github.com/snakemake/snakemake/commit/dd27209b4a600d3704cabc39776dfef718129197))
* small bug in snakemake.executors ([#1440](https://www.github.com/snakemake/snakemake/issues/1440)) ([6e64292](https://www.github.com/snakemake/snakemake/commit/6e64292cfa7d5bd9f6cb786681b3710ee51abc43))


### Documentation

* fix list display in docs ([3724367](https://www.github.com/snakemake/snakemake/commit/372436747a97496466b60dc60ee0ebe4cfef1016))
* fix list display in docs ([2dd0e91](https://www.github.com/snakemake/snakemake/commit/2dd0e91f8b7e13d0ffcebe4ed11024a39357ebc7))
* Fix typo and grammar mistake in scatter-gather section. ([#1441](https://www.github.com/snakemake/snakemake/issues/1441)) ([f218aaa](https://www.github.com/snakemake/snakemake/commit/f218aaad1b9b80074ea602cde0352c34c18e70b5))

### [7.0.1](https://www.github.com/snakemake/snakemake/compare/v7.0.0...v7.0.1) (2022-02-26)


### Bug Fixes

* avoid incomplete remote files in case of errors and automatically retry download and upload ([#1432](https://www.github.com/snakemake/snakemake/issues/1432)) ([8fc23ed](https://www.github.com/snakemake/snakemake/commit/8fc23ed09f9c6de7519160797584ff9df3104939))
* do not apply module prefix in case of remote files ([5645b3f](https://www.github.com/snakemake/snakemake/commit/5645b3f75066d8d1d8841b6e6732cd8ad098f67f))
* do not require --cores or --jobs to be set when --cleanup-metadata is used. ([#1429](https://www.github.com/snakemake/snakemake/issues/1429)) ([9c73907](https://www.github.com/snakemake/snakemake/commit/9c739079bf3e6340facaa03f88f757df36f6dd91))
* more robust place for runtime source file cache ([#1436](https://www.github.com/snakemake/snakemake/issues/1436)) ([2681f6f](https://www.github.com/snakemake/snakemake/commit/2681f6f163832dfa5214e10f5234d256f5a13407))
* provide details on error when failing to evaluate default resources ([#1430](https://www.github.com/snakemake/snakemake/issues/1430)) ([04f39a9](https://www.github.com/snakemake/snakemake/commit/04f39a92f58c21265d859666bd63fe686e1d61f5))
* provide proper error when using immediate submit in combination with checkpoint jobs. ([#1437](https://www.github.com/snakemake/snakemake/issues/1437)) ([865cf0f](https://www.github.com/snakemake/snakemake/commit/865cf0f22656e22cd2450e5537421ce70c1705f9))


### Documentation

* explain relative path interpretation ([#1428](https://www.github.com/snakemake/snakemake/issues/1428)) ([add9a05](https://www.github.com/snakemake/snakemake/commit/add9a05eecf10f45113bc511a1d166e1708ff756))
* Fix problems with code blocks and broken internal link. ([#1424](https://www.github.com/snakemake/snakemake/issues/1424)) ([5d4e7d8](https://www.github.com/snakemake/snakemake/commit/5d4e7d8c4d7901c41bfb8f01c4b2c6551add59f7))
* temaplte rendering examples and available variables ([#1431](https://www.github.com/snakemake/snakemake/issues/1431)) ([5995e9e](https://www.github.com/snakemake/snakemake/commit/5995e9ebf7b037479b3f1317cb920773410bd2f2))
* update copyright year ([#1427](https://www.github.com/snakemake/snakemake/issues/1427)) ([6b9f5da](https://www.github.com/snakemake/snakemake/commit/6b9f5da0c986d5de444e00b45656dba85244a6c7))

## [7.0.0](https://www.github.com/snakemake/snakemake/compare/v6.15.5...v7.0.0) (2022-02-23)


### ⚠ BREAKING CHANGES

* require at least Python 3.7 ([fd5daae](https://www.github.com/snakemake/snakemake/commit/fd5daaeff070f9987dba411a0f5262c533a2f666))

### Features

* adding service jobs, i.e. the ability to define jobs that provide a resource for consumers (like a shared memory device or a database), and will be automatically terminated by Snakemake once all consumers are finished. (see [docs](https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#service-rules-jobs), [#1413](https://www.github.com/snakemake/snakemake/issues/1413)) ([a471adb](https://www.github.com/snakemake/snakemake/commit/a471adbb785e5ac7f0c854fd09781c502b577c65))
* support for group local jobs by enabling optional groupid consideration in input functions (see [docs](https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#group-local-jobs), [#1418](https://www.github.com/snakemake/snakemake/issues/1418)) ([5d45493](https://www.github.com/snakemake/snakemake/commit/5d45493db4485af2f4b288b5002605c87315d2b7))
* Adding --cluster-cancel and --cluster-cancel-nargs ([#1395](https://www.github.com/snakemake/snakemake/issues/1395)) ([0593de1](https://www.github.com/snakemake/snakemake/commit/0593de134499712929ba75e65f65df90491eac2e))
* cluster sidecar ([#1397](https://www.github.com/snakemake/snakemake/issues/1397)) ([b992cd1](https://www.github.com/snakemake/snakemake/commit/b992cd19dc1c011f536e3662a3ddffc8b1bb9f67))
* template rendering integration (yte and jinja2) ([#1410](https://www.github.com/snakemake/snakemake/issues/1410)) ([e1cbde5](https://www.github.com/snakemake/snakemake/commit/e1cbde5a378a29e3e7c7c16c73e08b35afa47a56))


### Bug Fixes

* bug in pipe group handling that led to multiple assignments of the same group id to different groups; bug that accidentally added already running groups of the list of ready jobs (issue [#1331](https://www.github.com/snakemake/snakemake/issues/1331)) ([#1332](https://www.github.com/snakemake/snakemake/issues/1332)) ([1a9b483](https://www.github.com/snakemake/snakemake/commit/1a9b483a6c675315d74bff791502c2bdd74609c1))
* display wrapper or external script code in report [#1393](https://www.github.com/snakemake/snakemake/issues/1393) ([#1404](https://www.github.com/snakemake/snakemake/issues/1404)) ([a007bd1](https://www.github.com/snakemake/snakemake/commit/a007bd11fb49a5765a643ec78e19f30b3a10dfab))
* do not pass SNAKEMAKE_PROFILE into cluster-submit ([#1398](https://www.github.com/snakemake/snakemake/issues/1398)) ([#1407](https://www.github.com/snakemake/snakemake/issues/1407)) ([7189183](https://www.github.com/snakemake/snakemake/commit/71891839e397130fb1af3e499d30fa9a953a93f7))
* issue with duplicated prefix for checkpoints on cloud ([#1294](https://www.github.com/snakemake/snakemake/issues/1294)) ([8ed0c8c](https://www.github.com/snakemake/snakemake/commit/8ed0c8cb453b6ebf6df138391a0681ffc8442e09))
* keep flags with apply_wildcards on cloned IOFile ([#1416](https://www.github.com/snakemake/snakemake/issues/1416)) ([23c943f](https://www.github.com/snakemake/snakemake/commit/23c943f0e285f2dc725aa3e4a2e8798021085cb3))
* remove raise that limits using --config with dicts ([#1341](https://www.github.com/snakemake/snakemake/issues/1341)) ([bd65057](https://www.github.com/snakemake/snakemake/commit/bd65057a782355ede86ad0bb912e063ff25a97f5))
* Repair MREs from [#823](https://www.github.com/snakemake/snakemake/issues/823) ([#1203](https://www.github.com/snakemake/snakemake/issues/1203)) ([b007979](https://www.github.com/snakemake/snakemake/commit/b0079791718a390d1f920df15a405cf633314312))
* warn on non-file-modification-date changes like params, code, or input files ([#1419](https://www.github.com/snakemake/snakemake/issues/1419)) ([b5f53f0](https://www.github.com/snakemake/snakemake/commit/b5f53f09ae8c01e1223d2279c3a7f59819a8b44f))

### [6.15.5](https://www.github.com/snakemake/snakemake/compare/v6.15.4...v6.15.5) (2022-02-09)


### Bug Fixes

* convert conda env to string before checks ([#1382](https://www.github.com/snakemake/snakemake/issues/1382)) ([7a8da9f](https://www.github.com/snakemake/snakemake/commit/7a8da9fbf01a037a99ebaa3732fe25e87a96fcd2))
* fix pepfile handling in case of module usage ([#1387](https://www.github.com/snakemake/snakemake/issues/1387)) ([f097a76](https://www.github.com/snakemake/snakemake/commit/f097a761472248d779113cdb22b5274395828bcb))

### [6.15.4](https://www.github.com/snakemake/snakemake/compare/v6.15.3...v6.15.4) (2022-02-09)


### Bug Fixes

* fix issue when generating unit tests for rules with directory output ([#1385](https://www.github.com/snakemake/snakemake/issues/1385)) ([7db614f](https://www.github.com/snakemake/snakemake/commit/7db614fa1753179d2cdc20095df17d5ac2885ad0))


### Documentation

* fix tutorial setup instructions for MacOS. ([#1383](https://www.github.com/snakemake/snakemake/issues/1383)) ([b57b749](https://www.github.com/snakemake/snakemake/commit/b57b7493d372605323204122af859ede38864e4d))

### [6.15.3](https://www.github.com/snakemake/snakemake/compare/v6.15.2...v6.15.3) (2022-02-07)


### Bug Fixes

* skip global report caption when using a module ([#1379](https://www.github.com/snakemake/snakemake/issues/1379)) ([a755cee](https://www.github.com/snakemake/snakemake/commit/a755ceefa478d51070f926beed9090067771edf1))

### [6.15.2](https://www.github.com/snakemake/snakemake/compare/v6.15.1...v6.15.2) (2022-02-05)


### Bug Fixes

* avoid mutable default argument ([#1330](https://www.github.com/snakemake/snakemake/issues/1330)) ([978cc93](https://www.github.com/snakemake/snakemake/commit/978cc9327ce7deb517ad609977e1ce432c58c5e2))
* don't raise WorkflowError when entry is empty ([#1368](https://www.github.com/snakemake/snakemake/issues/1368)) ([1fc6f7b](https://www.github.com/snakemake/snakemake/commit/1fc6f7b5d7e7d7f40baab961db89c4b59c950bf7))
* fix assertion error in conda env file spec when applying wildcards (thanks [@ddesvillechabrol](https://www.github.com/ddesvillechabrol)) ([#1377](https://www.github.com/snakemake/snakemake/issues/1377)) ([6200652](https://www.github.com/snakemake/snakemake/commit/6200652b9aff2362a63581cee58eb9f9cae189da))
* fix None type error when invoking Workflow object manually ([#1366](https://www.github.com/snakemake/snakemake/issues/1366)) ([fca3895](https://www.github.com/snakemake/snakemake/commit/fca3895430c206fc159e71622ee567f77566980d))
* XRootDHelper.exists supports non posix filesystem (object store) ([#1348](https://www.github.com/snakemake/snakemake/issues/1348)) ([7a3ad2f](https://www.github.com/snakemake/snakemake/commit/7a3ad2f438586690dd40e4c8ec591d8c10b22b00))


### Documentation

* add sentence about workflow template to docs ([#1369](https://www.github.com/snakemake/snakemake/issues/1369)) ([5fabffb](https://www.github.com/snakemake/snakemake/commit/5fabffbb4af8e9e122677e5adeaebf2d6bd0eeb3))
* fix typo in installation.rst ([#1344](https://www.github.com/snakemake/snakemake/issues/1344)) ([c45d47a](https://www.github.com/snakemake/snakemake/commit/c45d47a79b78a1afed3b1319e6cafd1b2525fe43))

### [6.15.1](https://www.github.com/snakemake/snakemake/compare/v6.15.0...v6.15.1) (2022-01-31)


### Bug Fixes

* consider post-deploy script for env hashing   ([#1363](https://www.github.com/snakemake/snakemake/issues/1363)) ([d50efd9](https://www.github.com/snakemake/snakemake/commit/d50efd9d16d029fb0e5b14b182882c71a20552bb))

## [6.15.0](https://www.github.com/snakemake/snakemake/compare/v6.14.0...v6.15.0) (2022-01-29)


### Features

* adding default_target directive for declaring default target rules that are not the first rule in the workflow. ([#1358](https://www.github.com/snakemake/snakemake/issues/1358)) ([638ec1a](https://www.github.com/snakemake/snakemake/commit/638ec1a983741cd7ba8faaf1a9dc76ae43d012e5))


### Bug Fixes

* Draft notebook filename with wildcards and params. ([#1352](https://www.github.com/snakemake/snakemake/issues/1352)) ([11d4dc8](https://www.github.com/snakemake/snakemake/commit/11d4dc88598ffb901450bd4e076b91f4e27d37b0))
* proper error message when defining cache eligibility for rules with multiple output files and no multiext declaration. ([#1357](https://www.github.com/snakemake/snakemake/issues/1357)) ([47b5096](https://www.github.com/snakemake/snakemake/commit/47b5096ebbdd3d94a9c99b443064b1b0de389c64))


### Documentation

* Command line arguments for configuration files ([#1343](https://www.github.com/snakemake/snakemake/issues/1343)) ([ad8aaa4](https://www.github.com/snakemake/snakemake/commit/ad8aaa4853a150211513baecc474956575d326eb))
* fix broken link in executor_tutorial/tutorial.rst ([#1360](https://www.github.com/snakemake/snakemake/issues/1360)) ([c9be764](https://www.github.com/snakemake/snakemake/commit/c9be76482d05577c4b1528b0e52ba15fc17a1dd5))

## [6.14.0](https://www.github.com/snakemake/snakemake/compare/v6.13.1...v6.14.0) (2022-01-26)


### Features

* Added timestamp to each log message ([#1304](https://www.github.com/snakemake/snakemake/issues/1304)) ([a5769f0](https://www.github.com/snakemake/snakemake/commit/a5769f0baeaa829b7813dee8c78902edbb42cf4b))
* implement support for removing GFAL remote files ([#1103](https://www.github.com/snakemake/snakemake/issues/1103)) ([25943e5](https://www.github.com/snakemake/snakemake/commit/25943e5630ff6d83afa5cba28edf473ce2ca87da))
* specify conda environments via their name ([#1340](https://www.github.com/snakemake/snakemake/issues/1340)) ([735ab23](https://www.github.com/snakemake/snakemake/commit/735ab2301d0905ea054ad6efa3150acb296d0e78))
* support for post deploy scripts ([#1325](https://www.github.com/snakemake/snakemake/issues/1325)) ([e5dac4f](https://www.github.com/snakemake/snakemake/commit/e5dac4ff297b7aeeb1e1a0bbdd03cb967cee3011))


### Documentation

* link to list of dependencies from installation ([#1336](https://www.github.com/snakemake/snakemake/issues/1336)) ([99d7bfe](https://www.github.com/snakemake/snakemake/commit/99d7bfef1285f131d0e60331511bc4833e7e414a))
* update URL to emacs snakemake-mode ([#1339](https://www.github.com/snakemake/snakemake/issues/1339)) ([dae7b8f](https://www.github.com/snakemake/snakemake/commit/dae7b8fb0e580a1878d36881cfb5ffc8adeaeb9f))

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

-   Fixed a bug occurring when using --resources in the command line
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
-   Fixed metadata migration errors occurring with long output file
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
-   Made hardcoded bucket in google cloud executor configurable.
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
-   Persistence dir "aux" has been renamed to "auxiliary" in order to
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
    not anymore embedded but rather are stored in an auxiliary folder,
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
    partial expansion by skipping wildcards for which no keywords are
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
-   Allow custom log handler definitions via --log-handler-script (e.g.
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
-   Fix ImportError occurring with script directive.
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
    NCBI resource (reference genome, gene/protein sequences, ...) as
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
-   --summary now displays the log files corresponding to each output
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
    This is useful when dealing with filenames that contain
    whitespaces.

#### Changed

-   Snakemake now deletes output files before job execution. Further, it
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
    reusable [wrapper
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
