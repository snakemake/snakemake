import sys
import os
from os.path import join
from subprocess import call
from tempfile import mkdtemp
import hashlib
import urllib
from snakemake import snakemake

__author__ = "Tobias Marschall, Marcel Martin"

def dpath(path):
	"""get path to a data file (relative to the directory this
	test lives in)"""
	return join(os.path.dirname(__file__), path)

SCRIPTPATH = dpath("../bin/snakemake")

def md5sum(filename):
	data = open(filename, 'rb').read()
	return hashlib.md5(data).hexdigest()


def is_connected():
    try:
        urllib.request.urlopen("http://www.google.com", timeout=1)
        return True
    except urllib.request.URLError:
        return False


def run(path, shouldfail=False, needs_connection=False, snakefile="Snakefile", **params):
    """
    Test the Snakefile in path.
    There must be a Snakefile in the path and a subdirectory named
    expected-results.
    """
    if needs_connection and not is_connected():
        print("Skipping test because of missing internet connection", file=sys.stderr)
        return False

    results_dir = join(path, 'expected-results')
    snakefile = join(path, snakefile)
    assert os.path.exists(snakefile)
    assert os.path.exists(results_dir) and os.path.isdir(results_dir), '{} does not exist'.format(results_dir)
    tmpdir = mkdtemp()
    try:
        call('cp `find {} -maxdepth 1 -type f` {}'.format(path, tmpdir), shell=True)
        success = snakemake(snakefile, cores=3, workdir=tmpdir, stats = "stats.txt", snakemakepath = SCRIPTPATH, **params)
        if shouldfail:
            assert not success, "expected error on execution"
        else:
            assert success, "expected successful execution"
            for resultfile in os.listdir(results_dir):
                if not os.path.isfile(os.path.join(results_dir, resultfile)):
                    # FIXME: this means tests cannot use directories
                    continue # skip .svn dirs etc.
                targetfile = join(tmpdir, resultfile)
                expectedfile = join(results_dir, resultfile)
                assert os.path.exists(targetfile), 'expected file "{}" not produced'.format(resultfile)
                assert md5sum(targetfile) == md5sum(expectedfile), 'wrong result produced for file "{}"'.format(resultfile)
    finally:
        call(['rm', '-rf', tmpdir])


def test01():
    run(dpath("test01"))

def test02():
    run(dpath("test02"))

def test03():
    run(dpath("test03"), targets=['test.out'])

def test04():
    run(dpath("test04"), targets=['test.out'])

def test05():
    run(dpath("test05"))

def test06():
    run(dpath("test06"), targets=['test.bla.out'])

def test07():
    run(dpath("test07"), targets=['test.out', 'test2.out'])

def test08():
    run(dpath("test08"), targets=['test.out', 'test2.out'])

def test09():
    run(dpath("test09"), shouldfail=True)

def test10():
    run(dpath("test10"))

def test11():
    run(dpath("test11"))

def test12():
    run(dpath("test12"))

def test13():
    run(dpath("test13"))

def test14():
    run(dpath("test14"), snakefile="Snakefile.nonstandard", cluster="./qsub")

def test15():
    run(dpath("test15"))

def test16():
    run(dpath("test16"))

def test_dynamic():
    run(dpath("test_dynamic"))

def test_params():
    run(dpath("test_params"))

def test_same_wildcard():
    run(dpath("test_same_wildcard"))

def test_conditional():
    run(dpath("test_conditional"), targets="test.out test.0.out test.1.out test.2.out".split())

def test_shell():
    run(dpath("test_shell"))

def test_temp():
    run(dpath("test_temp"), cluster="./qsub", targets="test.realigned.bam".split())

def test_keyword_list():
    run(dpath("test_keyword_list"))

def test_subworkflows():
    run(dpath("test_subworkflows"))

def test_globwildcards():
    run(dpath("test_globwildcards"))

def test_local_import():
    run(dpath("test_local_import"))

def test_ruledeps():
    run(dpath("test_ruledeps"))

def test_persistent_dict():
    run(dpath("test_persistent_dict"))

def test_url_include():
    run(dpath("test_url_include"), needs_connection=True)

def test_touch():
    run(dpath("test_touch"))

def test_config():
    run(dpath("test_config"))

def test_benchmark():
    run(dpath("test_benchmark"))

def test_temp_expand():
    run(dpath("test_temp_expand"))

def test_cluster_dynamic():
    run(dpath("test_cluster_dynamic"), cluster="./qsub")
