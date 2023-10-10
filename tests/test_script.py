from textwrap import dedent

from snakemake.io import InputFiles
from snakemake.script import RustScript, BashEncoder


class TestRustScriptExtractManifest:
    def test_single_line_manifest_with_shebang_and_second_manifest(self):
        source = dedent(
            """#!/usr/bin/env rust-script
// cargo-deps: time="0.1.25", serde="*"
// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = '// cargo-deps: time="0.1.25", serde="*"\n'

        assert manifest == expected_manifest

        expected_remaining_src = dedent(
            """// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        assert remaining_src == expected_remaining_src

    def test_single_line_manifest_not_at_start_with_shebang(self):
        source = dedent(
            """#!/usr/bin/env rust-script
// this is where cargo-deps should be
// cargo-deps: time="0.1.25", serde="*"
// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = ""

        assert manifest == expected_manifest

        expected_remaining_src = dedent(
            """// this is where cargo-deps should be
// cargo-deps: time="0.1.25", serde="*"
// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        assert remaining_src == expected_remaining_src

    def test_single_line_manifest_not_at_start_without_shebang(self):
        source = dedent(
            """// this is where cargo-deps should be
// cargo-deps: time="0.1.25", serde="*"
// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = ""

        assert manifest == expected_manifest

        expected_remaining_src = dedent(
            """// this is where cargo-deps should be
// cargo-deps: time="0.1.25", serde="*"
// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        assert remaining_src == expected_remaining_src

    def test_single_line_manifest_with_empty_line_without_shebang(self):
        source = dedent(
            """
// cargo-deps: time="0.1.25", serde="*"
// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = '\n// cargo-deps: time="0.1.25", serde="*"\n'

        assert manifest == expected_manifest

        expected_remaining_src = dedent(
            """// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        assert remaining_src == expected_remaining_src

    def test_single_line_manifest_is_case_insensitive(self):
        source = dedent(
            """
// Cargo-deps: time="0.1.25", serde="*"
// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = '\n// Cargo-deps: time="0.1.25", serde="*"\n'

        assert manifest == expected_manifest

        expected_remaining_src = dedent(
            """// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        assert remaining_src == expected_remaining_src

    def test_single_line_manifest_spacing_has_no_impact(self):
        source = dedent(
            """
// cargo-deps : time="0.1.25", serde="*"
// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = '\n// cargo-deps : time="0.1.25", serde="*"\n'

        assert manifest == expected_manifest

        expected_remaining_src = dedent(
            """// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        assert remaining_src == expected_remaining_src

    def test_single_line_manifest_formatting_not_touched_even_if_wrong(self):
        """The dependency delimiter is wrong, but we let rust-script deal with it"""
        source = dedent(
            """
// cargo-deps: time="0.1.25"; serde="*"
// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = '\n// cargo-deps: time="0.1.25"; serde="*"\n'

        assert manifest == expected_manifest

        expected_remaining_src = dedent(
            """// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        assert remaining_src == expected_remaining_src

    def test_single_line_manifest_spelt_wrong(self):
        source = dedent(
            """
// cargo-dependencies: time="0.1.25"; serde="*"
// You can also leave off the version number, in which case, it's assumed
// to be "*".  Also, the `cargo-deps` comment *must* be a single-line
// comment, and it *must* be the first thing in the file, after the
// shebang.
// This second dependency line should be ignored
// cargo-deps: time="0.1.25", libc="0.2.5"
fn main() {
    println!("{}", time::now().rfc822z());
}
"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = ""

        assert manifest == expected_manifest

        expected_remaining_src = source

        assert remaining_src == expected_remaining_src

    def test_code_block_manifest_with_shebang(self):
        source = dedent(
            """#!/usr/bin/env rust-script
//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! [dependencies]
//! time = "0.1.25"
//! ```
fn main() {
    println!("{}", time::now().rfc822z());
}

"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = dedent(
            """//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! [dependencies]
//! time = "0.1.25"
//! ```
"""
        )

        assert manifest == expected_manifest

        expected_remaining_src = dedent(
            """fn main() {
    println!("{}", time::now().rfc822z());
}

"""
        )

        assert remaining_src == expected_remaining_src

    def test_code_block_manifest_without_shebang(self):
        source = dedent(
            """
//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! [dependencies]
//! time = "0.1.25"
//! ```
fn main() {
    println!("{}", time::now().rfc822z());
}

"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = dedent(
            """\n//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! [dependencies]
//! time = "0.1.25"
//! ```
"""
        )

        assert manifest == expected_manifest

        expected_remaining_src = dedent(
            """fn main() {
    println!("{}", time::now().rfc822z());
}

"""
        )

        assert remaining_src == expected_remaining_src

    def test_code_block_manifest_spacing_around_language(self):
        source = dedent(
            """
//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```  cargo
//! [dependencies]
//! time = "0.1.25"
//! ```
fn main() {
    println!("{}", time::now().rfc822z());
}

"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = dedent(
            """\n//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```  cargo
//! [dependencies]
//! time = "0.1.25"
//! ```
"""
        )

        assert manifest == expected_manifest

        expected_remaining_src = dedent(
            """fn main() {
    println!("{}", time::now().rfc822z());
}

"""
        )

        assert remaining_src == expected_remaining_src

    def test_code_block_manifest_has_non_cargo_block(self):
        source = dedent(
            """
//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```rust
//! [dependencies]
//! time = "0.1.25"
//! ```
fn main() {
    println!("{}", time::now().rfc822z());
}

"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = ""

        assert manifest == expected_manifest

        expected_remaining_src = source

        assert remaining_src == expected_remaining_src

    def test_code_block_manifest_missing_closing_fence(self):
        source = dedent(
            """
//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! [dependencies]
//! time = "0.1.25"
//! 
fn main() {
    println!("{}", time::now().rfc822z());
}

"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = ""

        assert manifest == expected_manifest

        expected_remaining_src = source

        assert remaining_src == expected_remaining_src

    def test_code_block_manifest_not_in_first_comment_block(self):
        source = dedent(
            """//! crate comment
static FOO: &str = "foo";
//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! [dependencies]
//! time = "0.1.25"
//! ```
//! 
fn main() {
    println!("{}", time::now().rfc822z());
}

"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = ""

        assert manifest == expected_manifest

        expected_remaining_src = source

        assert remaining_src == expected_remaining_src

    def test_code_block_manifest_with_outer_line_doc_comment(self):
        source = dedent(
            """#!/usr/bin/env rust-script
/// This is a regular crate doc comment, but it also contains a partial
/// Cargo manifest.  Note the use of a *fenced* code block, and the
/// `cargo` "language".
///
/// ```cargo
/// [dependencies]
/// time = "0.1.25"
/// ```
fn main() {
    println!("{}", time::now().rfc822z());
}

"""
        )

        manifest, remaining_src = RustScript.extract_manifest(source)

        expected_manifest = dedent(
            """/// This is a regular crate doc comment, but it also contains a partial
/// Cargo manifest.  Note the use of a *fenced* code block, and the
/// `cargo` "language".
///
/// ```cargo
/// [dependencies]
/// time = "0.1.25"
/// ```
"""
        )

        assert manifest == expected_manifest

        expected_remaining_src = dedent(
            """fn main() {
    println!("{}", time::now().rfc822z());
}

"""
        )

        assert remaining_src == expected_remaining_src


class TestBashEncoder:
    def test_named_list_one_named_one_str(self):
        """InputFiles is a subclass of snakemake.io.NamedInput
        ierate over input and store each with the integer index - i.e 0, 1, 2
        then use input.items() to iterate over the named files and store them as named also
        check how this works with named things being lists
        """
        named_list = InputFiles(["test.in", "named.in"])
        named_list._set_name("named", 1)

        actual = BashEncoder.encode_namedlist(named_list)
        expected = r"""( [0]="test.in" [1]="named.in" [named]="named.in" )"""

        assert actual == expected

    def test_named_list_named_is_list(self):
        """Named lists that are lists of files become a space-separated string as you
        can't nest arrays in bash"""
        named_list = InputFiles(["test1.in", ["test2.in", "named.in"]])
        named_list._set_name("named", 1)

        actual = BashEncoder.encode_namedlist(named_list)
        expected = r"""( [0]="test1.in" [1]="test2.in named.in" [named]="test2.in named.in" )"""

        assert actual == expected
