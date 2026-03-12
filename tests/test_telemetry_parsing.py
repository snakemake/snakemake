"""Unit tests for telemetry parsing utilities."""

import pytest
from snakemake.telemetry._parsing import get_file_type, parse_shell_tool


class TestGetFileType:
    """Test file type extraction from file paths."""

    # Standard cases
    def test_simple_extension(self):
        """Test basic file extension extraction."""
        assert get_file_type("data.txt") == ".txt"
        assert get_file_type("/path/to/file.csv") == ".csv"

    def test_compound_gz_extension(self):
        """Test compound .fastq.gz style extensions."""
        assert get_file_type("reads.fastq.gz") == ".fastq.gz"
        assert get_file_type("/data/sample.fastq.gz") == ".fastq.gz"

    def test_compound_bz2_extension(self):
        """Test compound .tar.bz2 style extensions."""
        assert get_file_type("archive.tar.bz2") == ".tar.bz2"

    def test_compound_xz_extension(self):
        """Test compound .tar.xz style extensions."""
        assert get_file_type("archive.tar.xz") == ".tar.xz"

    def test_compound_zst_extension(self):
        """Test compound .tar.zst style extensions."""
        assert get_file_type("archive.tar.zst") == ".tar.zst"

    def test_compound_zip_extension(self):
        """Test compound .zip extensions."""
        assert get_file_type("archive.zip") == ".zip"
        assert get_file_type("data.tar.zip") == ".tar.zip"

    def test_compound_other_compression(self):
        """Test other compression suffixes."""
        assert get_file_type("archive.tar.rar") == ".tar.rar"
        assert get_file_type("archive.tar.7z") == ".tar.7z"

    def test_no_extension(self):
        """Test files without extensions."""
        assert get_file_type("Makefile") is None
        assert get_file_type("README") is None
        assert get_file_type("/path/to/bin/noext") is None

    def test_empty_path(self):
        """Test empty string handling."""
        assert get_file_type("") is None

    def test_dotfile(self):
        """Test dotfiles (hidden files starting with .)."""
        # Dotfiles like .bashrc - the function treats the part after the first dot
        # as the extension. This is a design choice - could be debated either way.
        result = get_file_type(".bashrc")
        # Currently returns ".bashrc" (everything after the first dot)
        # Alternative interpretation: None (dotfiles have no extension)
        assert result in [".bashrc", None]

        # .gitignore is similar
        result = get_file_type(".gitignore")
        assert result in [".gitignore", None]

    def test_multiple_dots(self):
        """Test files with multiple dots but not compound extensions."""
        # file.name.with.dots.txt -> should return .txt
        assert get_file_type("file.name.with.dots.txt") == ".txt"

    # Edge cases and brittleness issues

    def test_space_separated_paths_returns_first_file_ext(self):
        """Test handling of space-separated paths (best-effort fallback).

        When multiple file paths are joined with spaces (as happens when
        Namedlist is accidentally converted to string), the function now
        extracts the extension from the first file only as a best-effort
        fallback, rather than returning garbage.
        """
        # Simulating what happens when Namedlist.__str__() is called:
        # InputFiles(["genome.fa", "genome.amb", ...]) becomes:
        # "genome.fa genome.amb genome.ann genome.pac genome.0123 genome.64"
        space_separated = (
            "genome.fa genome.amb genome.ann genome.pac genome.0123 genome.64"
        )
        result = get_file_type(space_separated)
        # Should extract extension from first file only
        assert result == ".fa"

    def test_space_separated_with_fastq_gz(self):
        """Test space-separated paths with compound extension."""
        space_sep = "reads.fastq.gz reads2.fastq.gz"
        result = get_file_type(space_sep)
        # Should extract .fastq.gz from the first file
        assert result == ".fastq.gz"

    def test_single_file_with_numeric_extension(self):
        """Test numeric file extensions like .0123, .64."""
        # These are legitimate extensions used by bioinformatics tools (e.g., BWA index files)
        assert get_file_type("genome.0123") == ".0123"
        assert get_file_type("genome.64") == ".64"

    def test_path_with_spaces(self):
        """Test paths containing spaces in directory names."""
        # This is different from space-separated file lists
        # A single file with spaces in the path
        assert get_file_type("/path with spaces/file.txt") == ".txt"
        assert get_file_type("/my data/results.csv") == ".csv"

    def test_weird_extensions(self):
        """Test unusual but valid extensions."""
        assert get_file_type("file.a") == ".a"  # single char extension
        assert get_file_type("file.AB") == ".AB"  # uppercase extension
        assert get_file_type("file.ext1.ext2") == ".ext2"  # not a compound extension

    def test_url_like_paths(self):
        """Test paths that look like URLs."""
        # These might appear in remote storage contexts
        assert get_file_type("s3://bucket/file.txt") == ".txt"
        assert get_file_type("https://example.com/data.csv") == ".csv"

    def test_trailing_slash(self):
        """Test paths with trailing slashes."""
        # Directories - should handle gracefully
        assert get_file_type("/path/to/dir/") is None
        assert get_file_type("dir/") is None


class TestParseShellTool:
    """Test shell command parsing."""

    def test_simple_command(self):
        """Test basic command extraction."""
        tool, params = parse_shell_tool("samtools sort -@ 4 input.bam")
        assert tool == "samtools"
        assert params == "sort -@ 4 input.bam"

    def test_command_with_interpreter(self):
        """Test that interpreters like python are skipped."""
        tool, params = parse_shell_tool("python scripts/process.py --input data")
        assert tool == "scripts/process.py"
        assert "--input" in params

    def test_command_with_python3(self):
        """Test python3 interpreter skipping."""
        tool, params = parse_shell_tool("python3 script.py arg1 arg2")
        assert tool == "script.py"
        assert params == "arg1 arg2"

    def test_empty_command(self):
        """Test empty string handling."""
        tool, params = parse_shell_tool("")
        assert tool is None
        assert params == ""

    def test_comment_only(self):
        """Test command with only comments."""
        tool, params = parse_shell_tool("# This is a comment")
        assert tool is None
        assert params == ""

    def test_multiline_with_comments(self):
        """Test multiline command with comments."""
        cmd = """# Comment line
        samtools view -b input.sam"""
        tool, params = parse_shell_tool(cmd)
        assert tool == "samtools"
        assert "view" in params

    def test_only_interpreter(self):
        """Test command that is only an interpreter."""
        tool, params = parse_shell_tool("python3")
        assert tool is None
        assert params == ""

    def test_bash_script(self):
        """Test bash command parsing."""
        tool, params = parse_shell_tool("bash script.sh arg1")
        assert tool == "script.sh"
        assert params == "arg1"

    def test_multiline_complex(self):
        """Test complex multiline shell command."""
        cmd = """# Index the genome
bwa index {input}

# Map reads
bwa mem {input} {output}"""
        tool, params = parse_shell_tool(cmd)
        assert tool == "bwa"
        assert params == "index {input}"
