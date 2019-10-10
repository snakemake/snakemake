from snakemake import get_argument_parser, available_cpu_count


class TestArgumentParser:
    def test_no_cores_passed_sets_cores_to_default(self):
        parser = get_argument_parser()
        cli_args = ["snakemake"]
        args = parser.parse_args(cli_args)

        actual = args.cores
        expected = 1

        assert actual == expected

    def test_no_cores_number_passed_sets_cores_to_number_available(self):
        parser = get_argument_parser()
        cli_args = ["snakemake", "--cores"]
        args = parser.parse_args(cli_args)

        actual = args.cores
        expected = available_cpu_count()

        assert actual == expected

    def test_zero_cores_passed_sets_cores_to_zero(self):
        parser = get_argument_parser()
        cli_args = ["snakemake", "--jobs", "0"]
        args = parser.parse_args(cli_args)

        actual = args.cores
        expected = 0

        assert actual == expected

    def test_two_cores_passed_sets_cores_to_two(self):
        parser = get_argument_parser()
        cli_args = ["snakemake", "--jobs", "2"]
        args = parser.parse_args(cli_args)

        actual = args.cores
        expected = 2

        assert actual == expected
