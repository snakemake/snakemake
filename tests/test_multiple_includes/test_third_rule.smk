rule test_third_rule:
    input: rules.test_second_rule.output
    output: 'test3.txt'
    shell: 'touch {output}'