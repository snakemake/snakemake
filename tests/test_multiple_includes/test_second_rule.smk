rule test_second_rule: 
    input: rules.test_rule.output
    output: 'test2.txt'
    shell: 'touch {output}'