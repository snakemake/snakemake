rule test_rule: 
    output: 'test1.txt'
    shell: 'touch {output}'