rule:
	input: 'test.in'
	output: 'test.inter'
	shell: 
		'cp {input} {output}'
