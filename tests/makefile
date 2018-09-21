test-lambda:
	cd lambda && fc_run fc_run.cfg
	python ../scripts/get_asm_stats.py lambda/2-asm-falcon lambda/2-asm-falcon
	nosetests test_lambda.py --with-xunit --xunit-file=results/nosetests_lambda.xml
test-ecoli:
	cd ecoli && fc_run fc_run.cfg
	python ../scripts/get_asm_stats.py ecoli/2-asm-falcon ecoli/2-asm-falcon
	nosetests test_ecoli.py --with-xunit --xunit-file=results/nosetests_ecoli.xml
test-yeast:
	cd yeast && fc_run fc_run.cfg
	python ../scripts/get_asm_stats.py yeast/2-asm-falcon yeast/2-asm-falcon
	nosetests test_yeast.py --with-xunit --xunit-file=results/nosetests_yeast.xml

