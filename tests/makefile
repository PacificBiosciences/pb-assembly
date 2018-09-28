test-ecoli:
	cd ecoli && fc_run fc_run.cfg
	python ../scripts/get_asm_stats.py ecoli/2-asm-falcon
	nosetests test_ecoli.py --with-xunit --xunit-file=results/nosetests_ecoli.xml
	python benchmark.py ecoli
test-yeast:
	cd yeast && fc_run fc_run.cfg
	python ../scripts/get_asm_stats.py yeast/2-asm-falcon
	nosetests test_yeast.py --with-xunit --xunit-file=results/nosetests_yeast.xml
	python benchmark.py yeast
