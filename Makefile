SRC=*.d **/*.d

all:
	wc -l $(SRC)
	dmd -unittest -oftest -main -debug -Dddocs $(SRC) && ./test
