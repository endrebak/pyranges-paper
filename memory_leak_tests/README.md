# Memory Leak Tests

`valgrind --log-file=valgrind_pandas.txt --tool=memcheck --leak-check=full --num-callers=30 --suppressions=tests/valgrind-python.supp python tests/test_pandas.py`
