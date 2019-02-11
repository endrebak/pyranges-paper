# Jaccard counterexamples:

```
test_jaccard(
gr=+--------------|-----------|-----------|------------|-----------|----------+
| Chromosome   |     Start |       End | Name       |     Score | Strand   |
| (int8)       |   (int64) |   (int64) | (object)   |   (int64) | (int8)   |
|--------------|-----------|-----------|------------|-----------|----------|
| chr1         |         1 |         2 | a          |         0 | +        |
+--------------|-----------|-----------|------------|-----------|----------+
PyRanges object has 1 sequences from 1 chromosomes., gr2=+--------------|-----------|-----------|------------|-----------|----------+
| Chromosome   |     Start |       End | Name       |     Score | Strand   |
| (int8)       |   (int64) |   (int64) | (object)   |   (int64) | (int8)   |
|--------------|-----------|-----------|------------|-----------|----------|
| chr1         |         1 |         4 | a          |         0 | +        |
+--------------|-----------|-----------|------------|-----------|----------+
PyRanges object has 1 sequences from 1 chromosomes., strandedness=False)
Traceback (most recent call last):
  File "/mnt/work/endrebak/software/anaconda/lib/python3.6/site-packages/hypothesis/core.py", line 811, in run
    falsifying_example.__expected_traceback,
  File "/mnt/work/endrebak/software/anaconda/lib/python3.6/site-packages/hypothesis/core.py", line 584, in execute
    result = self.test_runner(data, run)
  File "/mnt/work/endrebak/software/anaconda/lib/python3.6/site-packages/hypothesis/executors.py", line 58, in default_new_style_executor
    return function(data)
  File "/mnt/work/endrebak/software/anaconda/lib/python3.6/site-packages/hypothesis/core.py", line 576, in run
    return test(*args, **kwargs)
  File "/home/endrebak/code/pyranges/tests/test_binary.py", line 234, in test_jaccard
    assert abs(result - bedtools_jaccard) < 0.001
AssertionError: assert 0.166667 < 0.001
 +  where 0.166667 = abs((0.5 - 0.333333))
```

```
You can reproduce this example by temporarily adding @reproduce_failure('3.59.0', b'AXicY2RAA4woPCYGBgAAWwAF') as a decorator on your test case
Falsifying example: test_jaccard(gr=+--------------|-----------|-----------|------------|-----------|----------+
| Chromosome   |     Start |       End | Name       |     Score | Strand   |
| (int8)       |   (int64) |   (int64) | (object)   |   (int64) | (int8)   |
|--------------|-----------|-----------|------------|-----------|----------|
| chr1         |         1 |         2 | a          |         0 | +        |
+--------------|-----------|-----------|------------|-----------|----------+
PyRanges object has 1 sequences from 1 chromosomes., gr2=+--------------|-----------|-----------|------------|-----------|----------+
| Chromosome   |     Start |       End | Name       |     Score | Strand   |
| (int8)       |   (int64) |   (int64) | (object)   |   (int64) | (int8)   |
|--------------|-----------|-----------|------------|-----------|----------|
| chr1         |         2 |         4 | a          |         0 | +        |
+--------------|-----------|-----------|------------|-----------|----------+
PyRanges object has 1 sequences from 1 chromosomes., strandedness=False)
Traceback (most recent call last):
  File "/mnt/work/endrebak/software/anaconda/lib/python3.6/site-packages/hypothesis/core.py", line 811, in run
    falsifying_example.__expected_traceback,
  File "/mnt/work/endrebak/software/anaconda/lib/python3.6/site-packages/hypothesis/core.py", line 584, in execute
    result = self.test_runner(data, run)
  File "/mnt/work/endrebak/software/anaconda/lib/python3.6/site-packages/hypothesis/executors.py", line 58, in default_new_style_executor
    return function(data)
  File "/mnt/work/endrebak/software/anaconda/lib/python3.6/site-packages/hypothesis/core.py", line 576, in run
    return test(*args, **kwargs)
  File "/home/endrebak/code/pyranges/tests/test_binary.py", line 232, in test_jaccard
    assert result == 1
AssertionError: assert 0.0 == 1
```

```
test_jaccard(
gr=+--------------|-----------|-----------|------------|-----------|----------+
| Chromosome   |     Start |       End | Name       |     Score | Strand   |
| (int8)       |   (int64) |   (int64) | (object)   |   (int64) | (int8)   |
|--------------|-----------|-----------|------------|-----------|----------|
| chr1         |         1 |         2 | a          |         0 | +        |
+--------------|-----------|-----------|------------|-----------|----------+
PyRanges object has 1 sequences from 1 chromosomes., gr2=+--------------|-----------|-----------|------------|-----------|----------+
| Chromosome   |     Start |       End | Name       |     Score | Strand   |
| (int8)       |   (int64) |   (int64) | (object)   |   (int64) | (int8)   |
|--------------|-----------|-----------|------------|-----------|----------|
| chr1         |         1 |         3 | a          |         0 | +        |
+--------------|-----------|-----------|------------|-----------|----------+
PyRanges object has 1 sequences from 1 chromosomes., strandedness=False)
```
