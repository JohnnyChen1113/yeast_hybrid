# yeast_hybrid
scripts for yeast hybrid data analysis

## list of scripts
- select_bam.py

### Usage of select_bam.py
```
python select_bam.py --help
```
Here is the help information:
```
Usage: select_bam.py [OPTIONS]

  Filter BAM file based on specific criteria: 1. Filter out reads with flag=4
  2. Keep only reads with NH:i:1 or NH:i:2 3. For NH:i:2, check if both
  suffixes have one record each, modify quality if true

Options:
  -i, --input PATH   Input BAM file  [required]
  -o, --output PATH  Output BAM file  [required]
  -s, --suffix TEXT  Two chromosome suffixes separated by comma (e.g.,
                     "_A,_B")  [required]
  --help             Show this message and exit.
```
You can use it like:
```
python select_bam.py --input data/YJF1453_rep2.oak1_N17.bam --output data/YJF1453_rep2.oak1_N17.filtered.bam -s '_HN6,_oak1'
```

We also provide
