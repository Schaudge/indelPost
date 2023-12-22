# indelPost

[![Documentation Status](https://readthedocs.org/projects/indelpost/badge/?version=latest)](https://indelpost.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/indelpost.svg)](https://badge.fury.io/py/indelpost)

indelPost is a Python library for indel processing via realignment and read-based phasing to resolve alignment ambiguities. By importing the library, 
users write their own scripts to solve alignment-sensitive problems such as:
* compare/integrate indels that are differently called by multiple variant callers (e.g., complex indels).
* compare indel alignments in multiple mappings (e.g., match DNA indels to RNA-Seq for expression check).  
* construct a complex indel from a simple indel by read-based phasing.    
* extract reads supporting the target indel from BAM file.
* pull variant records matching the target indel from VCF file.
* count reads supporting indels by realignment.

Visit [documentation](https://indelpost.readthedocs.io/en/latest) for detail.

To install (require Linux with Python>=3.6 pre-installed):
```
pip3 install indelpost
```
or from source for develop or debug
```
git clone https://github.com/stjude/indelPost.git
cd indelPost 
python3 -m pip install --prefix=/home/schaudge/.local/lib/python3.10/site-packages/ --editable .  # recommand for debug mode
```
or
```
pip3 install --prefix=/home/schaudge/.local/lib/python3.10/site-packages/ --editable .
```

## Troubleshoot
If you get something like:
```
... may indicate binary incompatibility. Expected 88 from C header, got 72 from PyObject
```
or
```
AttributeError: module 'pysam.libcalignmentfile' has no attribute 'IteratorColumnAll'
```
try:
```
pip3 uninstall cython pysam indelpost
pip3 install cython
pip3 install pysam
pip3 install indelpost --no-binary indelpost --no-build-isolation
```

## Debug
Make sure the python3-dbg and cygdb were installed! And recompile the cython pyx file by python3-gdb (add gdb_debug=True parameter in cythonize!)
```
python3-dbg setup.py build_ext --inplace
```
then, begin the debug process (with source code script annotate_indels.py):
```
cygdb . -- --args python3-dbg annotate_indels.py
```
A text 
```
Type "apropos cy" to search for commands related to "cy"...         # cy ---> cython debug command
Reading symbols from python3-dbg...
(gdb) cy break count_alleles       # function name, other format: indelpost.pileup.fetch_reads, indelpost.varaln.VariantAlignment.fetch_reads
(gdb) cy run
(gdb) cy list
(gdb) cy locals
(gdb) cy next (step)
(gdb) 
```
some debug tips for review:

cy break indelpost.varaln:610   # some gdb warning tips was not right!
cy break indelpost.varaln.VariantAlignment.count_alleles

If    "cy print variable"    will be broken, try:   "cy exec print(variable)"   for a replacement.


## Reference
Hagiwara K et al. (2022) indelPost: harmonizing ambiguities in simple and complex indel alignments. [Bioinformatics](https://doi.org/10.1093/bioinformatics/btab601)

## Contact
* kohei.hagiwara[AT]stjude.org 
