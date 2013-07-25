zsniff
======

zsniff deflate sniffer/extractor

zsniff 0.2 by Candice Quates, Vassil Roussev 07/2013

Usage: ./zsniff [--extract] [--block 65536] filename(s)

Build using a c++11 capable compiler.  

Pass small files or fragments as arguments. 

Files may be broken into blocks if desired using [--block bytes],however
DEFLATE reconstruction may be affected by cut-off structures.   

Decoded data will be extracted if found using [--extract].
