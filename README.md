# Hashing the human genome: a Bloom filter for all DNA k-mers
In order to detect whether a piece of DNA comes from a specific organism, we want to check how similar it is to the organisms genome. However, this is too computationally demanding. A good way around this problem is to use hashing to construct a Bloom filter for all 30-mers in the organisms genome. A Bloom filter is a data structure designed to find, in a quick and memory-efficient way, whether an element is present in a set (in our case, whether a piece of DNA is part of a genome). It uses map reduction and binary representation to be able to handle large amounts of data.

Once the filter is built, any piece of DNA can be queried against this filter and a probability of it belonging to the organism can be retrieved.... Proper detective work! (In a computationally efficient way that can deal with the enourmous amounts of data that a genome represents :P ) 

In order to run the program, bloom.py must be run on the whole genome and then query.py can be used to query samples.
