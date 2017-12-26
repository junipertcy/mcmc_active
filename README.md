# mcmc_active [![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com) [![Build Status](https://travis-ci.org/junipertcy/mcmc_active.svg?branch=master)](https://travis-ci.org/junipertcy/mcmc_active)

**mcmc_active** implements the MCMC algorithm to actively learn the node labels in a stochastic block model in an optimal manner.

This program is an object-oriented re-implemetation of the code by Yaojia Zhu. Please refer to the [original release page](https://github.com/everyxs/gephi-plugins/releases/tag/%40active).

Documentation will be updated soon.

## Usage
Depends on `boost::program_options` and `cmake`.

Compilation:
```
cmake .
make
```
The binaries are built in `bin/`.



## Example call and output
If the build was successful, we can run it via:
```
bin/mcmc_active -i dataset/karate.gml -r 3 -g 2 -p 5 -q 1
```

the output should look like:

```
gml_path: dataset/karate.gml
Run #1 (out of 3) is running...
0.685983,0.595556,0.604112,0.448664, ...,
The indexes of the selected node(s) are: 14,
...
Run #2 (out of 3) is running...
...

```

where `0.685983,0.595556,0.604112,0.448664, ...,` is in std::cout whilst the others appear in std::clog (thereby allowing for easy redirection of the output).
An float on the output line corresponds to the expected mutual information gain of the vertex `v_0, v_1,..., v_n`, assuming the correct label of which vertex is known and fixed. 


