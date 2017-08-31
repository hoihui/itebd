# iTEBD

A **single** header file for itebd calculations based on the [ITensor library](http://itensor.org).

# Installation and Linking

Refer to [ITensor's instruction](http://itensor.org/docs.cgi?page=install) for using iTensor. You just need to add

`#include "itebd.h"`

and make sure the header file is within path.

# Examples

I have created a few examples to reproduce figures of some arXiv papers that utilized iTEBD:

* [arXiv:cond-mat/0605597](https://arxiv.org/abs/cond-mat/0605597) (figs. 5 and 6)

* [arXiv:0907.3206](https://arxiv.org/abs/0907.3206) (fig.3b)

* [arXiv:1302.4460](https://arxiv.org/abs/1302.4460) (fig.2a)

* [arXiv:1503.02010](https://arxiv.org/abs/1503.02010) (fig.1)

* [arXiv:1604.03571](https://arxiv.org/abs/1604.03571) (fig.2a)


Please `make` the corresponding file's name you want to compile. For example, `make 1503.02010` generates executable `1503.02010` that produces figure 1 of [arXiv:1503.02010](https://arxiv.org/abs/1503.02010).

By default, it is now set to reproduce figures 5 and 6 of [Vidal's original paper](https://arxiv.org/abs/cond-mat/0605597).


# Basic Usage

Henceforth I assume we have `using namespace itensor;`

## Instantiation

```C++
itebd<IndexType,z> sys(ampo,QNA,QNB,Chi0);
```

* `IndexType` could be `Index` or `IQIndex`, depending on whether you want to use symmetries.

* `z` optional and defaults to 2. The coordination number. Could be useful for exotic system with >2 neighbors e.g. [arXiv:0712.1806](https://arxiv.org/abs/0712.1806) (I have not tried).

* `ampo` is a two-site [AutoMPO](http://itensor.org/docs.cgi?page=classes/autompo) object of ITensor. Hamiltonian is specified by `ampo`. The [sites](http://itensor.org/docs.cgi?page=classes/siteset) must be of type `IQIndex`.

* `QNA` and `QNB` are optional and defaults to 1. The initial quantum numbers of the two sites (site A and site B). Physical meanings of these numbers are determined by the sites which are IQIndex .

* `Chi0` is optional and defaults to 1. The initial bond dimension **of each set of quantum numbers**.

## Propagation (step)

```C++
double itebd::step(std::complex<double> dt, size_t steps = 1, double thres = 1E-10, int maxm = 0);
```

* `dt` is the timestep **on real axis**. When it is real, it propagates in real time.

* `steps` number of steps

* `thres` is the `Cutoff` value used for [SVD](http://itensor.org/docs.cgi?page=classes/svdalgs).

* `maxm` is the `Maxm` value used for [SVD](http://itensor.org/docs.cgi?page=classes/svdalgs).

* **returns**: Energy estimated based on the log of the expection value of the evolution operator


```C++
double itebd::step_imag(std::complex<double> dt, size_t steps = 1, double thres = 1E-10, int maxm = 0);
```

* identical to above but in the imaginary axis.


## Local Measurements

```C++
std::complex<double> itebd::measure(const T &op)
```

*  `op` is a two-neighboring sites operator with indices the same as the site used to construct `ampo`.

* **returns**: Expectation value of the operator at the current state

## Long-Range two-opeartor correlation measurements

```C++
std::vector<std::complex<double>> itebd::measure(const T &op1A, const T &op2, const std::vector<int> &rv);
```

* `op1A` is the operator at the origin and must have the same Index as the first site in the sites used to contruct `ampo`

* `opt2` is the opeartor at the second site. Can have either one of the site indices.

* `rv` is the list of integer distances between sites 1 and 2 to compute <opt1A*opt2>

* **returns**: a vector of complex values of <opt1A*opt2>

## Helpful member variables / functions

* `.tv` (vector of doubles): the times (in absolute values) traversed by the system since it starts

* `.Sv` (vector of doubles): the entanglement entropies the system had for each time in `.tv`

* `.bonddim` (int): current bond dimension

* `.GA`, `.GB`, `.L[0]`, `.L[1]`: the Gammas and Lambdas of the state as used in [Vidal's paper](https://arxiv.org/abs/cond-mat/0605597)

* `.resetTime()` (void): resets both `.tv` and `.Sv`

* `.setH(const AutoMPO &ampo)` (void): sets a new Hamiltonian while keeping the current states unchanged.

* `.backup()` and `.restore()` (void): backup and restore the state

*  `.randomize()` (void): randomize all the tensors (e.g. before imaginary time propagation)

# Higher Dimensions?

By instantiating with `z=4` you can try to use this for 2-dimensional systems where the "environment" is never computed. Some call this "[simplified update](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.81.165104)" of iPEPS and saw good results **if the environment is at least computed for measurements**. I have not implemented this yet so it should never be used for serious work in 2D. cond-mat.0605597 implements simplified update without calculation of environment.
