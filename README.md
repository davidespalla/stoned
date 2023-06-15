# stoned
Software tools for neural data
---

Stoned provides a series of functions for to streamline the  analysis of neural data. The package is under developement.

# Installation

## Prerequisites
- Python 3.7 or above 

## Install from source

We recommend installing `stoned` in a virtual environment. For example, in conda:
```bash
conda create --name stoned_env
conda activate stoned_env
```

Then installing `stoned` with `pip` by cloning the github repository locally:
```bash
git clone https://github.com/davidespalla/stoned.git
cd stoned
python -m pip install .

```

That's it! After following these steps, you should be able to use the functions in `stoned` within the virtual environment.


# Structure

- `lfp` provides functions for signal processing and frequency analysis of local field potential data.
- `spatial` provides function for standard analysis of spatial coding such as ratemaps calculation and skaggs spatial information.