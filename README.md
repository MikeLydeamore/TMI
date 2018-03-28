# TMI: Two-state Markov Model Inference on interval-censored data

This package performs inference on interval-censored (panel) data, on a two-state Markov Chain.

To install, use the following snippet:

```
library(devtools)
install_github("MikeLydeamore/TMI")
```

This package requires a C++ compiler to install.
Also required is to compile with C++11 on. To do this, please execute:
```
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
```

Please note that package load time is quite long, due to compilation of the stan model.

To run the chains in parallel (highly recommended)
```
options(mc.cores = parallel::detectCores())
```
