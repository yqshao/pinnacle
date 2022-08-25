# Biased Sampling Methods

## Ensemble-biased Sampling

The `tips.bias.EnsembleBiasedCalculator` implements an ensemble-based sampling
methods that applies a biasing potential.

### Usage

The `EnsembleBiasedCalculator` is implemented in ASE and can be used like any
other calculators in ASE to run e.g. MD simulations. The calculator computes
several other properties named as:

- `*_bias`: the biasing potential or the corresponding force
- `*_avg`: the average potential or forces
- `*_std`: the standard deviation in energy or force components

```Python
import pinn
from ase.io import read
from tips.bias import EnsembleBiasedCalculator
models = [f'examples/ensemble/pinet-{i+1}/model' for i in range(5)]
calcs = [pinn.get_calc(model) for model in models]
calc = EnsembleBiasedCalculator(calcs, bias='heaviside', kb=100)
atoms = read("examples/ensemble/water.xyz")
atoms.get_property('forces_bias')
```

In the above example, we used a biasing potential to prevent the potential from
sampling areas with a high disagreement, the available biasing potentials are
listed below.

### Heaviside biasing

The `'heaviside'` biasing applies a potential defined in by Schran et al.
[@2020_SchranBrezinaMarsalek]

\begin{align}
E^{(b)} &=
\theta(\sigma_E -\sigma_0) \frac{1}{2} k^{b} (\sigma_E -\sigma_0)^2 \\
-\nabla_\alpha E^{(b)} &=
\theta(\sigma_E -\sigma_0) k^{b} \frac{\sigma_E -\sigma_0}{\sigma_E} \cdot
\frac{1}{n}\sum_{i=1}^{n} - \Delta E_i \nabla_\alpha \Delta E_i
\end{align}

where $E^{(b)}$ is the biasing potential,
$\nabla_\alpha E^{(b)}$ is the corresponding biasing force;
$\sigma_E$ and $\Delta E_i$ are the standard deviation and individual disagreement
in total energy predictions, defined as:

\begin{align}
\sigma_E &= \left[\frac{1}{n} \sum_{i=1}^n (\Delta E_i)^2 \right]^{1/2} \\
\Delta E_i &= E - E_i
\end{align}

??? "Source code"

    ```groovy
    --8<-- "python/tips/bias/ase.py"
    ```

\bibliography
