# -*- coding: utf-8 -*-

import numpy as np
from ase.calculators.calculator import Calculator

class EnsembleBiasedCalculator(Calculator):
    """ This is an implementation of ensemble calculator which
    combines multiple calcualtors and optionally baising the potential
    according the disagreement.
    """

    def __init__(self, calcs, bias=None, kb=0, sigma0=0):
        """
        Args:
            calcs: List of `ase.calcualtors`
            bias: None or 'heaviside', see the documentatin for details.
            kb: harmonic form factor in heaviside-style biasing
            sigma0: tolerance in heaviside-style biasing
        """
        super().__init__()
        self.n = len(calcs)
        self.calcs = calcs
        self.orig_properties = ['energy', 'forces', 'stress']
        self.implemented_properties = [
            'energy', 'energy_std', 'energy_avg', 'energy_bias',
            'forces', 'forces_std', 'forces_avg', 'forces_bias',
            'stress', 'stress_std', 'stress_avg', 'stress_bias'
        ]
        self.bias = bias
        self.kb = kb
        self.sigma0 = sigma0

    def get_orig_properties(self, properties, atoms):
        results = {}
        for prop in self.orig_properties:
            contrib = [calc.get_property(prop, atoms) for calc in self.calcs]
            results[f'{prop}_avg'] = np.mean(contrib, axis=0)
            results[f'{prop}_std'] = np.std(contrib, axis=0)
            results[f'{prop}_all'] = contrib
        return results

    def get_properties(self, properties, atoms):
        results = self.get_orig_properties(properties, atoms)
        if self.bias is None:
            for prop in self.orig_properties:
                results[prop] = results[f'{prop}_avg']
                results[f'{prop}_bias'] = np.zeros_like(results[f'{prop}_avg'])
        elif self.bias == 'heaviside':
            n = self.n
            kb = self.kb
            sigma0 = self.sigma0
            sigmaE = results['energy_std']
            theta = np.heaviside(sigmaE-sigma0, 0.5)
            DeltaE = results['energy_avg'] - results['energy_all']
            results['energy_bias'] = 0.5 * theta * kb * (sigmaE - sigma0)**2

            prefac = theta * kb * (sigmaE - sigma0) / n
            for prop in ['forces', 'stress']: # direvarive properties
                DeltaP = [results[f'{prop}_avg'] - P for P in results[f'{prop}_all']]
                results[f'{prop}_bias'] = prefac * sum(DE*DP for DE, DP in zip(DeltaE, DeltaP))
            for prop in self.orig_properties:
                results[prop] = results[f'{prop}_avg'] + results[f'{prop}_bias']
        else:
            raise ValueError(f"Unknown biasing method: {self.bias}")

        return results

    def calculate(self, atoms, properties, system_changes):
        """Calculates all the specific property for each calculator and
        returns with the summed value.

        """
        self.atoms = atoms.copy()  # for caching of results
        self.results = self.get_properties(properties, atoms)
