################################################################################
# Raffaele Cheula, cheula.raffaele@gmail.com
################################################################################

import numpy as np
from ase import units
from ase.thermochemistry import ThermoChem
from shape.ase_utils import get_moments_of_inertia_xyz

################################################################################
# HINDERED THERMO CLASS
################################################################################

class HinderedThermo(ThermoChem):

    def __init__(self, vib_energies, site_area, trans_barrier_energy = None,
                 rot_barrier_energy = None, rotationalminima = None,
                 potentialenergy = 0., atoms = None, mass = None,
                 inertia = None, rotation_center = None, symmetrynumber = 1,
                 with_std_conc_entropy = False, indices = None):

        if rot_barrier_energy is None and trans_barrier_energy is not None:
            if indices:
                if len(indices) > 2:
                    raise RuntimeError('you cannot specify more than '
                                       '2 translations')
                else: 
                    self.vib_energies = sorted(vib_energies,
                        key = lambda x: x.real) 
                    for i in sorted(indices, reverse = True):
                        del self.vib_energies[i]
            else:
                self.vib_energies = sorted(vib_energies,
                    key = lambda x: x.real, reverse = True)[:-2]
            self.trans_barrier_energy = trans_barrier_energy * units._e
            self.rot_barrier_energy = None            

        elif rot_barrier_energy is not None and trans_barrier_energy is None:
            if indices:
                if len(indices) > 1:
                    raise RuntimeError('you cannot specify use more than '
                                       '1 rotation')
                else: 
                    self.vib_energies = sorted(vib_energies,
                        key = lambda x: x.real)  
                    for i in sorted(indices, reverse = True):
                        del self.vib_energies[i]
            else:  
                self.vib_energies = sorted(vib_energies,
                    key = lambda x: x.real, reverse = True)[:-1]
            self.trans_barrier_energy = None  
            self.rot_barrier_energy = rot_barrier_energy * units._e
            self.rotationalminima = rotationalminima
            self.symmetry = symmetrynumber

        elif rot_barrier_energy is not None and \
            trans_barrier_energy is not None:
            if indices:
                if len(indices) > 3:
                    raise RuntimeError('you cannot specify more than 2 '
                                       'translations and 1 rotation')
                else: 
                    self.vib_energies = sorted(vib_energies,
                        key = lambda x: x.real)
                    for i in sorted(indices, reverse = True):
                        del self.vib_energies[i]
            else:    
                self.vib_energies = sorted(vib_energies,
                    key = lambda x: x.real, reverse = True)[:-3]
            self.trans_barrier_energy = trans_barrier_energy * units._e
            self.rot_barrier_energy = rot_barrier_energy * units._e
            self.rotationalminima = rotationalminima
            self.symmetry = symmetrynumber
        else:
            raise RuntimeError('Either rot_barrier_energy must be specified or '
                               'trans_barrier_energy must be specified.')

        
        self.area = site_area / units.m**2
        self.potentialenergy = potentialenergy
        self.atoms = atoms
        self.with_std_conc_entropy = with_std_conc_entropy

        if (mass or atoms) and (inertia or atoms):
            if mass:
                self.mass = mass * units._amu
            elif atoms:
                self.mass = np.sum(atoms.get_masses()) * units._amu
            if inertia:
                self.inertia = inertia * units._amu / units.m**2
            elif atoms:
                I = get_moments_of_inertia_xyz(atoms, center = rotation_center)
                self.inertia = (I[2] * units._amu / units.m**2)

        else:
            raise RuntimeError('Either mass and inertia of the '
                               'adsorbate must be specified or '
                               'atoms must be specified.')

        # Make sure no imaginary frequencies remain.
        if sum(np.iscomplex(self.vib_energies)):
            raise ValueError('Imaginary frequencies are present.')
        else:
            self.vib_energies = np.real(self.vib_energies)  # clear +0.j

        # Calculate hindered translational and rotational frequencies
        if self.trans_barrier_energy is not None:
            self.freq_t = np.sqrt(self.trans_barrier_energy / (2 * self.mass *
                                                           self.area))
        
        if self.rot_barrier_energy not in (0., None):
            self.freq_r = 1. / (2 * np.pi) * np.sqrt(self.rotationalminima**2 *
                          self.rot_barrier_energy / (2 * self.inertia))

    def get_internal_energy(self, temperature, verbose = True):
        """Returns the internal energy (including the zero point energy),
        in eV, in the hindered translator and hindered rotor model at a
        specified temperature (K)."""

        from scipy.special import iv

        self.verbose = verbose
        write = self._vprint
        fmt = '%-15s%13.3f eV'
        write('Internal energy components at T = %.2f K:' % temperature)
        write('=' * 31)

        U = 0.

        write(fmt % ('E_pot', self.potentialenergy))
        U += self.potentialenergy

        # Translational Energy
        if self.trans_barrier_energy is not None:
            T_t = units._k * temperature / (units._hplanck * self.freq_t)
            R_t = self.trans_barrier_energy / (units._hplanck * self.freq_t)
            dU_t = 2 * (-1. / 2 - 1. / T_t / (2 + 16 * R_t) + R_t / 2 / T_t -
                        R_t / 2 / T_t *
                        iv(1, R_t / 2 / T_t) / iv(0, R_t / 2 / T_t) +
                        1. / T_t / (np.exp(1. / T_t) - 1))
            dU_t *= units.kB * temperature
            write(fmt % ('E_trans', dU_t))
            U += dU_t

        # Rotational Energy
        if self.rot_barrier_energy == 0.:
            dU_r = units.kB * temperature
            write(fmt % ('E_rot', dU_r))
            U += dU_r
        elif self.rot_barrier_energy is not None:
            T_r = units._k * temperature / (units._hplanck * self.freq_r)
            R_r = self.rot_barrier_energy / (units._hplanck * self.freq_r)
            dU_r = (-1. / 2 - 1. / T_r / (2 + 16 * R_r) + R_r / 2 / T_r -
                    R_r / 2 / T_r *
                    iv(1, R_r / 2 / T_r) / iv(0, R_r / 2 / T_r) +
                    1. / T_r / (np.exp(1. / T_r) - 1))
            dU_r *= units.kB * temperature
            write(fmt % ('E_rot', dU_r))
            U += dU_r

        # Vibrational Energy
        dU_v = self._vibrational_energy_contribution(temperature)
        write(fmt % ('E_vib', dU_v))
        U += dU_v

        # Zero Point Energy
        dU_zpe = self.get_zero_point_energy()
        write(fmt % ('E_ZPE', dU_zpe))
        U += dU_zpe

        write('-' * 31)
        write(fmt % ('U', U))
        write('=' * 31)

        return U


    def get_zero_point_energy(self, verbose=True):
        """Returns the zero point energy, in eV, in the hindered
        translator and hindered rotor model"""
        if self.trans_barrier_energy is not None:
            zpe_t = 2 * (1. / 2 * self.freq_t * units._hplanck / units._e)
        else:
            zpe_t = 0.    
        if self.rot_barrier_energy not in (0., None):
            zpe_r = 1. / 2 * self.freq_r * units._hplanck / units._e
        else:
            zpe_r = 0.
        zpe_v = self.get_ZPE_correction()
        zpe = zpe_t + zpe_r + zpe_v

        return zpe


    def get_entropy(self, temperature, verbose = True):
        """Returns the entropy, in eV/K, in the hindered translator
        and hindered rotor model at a specified temperature (K)."""

        from scipy.special import iv

        self.verbose = verbose
        write = self._vprint
        fmt = '%-15s%13.7f eV/K%13.3f eV'
        write('Entropy components at T = %.2f K:' % temperature)
        write('=' * 49)
        write('%15s%13s     %13s' % ('', 'S', 'T*S'))

        S = 0.

        # Translational Entropy
        if self.trans_barrier_energy is not None:
            T_t = units._k * temperature / (units._hplanck * self.freq_t)
            R_t = self.trans_barrier_energy / (units._hplanck * self.freq_t)
            S_t = 2 * (-1. / 2 + 1. / 2 * np.log(np.pi * R_t / T_t) -
                       R_t / 2 / T_t *
                       iv(1, R_t / 2 / T_t) / iv(0, R_t / 2 / T_t) +
                       np.log(iv(0, R_t / 2 / T_t)) +
                       1. / T_t / (np.exp(1. / T_t) - 1) -
                       np.log(1 - np.exp(-1. / T_t)))
            S_t *= units.kB
            write(fmt % ('S_trans', S_t, S_t * temperature))
            S += S_t

        # Rotational Entropy
        if self.rot_barrier_energy == 0.:
            S_r = (1 + np.log((temperature * 8 * np.pi**2 * self.inertia * 
                    units._k) / (self.symmetry * units._hplanck**2)))
            S_r *= units.kB
            write(fmt % ('S_rot', S_r, S_r * temperature))
            S += S_r
        elif self.rot_barrier_energy is not None:
            T_r = units._k * temperature / (units._hplanck * self.freq_r)
            R_r = self.rot_barrier_energy / (units._hplanck * self.freq_r)
            S_r = (-1. / 2 + 1. / 2 * np.log(np.pi * R_r / T_r) -
                np.log(self.symmetry) -
                R_r / 2 / T_r * iv(1, R_r / 2 / T_r) / iv(0, R_r / 2 / T_r) +
                np.log(iv(0, R_r / 2 / T_r)) +
                1. / T_r / (np.exp(1. / T_r) - 1) -
                np.log(1 - np.exp(-1. / T_r)))
            S_r *= units.kB
            write(fmt % ('S_rot', S_r, S_r * temperature))
            S += S_r

        # Vibrational Entropy
        S_v = self._vibrational_entropy_contribution(temperature)
        write(fmt % ('S_vib', S_v, S_v * temperature))
        S += S_v

        # Concentration Related Entropy
        if self.with_std_conc_entropy is True:
            N_over_A = np.exp(1. / 3) * (10.0**5 /
                                        (units._k * temperature))**(2. / 3)
            S_c = 1 - np.log(N_over_A) - np.log(self.area)
            S_c *= units.kB
            write(fmt % ('S_con', S_c, S_c * temperature))
            S += S_c

        write('-' * 49)
        write(fmt % ('S', S, S * temperature))
        write('=' * 49)

        return S


    def get_helmholtz_energy(self, temperature, verbose = True):
        """Returns the Helmholtz free energy, in eV, in the hindered
        translator and hindered rotor model at a specified temperature
        (K)."""

        self.verbose = True
        write = self._vprint

        U = self.get_internal_energy(temperature, verbose = verbose)
        write('')
        S = self.get_entropy(temperature, verbose = verbose)
        F = U - temperature * S

        write('')
        write('Free energy components at T = %.2f K:' % temperature)
        write('=' * 23)
        fmt = '%5s%15.3f eV'
        write(fmt % ('U', U))
        write(fmt % ('-T*S', -temperature * S))
        write('-' * 23)
        write(fmt % ('F', F))
        write('=' * 23)

        return F

################################################################################
# END
################################################################################
