"""
Theoretical simulator for telescopes

Assuming NGTS-like performance, what precision can be reached
with what type of telescope array? E.g. what configurations
can reach 10 ppm and what are the costs?
"""
import numpy as np
import matplotlib.pyplot as plt
from publications.plots.defaults import (
    general,
    one_by_three
    )

# TODO: add costings in somehow

class Common:
    """
    Common site, sky and scaling values
    """
    photons_vega = 1000 # s^-1 cm^-2 A^-1
    sky_brightness = 21.3 # mag arcsec^-2
    h_obs = 2500 # m
    h_scale = 8000 # m
    airmass = 1.20 # ---
    scint_coeff_pao = 1.56 # ---
    exptime_b = 3600 # s
    n_telescopes = 12 # ---
    mag_range = np.arange(5, 22.25, 0.25)

class NGTS:
    """
    Data class for NGTS
    """
    telescope_id = "NGTS"
    diameter = 200 # mm
    focal_length = 280 # mm
    bandwidth = 3700 # A
    pixel_size = 13.5 # microns
    fwhm = 3 # pixels
    total_efficiency = 65 # %
    exptime = 10 # s
    readout_time = 2.5 # s
    dark_current = 1 # e pixel^-1 s^-1
    read_noise = 7 # e pixel^-1

class NGTS1m:
    """
    Data class for NGTS-like 1m telescopes
    """
    telescope_id = "NGTS1m"
    diameter = 1000 # mm
    focal_length = 2800 # mm
    bandwidth = 3700 # A
    pixel_size = 13.5 # microns
    fwhm = 3 # pixels
    total_efficiency = 65 # %
    exptime = 10 # s
    readout_time = 2.5 # s
    dark_current = 1 # e pixel^-1 s^-1
    read_noise = 7 # e pixel^-1

class NGTS2m:
    """
    Data class for NGTS-like 2m telescopes
    """
    telescope_id = "NGTS2m"
    diameter = 2000 # mm
    focal_length = 5600 # mm
    bandwidth = 3700 # A
    pixel_size = 13.5 # microns
    fwhm = 3 # pixels
    total_efficiency = 65 # %
    exptime = 10 # s
    readout_time = 2.5 # s
    dark_current = 1 # e pixel^-1 s^-1
    read_noise = 7 # e pixel^-1

if __name__ == "__main__":
    # create telescope or set of telescopes
    telescopes = [NGTS, NGTS1m, NGTS2m]

    # use some plot defaults
    general()
    one_by_three()
    # make some nice plots
    fig, ax = plt.subplots(1, len(telescopes))

    # loop over each telescope
    for t, telescope in enumerate(telescopes):
        # calculate the derived parameters
        print(f"{telescope.telescope_id} preamble:")
        signal_vega = Common.photons_vega * telescope.bandwidth * (telescope.total_efficiency/100) * np.pi * ((telescope.diameter/10/2)**2)
        print(f"\tVega: {signal_vega:.2f}")
        platescale = 206.264806 * telescope.pixel_size / telescope.focal_length
        print(f"\tPlatescale: {platescale:.3f}")
        focal_ratio = telescope.focal_length / telescope.diameter
        print(f"\tf: {focal_ratio:.2f}")
        n_pixels_in_aper = np.pi * telescope.fwhm**2
        print(f"\tN pixels in aper: {n_pixels_in_aper:.2f}")

        # calculate the noise model components

        # place holders for the model components
        target_signal = np.empty_like(Common.mag_range)
        sky_signal = np.empty_like(Common.mag_range)
        read_noise_signal = np.empty_like(Common.mag_range)
        dark_signal = np.empty_like(Common.mag_range)
        scint_signal = np.empty_like(Common.mag_range)
        total_noise = np.empty_like(Common.mag_range)
        # scintillation fractional error
        scint_frac = np.sqrt(0.00001 * (Common.scint_coeff_pao**2) * ((telescope.diameter/1000.)**(-4./3.)) * (1./telescope.exptime) * (Common.airmass**3) * np.exp(-2.*Common.h_obs/Common.h_scale))
        # loop to build up noise model
        for i, mag in enumerate(Common.mag_range):
            target_signal[i] = (signal_vega / 2.5**mag) * telescope.exptime
            sky_signal[i] = (signal_vega * 10**(-0.4 * Common.sky_brightness)) * (platescale**2) * telescope.exptime * n_pixels_in_aper
            read_noise_signal[i] = (telescope.read_noise**2) * n_pixels_in_aper
            dark_signal[i] = telescope.dark_current * telescope.exptime * n_pixels_in_aper
            scint_signal[i] = (scint_frac * target_signal[i])**2
        # get the combined noise model
        total_noise = np.sqrt(target_signal + sky_signal + read_noise_signal + dark_signal + scint_signal)
        # get fractional noise model in ppm per exposure
        total_noise_ppm = (total_noise / target_signal) * 1E6
        # get fractional noise model in ppm per binned exposure
        total_noise_ppm_binned = total_noise_ppm / (np.sqrt(Common.exptime_b / (telescope.exptime + telescope.readout_time)))
        # get fractional noise model in ppm per binned exposure per ntel
        total_noise_ppm_binned_ntel = total_noise_ppm_binned / np.sqrt(Common.n_telescopes)

        # get noise components as fractional errors for plotting
        target_noise_ppm = (np.sqrt(target_signal)/target_signal) * 1E6
        sky_noise_ppm = (np.sqrt(sky_signal)/target_signal) * 1E6
        read_noise_ppm = (np.sqrt(read_noise_signal)/target_signal) * 1E6
        dark_noise_ppm = (np.sqrt(dark_signal)/target_signal) * 1E6
        scint_noise_ppm = (np.sqrt(scint_signal)/target_signal) * 1E6

        ax[t].semilogy(Common.mag_range, total_noise_ppm, 'k-', label=f'Total noise', alpha=0.3)
        ax[t].semilogy(Common.mag_range, target_noise_ppm, 'g--', label='Photon noise', lw=1, alpha=0.3)
        ax[t].semilogy(Common.mag_range, sky_noise_ppm, 'b--', label='Sky noise', lw=1, alpha=0.3)
        ax[t].semilogy(Common.mag_range, read_noise_ppm, '--', label='Read noise', color='brown', lw=1, alpha=0.3)
        ax[t].semilogy(Common.mag_range, dark_noise_ppm, '--', label='Dark noise', color='purple', lw=1, alpha=0.3)
        ax[t].semilogy(Common.mag_range, scint_noise_ppm, '--', label='Scint noise', color='cyan', lw=1, alpha=0.3)
        ax[t].semilogy(Common.mag_range, total_noise_ppm_binned, 'k-', label=f'Total noise ({Common.exptime_b}s 1t)', alpha=0.5)
        ax[t].semilogy(Common.mag_range, total_noise_ppm_binned_ntel, 'k-', label=f'Total noise ({Common.exptime_b}s 12t)')
        ax[t].axhline(50, ls='--', color='red', label='PLATO (1h)', lw=1)
        ax[t].set_ylabel("Noise (ppm)")
        ax[t].set_xlabel("I-band (mag)")
        ax[t].set_title(f"{telescope.telescope_id}")
        ax[t].set_ylim(1, 100000)
        ax[t].legend(fontsize='x-small')
    fig.tight_layout()
    fig.savefig(f'noise_model_telescope_simulator.png')
