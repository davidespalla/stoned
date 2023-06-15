"""
Functions for local field potential (LFP) time-frequency analysis
"""

import numpy as np
from scipy.signal import butter
from scipy.signal import sosfilt
from pywt import scale2frequency
from pywt import cwt


def bandpass_filter(signal, low_f, high_f, sampling_rate=1.0, filter_order=5):
    """
    Apply a bandpass filter to the input signal using a Butterworth filter.

    Parameters
    ----------
    signal : array-like
        The input signal to filter.

    low_f : float
        The lower bound of the frequency band.

    high_f : float
        The upper bound of the frequency band.

    sampling_rate : float, optional
        The sampling rate of the signal. Defaults to 1.0 if not specified.

    filter_order : int, optional
        The order of the Butterworth filter. Defaults to 5 if not specified.

    Returns
    -------
    filtered_signal : ndarray
        The filtered signal.

    Notes
    -----
    The Butterworth filter is implemented using the `sosfilt` function from `scipy.signal`.

    Examples
    --------
    >>> signal = np.random.randn(1000)  # Random input signal
    >>> low_f = 1.0  # Lower bound of the frequency band
    >>> high_f = 10.0  # Upper bound of the frequency band
    >>> filtered_signal = bandpass_filter(signal, low_f, high_f, sampling_rate=100.0)

    """
    filter = butter(filter_order, [low_f, high_f],
                    btype='band', output='sos', fs=sampling_rate)
    filtered_signal = sosfilt(filter, signal)
    return filtered_signal


def morlet_transform(signal, low_f, high_f, n_freqs=20, sampling_rate=1.0):
    """
    Apply the Morlet transform to the input signal within a specified frequency band.

    Parameters
    ----------
    signal : array-like
        The input signal to transform.

    low_f : float
        The lower bound of the frequency band.

    high_f : float
        The upper bound of the frequency band.

    n_freqs : int, optional
        The number of frequencies to use for the transform. Defaults to 20.

    sampling_rate : float, optional
        The sampling rate of the signal. Defaults to 1.0 if not specified.

    Returns
    -------
    C : ndarray
        The complex Morlet transform coefficients.

    freq : ndarray
        The frequencies associated with the transform coefficients.

    Notes
    -----
    The Morlet transform is performed using the continuous wavelet transform (CWT)
    with the 'cmor1-0.5' wavelet.

    Examples
    --------
    >>> signal = np.random.randn(1000)  # Random input signal
    >>> low_f = 1.0  # Lower bound of the frequency band
    >>> high_f = 10.0  # Upper bound of the frequency band
    >>> C, freq = morlet_transform(signal, low_f, high_f, n_freqs=20, sampling_rate=100.0)

    """
    frequencies = np.linspace(low_f, high_f, n_freqs)/sampling_rate
    scales = scale2frequency('cmor1-0.5', frequencies)
    C, freq = cwt(signal, wavelet='cmor1-0.5', scales=scales,
                  sampling_period=1.0/sampling_rate)
    return C, freq


def compute_power(signal, low_f, high_f, sampling_rate=1.0, n_freqs=20):
    """
    Compute the instantaneous power of the signal within a specified frequency band.

    Parameters
    ----------
    signal : array-like
        The input signal.

    low_f : float
        The lower bound of the frequency band.

    high_f : float
        The upper bound of the frequency band.

    sampling_rate : float, optional
        The sampling rate of the signal. Defaults to 1.0 if not specified.

    n_freqs : int, optional
        The number of frequencies used for the Morlet transform used to compute the power.
        Defaults to 20.

    Returns
    -------
    power : ndarray
        The computed power timecourse.

    Notes
    -----
    The power is computed by performing a Morlet transform on the input signal
    within the specified frequency band. The average absolute value of the transform
    coefficients is taken along the frequency axis to obtain the power at each timepoint.

    Examples
    --------
    >>> signal = np.random.randn(1000)  # Random input signal
    >>> low_f = 1.0  # Lower bound of the frequency band
    >>> high_f = 10.0  # Upper bound of the frequency band
    >>> power_spectrum = compute_power(signal, low_f, high_f, sampling_rate=100.0, n_freqs=20)
    """

    C, freq = morlet_transform(
        signal, low_f, high_f, sampling_rate=sampling_rate, n_freqs=n_freqs)
    power = np.mean(abs(C), axis=0)
    return power
