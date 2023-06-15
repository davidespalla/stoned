import numpy as np


def compute_ratemap1d(s, x, t, bins):
    """
    Compute the rate map for 1-dimensional spatial data.

    Parameters
    ----------
    s : list
        Spike times.

    x : array-like
        Array of occupied positions.

    t : array-like
        Times associated with the position values, should be in the same units of s.

    bins : array-like
        Spatial bins for dividing the 1-dimensional space, should be in the same units of x.

    Returns
    -------
    ratemap : ndarray
        The computed rate map.

    Notes
    -----
    The rate map is a representation of the average firing rate as a function of position.
    The rate map is obtained by dividing the spike count histogram by the occupancy,
    ignoring any division errors (resulting in NaN values).

    Examples
    --------
    >>> s = [1,12.3,50.1,60.3] # Spike times
    >>> x = np.linspace(0,1) # Position values
    >>> t = np.linspace(0,100)  # Times
    >>> bins = np.linspace(0, 1, 10)  # Spatial bins
    >>> ratemap = compute_ratemap1d(s, x, t, bins)

    """

    dt = t[1]-t[0]  # compute dt from time axis
    spike_positions = np.interp(s, x, t)  # find spike positions
    # compute histogram of spike counts
    spike_hist = np.histogram(spike_positions, bins)[0]
    occupancy = np.histogram(x, bins)[0]*dt  # compute occupancy of positions

    # divide by occupancy, do not display division error (will be nan values)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratemap = spike_hist / occupancy

    return ratemap


def spatial_info_perspike(rate_map, occupancy_prob, epsilon=pow(10, -15)):
    """
    Compute the Skaggs information (bits per spike) based on the rate map and occupancy probability.

    Parameters
    ----------
    rate_map : ndarray
        The rate map representing the firing rate at each spatial location.

    occupancy_prob : ndarray
        The occupancy probability corresponding to each spatial location (must sum to 1 and be of the
        same shape as rate_map).

    epsilon : float, optional
        A small value added to the rate map to avoid division by zero. Defaults to 1e-15.

    Returns
    -------
    info_per_spike : float
        The computed Skaggs information in bits per spike.

    Notes
    -----
    The Skaggs information quantifies the spatial information content per spike.
    The rate map and occupancy probability should have the same shape.
    If the average rate is zero, the result will be NaN.

    Examples
    --------
    >>> rate_map = np.array([0.5, 1.2, 0.8, 0.3])  # Example rate map
    # Example occupancy probability
    >>> occupancy_prob = np.array([0.1, 0.3, 0.4, 0.2])
    >>> info_per_spike = spatial_info_perspike(rate_map, occupancy_prob, epsilon=1e-15)

    """

    rate_map = rate_map.flatten()
    occupancy_prob = occupancy_prob.flatten()

    assert(np.sum(occupancy_prob) == 1, "Sum of occupancy_prob must be = 1")

    avg_rate = np.mean(rate_map)
    if avg_rate > 0:
        return sum(rate_map*np.log2((rate_map+epsilon)/avg_rate)*occupancy_prob)/avg_rate
    else:
        return np.nan


def spatial_info_persec(rate_map, occupancy_prob, epsilon=pow(10, -15)):
    """
    Compute the Skaggs information (bits per second) based on the rate map and occupancy probability.

    Parameters
    ----------
    rate_map : ndarray
        The rate map representing the firing rate at each spatial location.

    occupancy_prob : ndarray
        The occupancy probability corresponding to each spatial location (must sum to 1 and be of the
        same shape as rate_map).

    epsilon : float, optional
        A small value added to the rate map to avoid division by zero. Defaults to 1e-15.

    Returns
    -------
    info_per_second : float
        The computed Skaggs information in bits per second.

    Notes
    -----
    The Skaggs information quantifies the spatial information content per spike.
    The rate map and occupancy probability should have the same shape.
    If the average rate is zero, the result will be NaN.

    Examples
    --------
    >>> rate_map = np.array([0.5, 1.2, 0.8, 0.3])  # Example rate map
    # Example occupancy probability
    >>> occupancy_prob = np.array([0.1, 0.3, 0.4, 0.2])
    >>> info_per_secodn = spatial_info_persec(rate_map, occupancy_prob, epsilon=1e-15)

    """

    rate_map = rate_map.flatten()
    occupancy_prob = occupancy_prob.flatten()

    assert(np.sum(occupancy_prob) == 1, "Sum of occupancy_prob must be = 1")

    avg_rate = np.mean(rate_map)
    if avg_rate > 0:
        return sum(rate_map*np.log2((rate_map+epsilon)/avg_rate)*occupancy_prob)
    else:
        return np.nan


############### TO DO ###########################

def compute_ratemap2d(spike_times, positions, pos_times, nbins=[20, 20], range=None):
    spike_positions = spike_positions_2d(spike_times, positions, pos_times)
    if not range:
        range = [[min(positions[:, 0]), max(positions[:, 0])], [
            min(positions[:, 1]), max(positions[:, 1])]]
    occupancy = np.histogram2d(
        positions[:, 0], positions[:, 1], bins=nbins, range=range)[0]
    spikemap = np.histogram2d(
        spike_positions[0], spike_positions[1], bins=nbins, range=range)[0]

    with np.errstate(divide='ignore', invalid='ignore'):
        rate_map = spikemap / occupancy

    rate_map[np.logical_or(np.isnan(rate_map), np.isinf(rate_map))] = 0

    i

    return ratemap


def compute_occupancy(positions, nbins=[20, 20], range=None, sigma=None):
    if not range:
        range = [[min(positions[:, 0]), max(positions[:, 0])], [
            min(positions[:, 1]), max(positions[:, 1])]]

    occupancy = np.histogram2d(
        positions[:, 0], positions[:, 1], bins=nbins, range=range)[0]

    if sigma:
        occupancy = gaussian_filter(occupancy, sigma)

    return occupancy


def compute_spikemap(spike_times, positions, pos_times, nbins=[20, 20], range=None, sigma=None):
    spike_positions = spike_positions_2d(spike_times, positions, pos_times)
    if not range:
        range = [[min(positions[:, 0]), max(positions[:, 0])], [
            min(positions[:, 1]), max(positions[:, 1])]]
    spikemap = np.histogram2d(
        spike_positions[0], spike_positions[1], bins=nbins, range=range)[0]
    if sigma:
        spikemap = gaussian_filter(spikemap, sigma)

    return spikemap
