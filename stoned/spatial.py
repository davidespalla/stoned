import numpy as np


def compute_ratemap1d(spike_list, x, t, bins):
    '''
    s: 2d list
    x: array-like
    t: array-like 
    bins: array-like 

    return: ratemaps

    '''
    # compute dt from time axis
    dt = t[1]-t[0]
    # find spike positions
    spike_positions = [np.interp(s, x, t) for s in spike_list]
    # compute histogram of spike counts and put it into matrix
    spikes_hist = [np.histogram(s, bins)[0] for s in spike_positions]
    spikes_hist = np.vstack(spikes_hist).astype(np.float64)

    # compute occupancy
    occupancy = np.histogram(x, bins)[0]*dt

    # divide by occupancy, do not display division error (will be nan values)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratemaps = spikes_hist / occupancy

    # return ratemaps as matrix (neurons x bins)
    return ratemaps


def spatial_info_perspike(rate_map, occupancy_prob, epsilon=pow(10, -15)):
    '''
    takes ratemap and occupancy probability, returns skaggs information (bit/spike)
    '''
    rate_map = rate_map.flatten()
    occupancy_prob = occupancy_prob.flatten()
    avg_rate = np.mean(rate_map*occupancy_prob)
    if avg_rate > 0:
        return sum(rate_map*np.log2((rate_map+epsilon)/avg_rate)*occupancy_prob)/avg_rate
    else:
        return np.nan
    
def spatial_info_persec(rate_map, occupancy_prob, epsilon=pow(10, -15)):
    '''
    takes ratemap and occupancy probability, returns skaggs information (bit/seconds)
    '''
    rate_map = rate_map.flatten()
    occupancy_prob = occupancy_prob.flatten()
    avg_rate = np.mean(rate_map*occupancy_prob)
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
