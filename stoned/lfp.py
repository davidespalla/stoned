import numpy as np
import random
from scipy.ndimage import gaussian_filter1d
from collections import Counter

def prepare_data(trajectory,spikes,times,resampling_factor):
    #builds spike count matrix
    start_time = times[0]
    end_time = times[-1]
    n_bins = int(np.floor(trajectory.shape[1]/resampling_factor))
    X = np.empty((len(spikes),n_bins))
    time_bins = np.linspace(start_time,end_time,n_bins+1)
    for cell in range(len(spikes)):
           X[cell][:] = np.histogram(spikes[cell],bins=time_bins)[0]
            
    X = X.T
        
            
    # resample trajectory at same timescale as spike count matrix
    resampled_trajectory = np.empty((trajectory.shape[0],n_bins))
    for i in range(n_bins):
        resampled_trajectory[0][i] = np.mean(trajectory[0][resampling_factor*i:resampling_factor*(i+1)])
        resampled_trajectory[1][i] = np.mean(trajectory[1][resampling_factor*i:resampling_factor*(i+1)])
        
               
    return resampled_trajectory,X

def resample_ethogram(ethogram,times,resampling_factor):
    start_time = times[0]
    end_time = times[-1]
    n_bins=int(np.floor(len(ethogram)/resampling_factor))
            
    # use majority rule to resample ethogram
    resampled_ethogram=np.empty(n_bins)
    for i in range(n_bins):
        counter = Counter(ethogram[resampling_factor*i:resampling_factor*(i+1)])
        resampled_ethogram[i]=counter.most_common(1)[0][0]
      
    return resampled_ethogram


def bin_spikes(spikes,times,resampling_factor):
    #builds spike count matrix
    start_time = times[0]
    end_time = times[-1]
    n_bins = int(np.floor(times.shape[0]/resampling_factor))
    X = np.empty((len(spikes),n_bins))
    time_bins = np.linspace(start_time,end_time,n_bins+1)
    for cell in range(len(spikes)):
           X[cell][:] = np.histogram(spikes[cell],bins=time_bins)[0]
            
    X = X.T    
               
    return X

def calculate_speed(x,y,dt,smooth = False,sigma =1 ):
    v=np.zeros(len(x))
    for i in range(1,len(v)):
        v[i]=np.sqrt(pow(x[i]-x[i-1],2)+pow(y[i]-y[i-1],2))/dt
    if smooth:
        v = gaussian_filter1d(v,sigma)

    return v

def exlude_high_firing(spikes,tot_time,rate_th):
    filtered_spikes=[]
    for cell in range(len(spikes)):
        if len(spikes[cell])/tot_time < rate_th:
            filtered_spikes.append(spikes[cell])
        else:
            filtered_spikes.append(np.nan)
    return filtered_spikes

def speed_mask(x,y,frame_dt,pixel_to_cm,velocity_th=3):
    conversion_factor = pixel_to_cm/frame_dt
    v=np.zeros(len(x))
    mask=np.zeros(len(x))
    for i in range(1,len(v)):
        v[i]=np.sqrt(pow(x[i]-x[i-1],2)+pow(y[i]-y[i-1],2))*conversion_factor
        if v[i]>=velocity_th:
            mask[i]=1
    return mask

def box_mask(x,y,box_range):
    xmin = box_range[0],
    xmax = box_range[1]
    ymin = box_range[2]
    ymax = box_range[3]
    #1120,ymin=200,ymax=870
    mask=np.zeros(len(x))
    for i in range(len(mask)):
        if x[i]>xmin and x[i]<xmax and y[i]>ymin and y[i]<ymax:
            mask[i]=1
    return mask

def align_trajectory(trajectory,recording_mask):
    out_x = trajectory[0][:-1][recording_mask]
    out_y = trajectory[1][:-1][recording_mask]
    return np.asarray([out_x,out_y])

def align_facing(facing,recording_mask):
    return facing[:-1][recording_mask]

def align_ethogram(ethogram,recording_mask):
    return ethogram[:-1][recording_mask]

def enhance_box_range(box_range,pixels):
    out = box_range
    out[0] = out[0] -pixels
    out[1] = out[1] +pixels
    out[2] = out[2] -pixels
    out[3] = out[3] +pixels
    return out

def shuffle_activity(X):
    Xs = np.copy(X)
    Xs = Xs.T
    for i in range(len(Xs)):
        Xs[i][:] = random.sample(list(Xs[i]),len(Xs[i]))
    return Xs.T


def skaggs_info_persec(rate_map,occupancy_prob,epsilon=pow(10,-15)):
    rate_map = rate_map.flatten()
    occupancy_prob = occupancy_prob.flatten()
    avg_rate = np.mean(rate_map)
    if avg_rate >0:
        return sum(rate_map*np.log2((rate_map+epsilon)/avg_rate)*occupancy_prob)
    else:
        return np.nan

def skaggs_info_perspike(rate_map,occupancy_prob,epsilon=pow(10,-15)):
    rate_map = rate_map.flatten()
    occupancy_prob = occupancy_prob.flatten()
    avg_rate = np.mean(rate_map)
    if avg_rate >0:
        return sum(rate_map*np.log2((rate_map+epsilon)/avg_rate)*occupancy_prob)/avg_rate
    else:
        return np.nan
    
def cosine_similarity(A,B):
    A[isnan(A)]=0
    B[isnan(B)]=0
    return np.dot(A.flatten(),B.flatten())/(linalg.norm(A.flatten())*linalg.norm(B.flatten()))