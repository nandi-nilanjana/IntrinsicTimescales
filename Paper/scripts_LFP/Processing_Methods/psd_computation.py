"""
Run this script after the bipolar LFP computation
This script computes the psd  using psd_array_welch for each individual laminar site at every channel and each area for each epoch
Main function takes inputs - monkey name, epoch name (preSC1/preGo/full_trial), session names


PS : Here I am computing the epoch from preSC1/preGo - 1 seconds instead of 2 sec before preSC1/preGo in neurons. I am avoiding the touch epoch since it 
will be contaminated with beta ?  Since i take 1 sec before preGo i am including pink trials too in the computation which was not the case for neurons 
This is because when i compute 1 sec before preGo, all trials have the direction information since the SC3 offset has happened. But when I do 2 sec before Go,
like neurons, the SC3 onset is also included, so to avoid the bias that neurons have the direction information for blue/green trials but only creating for the 
pink ones


"""

# import packages
import os
import mne
import re
import h5py
from mne.time_frequency import psd_array_welch
import numpy as np
import xarray as xr
#from frites import io
import warnings
# to remove the spam of pandas FutureWarning with iteritems
warnings.simplefilter(action='ignore', category=FutureWarning)

#subfunctions 
def get_mixed_layers(layers_site, depths_site, spacing=400):
    """
    Get the channels which have info from mixed layers. Assuming we are doing a bipolar with spacing
    s, then at a layer change depth_i, we know that depth_i+s/2 and depth_i-s/2 are from different
    layers. This function will give back a vector of the same shape as layers in which the values of
    the mixed channels are marked as 'mixed'.
    :param layers_site: vector with the layers of each channel
    :param depths_site: vector with the depths of each channel
    :param spacing: spacing between the bipolar channels
    ----
    :return:
    mixed_layers: vector with the same shape as layers_site with the values of the mixed channels
    """
    # get the change in layers
    change_layers = (np.where(np.array(layers_site[1:]) != np.array(layers_site[:-1]))[0])
    depth_changes = (np.array(depths_site)[change_layers] + np.array(depths_site)[change_layers+1])/2

    # get the limits of mixed layers
    spacing_mm = spacing/(2*1000)  # in mm and half of the spacing

    # get the mixed layers
    mixed_layers = np.array(layers_site, dtype=object)

    # iterate over the depths and layers
    for i, (layer, depth) in enumerate(zip(layers_site, depths_site)):
        for depth_change in depth_changes:
            if (depth_change - spacing_mm <= depth) & (depth_change + spacing_mm >= depth):
                mixed_layers[i] = 'mixed'
                # in case it already belongs to the first mixed layer, then it should be mixed
                break

    return mixed_layers


def redefine_layers(layers_site):
    """ Redefine the labels of the layers to ignore the split between superficial L23 and deep.
    This way we only have L23, L5 and L6. We know that the 'DEEP' and 'SUP' comes after a '-'.
    :param layers_site: vector with the layers of each channel
    ----
    :return:
        new_layers: vector with the new layers
    """
    new_layers = []
    for layer_i in layers_site:
        new_layers.append(re.split('-', layer_i)[0])

    return new_layers



if __name__ == "__main__":
    
    monkey = 'Tomy'
    psd_freq = [0,100] #frequencies to compute the psd
    epoch = ['preSC1', 'preGo', 'full_trial']
    
    # list all the sessions that we want to preprocess
    if monkey == 'Mourad': 
       LAMINAR_SESS = ['Mo180330001','Mo180405001','Mo180405004','Mo180411001','Mo180412002',
                    'Mo180418002','Mo180419003','Mo180426004','Mo180503002', 'Mo180523002','Mo180524003', 
                    'Mo180525003','Mo180531002','Mo180614002','Mo180614006',
                    'Mo180615002','Mo180615005', 'Mo180619002','Mo180620004','Mo180622002',
                    'Mo180626003', 'Mo180627003','Mo180629005', 'Mo180703003','Mo180704003', 'Mo180705002'
                      'Mo180706002', 'Mo180710002','Mo180711004','Mo180712006']
    elif monkey == 'Tomy':
    
        LAMINAR_SESS =['t140924003','t140925001','t140926002','t140930001','t141001001','t141008001','t141010003','t150122001',
                    't150123001','t150204001','t150205004','t150303002','t150319003','t150320002','t150327002',
                    't150327003','t150415002','t150416002','t150423002','t150430002','t150520003','t150716001','t150702002']
    
    
    for e in epoch:
        for session in LAMINAR_SESS:
            #check where are we running the code
            current_path = os.getcwd()
        
            if current_path.startswith('/Users'):
                server = 'Volumes'  # local w VPN
            elif current_path.startswith('/envau'):
                server = 'envau'  # niolon
            elif current_path.startswith('/hpc/'):
                server = '/envau/work/'  # niolon
                
            # # set the path       
            path = f'/{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/LFP/{monkey}/Bipolar_sites/{session}'
            path_anot = f'/{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/LFP/{monkey}/Unipolar_sites/{session}'
            path_output = f'/{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/LFP/{monkey}/PSD/Epo_{e}/{session}_{psd_freq}'
            
            if not os.path.exists(path_output):
                os.makedirs(path_output)
        
            for i in os.listdir(path):
                if os.path.isfile((os.path.join(path,i))):
                    if f'{session}' in i and '-epo.fif' in i and 'LFP' in i and 'bipolar' in i:
                        file_name= [i]
                        break
                    else:
                        file_name = []
                        
              
            for j in os.listdir(path_anot):
                if os.path.isfile((os.path.join(path_anot,j))):
                    if f'{session}' in j and '-annot.fif' in j and 'LFP' in j:
                        fl_name= [j]
                        break
                    else:
                        fl_name = []
                        
            if not (len(file_name)) | (len(fl_name))>0:
                raise FileNotFoundError(f"No files matching the bipolar session {session} was found.")
        
            # load the data
            LFP_epochs = mne.read_epochs(os.path.join(path, file_name[0]), preload=False)
            
            #load the annotations file for that session to get the Go onset timing 
            annot = mne.read_annotations(os.path.join(path_anot,fl_name[0]))
            times = LFP_epochs.times
            
            #annot.description has the event names and anot.onset the corresponding event onset times
            if e == 'preSC1':
                #get the GO onset time from annotation file
                SC1_time = np.round(annot.onset[annot.description=='SC1'],3)
                idx_SC1 = np.where(times == SC1_time)[0][0]
                idx_before_SC1 =  np.argmin(np.abs(times - (SC1_time - 1)))#find the index of the time when the difference between between the times and the Sc1tim-1 is min
                #not using np.where here since the time SC1_time-1 will never exactly match the time value that i have in the array , there would be some decimal difference
                    
                #extract lfp
                lfp_segment = LFP_epochs.get_data()[:, :, idx_before_SC1:idx_SC1]
                
            elif e == 'preGo':
                #get the GO onset time from annotation file
                GO_time = np.round(annot.onset[annot.description=='GO'],3)
                idx_go = np.where(times == GO_time)[0][0]
                idx_before_go =  np.argmin(np.abs(times - (GO_time - 1)))#find the index of the time when the difference between between the times and the Gotim-1 is min
                #not using np.where here since the time Go_time-1 will never exactly match the time value that i have in the array , there would be some decimal difference
                    
                #extract lfp
                lfp_segment = LFP_epochs.get_data()[:, :, idx_before_go:idx_go]
                
            elif e == 'full_trial':
                #computing over the entire trial that is touch to Go              
                #get touch time
                t_time = np.round(annot.onset[annot.description=='touch'],3)
                
                #get Go time
                GO_time = np.round(annot.onset[annot.description=='GO'],3)
                
                #get the indices
                idx_touch = np.where(times == t_time)[0][0]
                idx_GO = np.where(times == (GO_time - 1))[0][0]
                
                #extract LFP segment
                lfp_segment = LFP_epochs.get_data()[:, :, idx_touch:idx_GO]
            
            # compute the psd on the segmented data 
            n_per_seg = int(LFP_epochs.info['sfreq'])
            psd, freqs = psd_array_welch(lfp_segment, sfreq=LFP_epochs.info['sfreq'],
                                          average='median',fmin = psd_freq[0], fmax = psd_freq[1],
                                          n_per_seg=n_per_seg,
                                          n_overlap=int(n_per_seg/2), n_fft=n_per_seg)       
            # save the data
            # build an xarray with the power in each site
            power_xr = xr.DataArray(psd,
                                    dims=('trials', 'channels', 'freqs'),
                                    coords={'trials': range(LFP_epochs.get_data().shape[0]),
                                            'channels': LFP_epochs.ch_names,
                                            'freqs': freqs})
        
            # check if there are two areas or one
            areas = []
            for area in LFP_epochs.metadata.keys():
                if 'area' in area:
                    areas.append(LFP_epochs.metadata[area][0])
               
            # save the xarray per area (i.e. per site)
            for i_area, area in enumerate(areas):
                area_mask = []                
                for i_ch, channel in enumerate(LFP_epochs.ch_names):
                    if area in channel:
                        area_mask.append(i_ch)
                
                # create the new array
                power_area = power_xr[:, area_mask, :]
                
                layers = LFP_epochs.metadata[f'layers_{+(i_area+1)}'][0] #get the layers from teh metadata
        
                # add the layers and the depths to the xarray. Change the layers to not split L23
                #layers = redefine_layers(LFP_epochs.metadata[f'layers_{+(i_area+1)}'][0]) #i do not need this since I don't split superficial and deep like Laura
                
                depths = LFP_epochs.metadata[f'depths_{+(i_area+1)}'][0]
        
                # get the mixed layers
                new_mix_layers = get_mixed_layers(layers, depths, spacing=400)        
                power_area = power_area.assign_coords(layers=('channels', layers),
                                                      depths=('channels', depths),
                                                      corrected_layers=('channels', new_mix_layers))
        
                # save the xarray
                power_area.to_netcdf(os.path.join(path_output,
                                                  f'{session}-bipolar_PSD-area_{area}-power.nc'))
