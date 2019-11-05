import os

for what in ['wp_mu','wm_mu','wp_e','wm_e']:
    for run in ['all']:
        cmd = 'python photos2pythia8.py {channel}nu_photos.root {channel}nu_pythia8.root {channel}nu_ratio.root --name {channel} --make {run}'.format(channel=what,run=run)
        print cmd
        os.system(cmd)

os.system('hadd -f qed_weights.root qed_weights_wp_mu.root qed_weights_wm_mu.root qed_weights_wp_e.root qed_weights_wm_e.root')
