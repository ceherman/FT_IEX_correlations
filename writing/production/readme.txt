◦ Figure directories
    ▪ ~/my_work/cadet/aex_flowthrough/images_manuscript
    ▪ ~/my_work/cadet/flowthrough_lysozyme/images/manuscript_images
    ▪ ~/my_work/cadet/flowthrough_FITC_lysozyme/images/manuscript_images
    ▪ ~/my_work/previous/windows_machines/2014_Latitude/Python/AEX pH 7 model parameters and k_prime
    ▪ ~/my_work/correlation_aex_data/images/manuscript_images


    ◦ Retracing the  k’ → K_eq calculations, phase ratio estimation, etc.
        ▪ Relevant folders
            • ~/my_work/previous/windows_machines/workstation/1_my_work/2_AEX_1
            • ~/my_work/previous/windows_machines/workstation/Documents/reference_info
            	pH_7_Keq_data_and_resin_info.xlsx has phase ratios estimated from ISEC data
            • ~/my_work/previous/windows_machines/2014_Latitude/Python/AEX pH 7 model parameters and k_prime
            • ~/my_work/correlation_aex_data
                ◦ aex_pH_7_all_data doesn’t have ADH – CAQ (see instead ~/my_work/previous/windows_machines/2014_Latitude/Python/AEX pH 7 model parameters and k_prime )
                ◦ kprime_data has k’ data, but experimental_data_kprime has dimensional Keq obtained by converting k’ data with phase ratios from ISEC data
                ◦ dimensionless_Keq_1 has dimensionless Keq estimated by converting k’ data
            • ~/my_work/jupyter/exp_data/0_previous
                ◦ Check ./summary_isocratic_data/ for readme on what data were dropped
                ◦ Check ./BMS_mAb_D_retention/
                    ▪ Has baseline interpolation
    ◦ Note:  ~/my_work/jupyter/misc_notebooks/legacy_windows
        ▪ Contains a duplicate of ~/my_work/previous/windows_machines/2014_Latitude/Python
            • AEX pH 7 model parameters and k_prime (modified recently at ~/my_work/previous), Hubbuch and Morbidelli models, lysozyme gradients and AEX pH 9 gradients
    ◦ Note:  conductivity – IS calibrations were at:  ~/my_work/previous/windows_machines/workstation/Documents/1_sp_sepharose_ff

Currently using in updating the the mAb D retention estimates
    /cadet/flowthrough_lysozyme for v_extra
    /jupyter/exp_data/0_previous/BMS_mAb_D_retention for updating estimate with t_extra
    /summary_isocratic_data for converting k' to dimensionless Keq and making table of phase ratios I used
    /correlation_aex_data (in yamamoto_comparison) for updating the Keq and correlatio regression plots

Plots of AEX k' data are at:
~/my_work/previous/windows_machines/2014_Latitude/Python/AEX pH 7 model parameters and k_prime

Plots of meta-correlation - my Keq data -are at:
~/my_work/correlation_aex_data/correlations.ipynb
(Phase ratio data may be here, too)

(Unused)
CADET fits to isocratic AEX data are at:
~/my_work/previous/windows_machines/workstation/1_my_work/2_AEX_1

Note:  check out the 2020-02-26 ACS talk for some neat mean-field electrostatics results
