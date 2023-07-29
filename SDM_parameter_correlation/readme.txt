The first set of work was performed in correlations.ipynb
The yamamoto_comparison.ipynb was created second
The buried_charges.ipynb was created third to identify buried charges (I had already identified 	numbers of surface charges in correlations.ipynb) - In progress
lit_comparison.ipynb was added fourth, after I had already compared to DePhillip's data in 		correlations.ipny and Saleh's gradient data in yamamoto_comparison.ipynb

The aex_pH_7_all_data.xlsx contains dimensionless Keq values fit in CADET, whereas the experimental_data_kprime folder contains Keq data in [m] that were obtained by converting k' values to Keq using a dimensional phase ratio (based on literature data with the cylindrical pore model). Apologies for the poor nomenclature; the folder is called kprime, but it's actually just Keq data.

The kprime_data folder actually contains experimental kprime data.

The ns_bim_spherical_pqr > Keq_predicitons directory contains dimensional Keq [m] predicted from the numerical Derjaguin model (that uses Qiang's NS-BIM code).

experimental_Keq_IS_correlation_params.csv contains a and b parameters fit to Keq = a (IS)^b for the (dimensional) experimental data.

ns_bim_spherical_pqr contains the folders charges and Keq_predictions. Keq_predictions is as explained above. The charges folder are the charges used in Qiang's code (i.e. discrete charge states extracted from .pqr files and projected onto spheres.

pdb_structure_files just contains pdb files - I don't believe I used these in this analysis.

pqr_structure_files contains pqr.pqr files that were generated with APBS at some point previously (see the retro description on 2021-01-06 in my computational notebook). The other files, including output.txt, pqr_atom.dep, and pqr_res.dep were generated with EDTSurf. These contain information about the depth of atoms beneath the protein surface and the protein surface area (the output.txt file is just the shell output of edtsurf)

In the images folder, models are referred to as linear, because I fit a 2-parameter model to the data. Keq_dimensional_lin_correlation corresponds to the experimental Keq from k' data, Keq_dimensionless_data_lin_correlation corresponds to Keq values fit in CADET, Keq_pred_lin_correlation corresponds to dimensional Keq predicted with the numerical Derjaguin model, meta_lin_correlation contains meta correlations, and protein characteristics has charge numbers vs depth beneath protein surface. Note the jupyter notebook saves some plots just to the images folder, and I manually sorted them.


The lit_data folder contains data from DePhillips' "Determinants of protein retention characteristics on cation-exchange adsorbents" and Saleh's "Straightforward method for calibration of mechanistic cation exchange chromatography models for industrial applications"

The lysozyme_seph_gradient_GH_results contains log10(GH [mM])-log10(IS [mM]) I obtained from linear gradient elution of lysozyme on SP Sepharose FF. I used the total phase ratio in the GH computation.

Used in yamamoto_comparison.ipynb___________________________________________________

dimensionless_Keq_1 contains estimates for dimensionless Keq values obtained by converting k' to Keq with solute-resin-dependent phase ratios (i.e. (1-eps_t_protein_resin)/eps_t_protein_resin). solute-resin-dependent total porosities were estimated with lysosyme measurments for all lysozyme data, fit column and particle porosities for AEX data, and need to be remeasured for mAb D. These files were generated in extra_ubuntu_spa/summary_isocratic_data/


Comparison to lit data________________________________________________________
Started with comparison to DePhillips' Keq [m] in correlations.ipynb and Saleh's gradient data in yamamoto_comparison.ipynb
Then went into lit_comparison.ipynb
General notes:  x in the DePhillips cleaned data is the NaCl concentraiton [M]
The DePhillips phase ratio .csv file contains phase ratios in m2/ml (see his table 4)
Pulled and used Staby data, in k'
Correlation holds up generally with all data, but it's cleaner for individual data sets
I performed a multiple linear regression with of the (a, b) data on my laptop in Minitab (Documents > Research > misc)
	The protein explains more of the variability in the data than the resin, although both 		are significant.

Noise problem________________________________________________________________________
Identification of 95% confidence interval and prediction interval for augmented k' data set

(95% prediction interval estimates on ln(a) propagated directly to a for constructing the 
bands on the k'-IS plots)

Conversion of electrostatics model prediction Keq to k' and comparison of meta-correlation
with that of the augmented data set. Only data shown are for POROS resins (the phase ratio is
unknown for Capto Q, and I didn't have the Keq estimates for lysozyme on hand)


























































