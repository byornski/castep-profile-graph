+-----------------------------------------------------------------------------+
|     Subroutine                              Total      Profiled   Time      |
|                                             calls      calls      (incl.)   |
+-----------------------------------------------------------------------------+
|   o-- <parent(s) not traced>                      1                         |
|  /                                                                          |
| O-> castep                                        1           1    197.58s  |
|    /                                                                        |
|   o-> memory_system_initialise                    1           1      0.00s  |
|   o-> comms_gcopy_real                            1           1      0.00s  |
|   o-> cell_read_wrapped                           1           1      0.13s  |
|   o-> ion_read                                    1           1      0.06s  |
|   o-> parameters_read                             1           1      0.11s  |
|   o-> bib_add                                     5           5      0.00s  |
|   o-> tddft_set_tddft_on                          1           1      0.00s  |
|   o-> comms_parallel_strategy                     1           1      0.00s  |
|   o-> cell_distribute_kpoints_wrapped             1           1      0.00s  |
|   o-> ion_initialise                              1           1      0.30s  |
|   o-> basis_initialise                            1           1      0.35s  |
|   o-> ion_real_initialise                         1           1      0.00s  |
|   o-> model_initialise                            1           1      0.09s  |
|   o-> nlxc_initialise                             1           1      0.00s  |
|   o-> parameters_output                           1           1      0.01s  |
|   o-> cell_output_wrapped                         1           1      0.01s  |
|   o-> dftd_sedc_initialize                        1           1      0.00s  |
|   o-> check_elec_ground_state                     1           1    191.16s  |
|   o-> check_forces_stresses                       1           1      1.92s  |
|   o-> popn_calculate_mulliken                     1           1      1.51s  |
|   o-> write_eigenvalues                           1           1      1.48s  |
|   o-> write_local_pot_and_den                     1           1      0.03s  |
|   o-> hirshfeld_calculate                         1           1      0.22s  |
|   o-> hirshfeld_output                            1           1      0.00s  |
|   o-> model_write                                 1           1      0.14s  |
|   o-> bib_output                                  1           1      0.00s  |
|   o-> model_deallocate                            1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> check_elec_ground_state                       1           1    191.16s  |
|    /                                                                        |
|   o-> castep_calc_storage                         1           1      0.01s  |
|   o-> phonon_store_mem_estimate                   3           3      0.00s  |
|   o-> electronic_calc_storage                     1           1      0.00s  |
|   o-> firstd_calc_storage                         2           2      0.00s  |
|   o-> castep_report_storage                       1           1      0.00s  |
|   o-> electronic_minimisation                     1           1    190.83s  |
|   o-> model_write                                 1           1      0.32s  |
+-----------------------------------------------------------------------------+
|   o-- check_elec_ground_state                     1                         |
|  /                                                                          |
| O-> electronic_minimisation                       1           1    190.83s  |
|    /                                                                        |
|   o-> comms_gcopy_logical                         1           1      0.00s  |
|   o-> comms_gcopy_character                       1           1      0.00s  |
|   o-> electronic_initialise                       1           1      0.04s  |
|   o-> electronic_prepare_H                       10          10      0.68s  |
|   o-> electronic_apply_H_energy_local             1           1      0.00s  |
|   o-> electronic_write_energies                  12          12      0.00s  |
|   o-> electronic_store_energy                    12          12      0.00s  |
|   o-> electronic_scf_banner                       1           1      0.00s  |
|   o-> electronic_write_scf_energies              12          12      0.00s  |
|   o-> hamiltonian_diagonalise_ks                 11          11      8.34s  |
|   o-> electronic_find_occupancies                11          11      0.14s  |
|   o-> electronic_check_occupancies               11          11      0.00s  |
|   o-> electronic_apply_H_energy_eigen            11          11      0.30s  |
|   o-> electronic_write_spin_density              11          11      0.00s  |
|   o-> electronic_dump                            11          11      0.00s  |
|   o-> hubbard_calculate_occ_matrix                9           9      0.00s  |
|   o-> hubbard_mix_occ_matrix                      9           9      0.00s  |
|   o-> density_calculate_soft_wvfn                 9           9      0.15s  |
|   o-> density_symmetrise                          9           9      0.00s  |
|   o-> density_augment                             9           9      0.14s  |
|   o-> dm_mix_density                              9           9      0.01s  |
|   o-> ewald_dipole_corr                           1           1      0.01s  |
|   o-> electronic_scf_footer                       1           1      0.00s  |
|   o-> electronic_write_occupancies                1           1      0.00s  |
|   o-> dftd_sedc_calculate_energy_cell             1           1    181.02s  |
|   o-> dftd_sedc_print_corr_energies               1           1      0.00s  |
|   o-> electronic_finalise                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_sedc_calculate_energy_cell             1                         |
|  /o-- dftd_sedc_calculate_forces_cell             1                         |
|  /o-- dftd_sedc_calculate_stress_cell             1                         |
| |/                                                                          |
| O-> dftd_sedc_calculate_cell                      3           3    181.03s  |
|    /                                                                        |
|   o-> dftd_sedc_update_cell                       3           3      0.02s  |
|   o-> comms_gcopy_logical                         4           4      0.00s  |
|   o-> dftd_sedc_update_TS                         1           1      0.22s  |
|   o-> sedc                                        1           1    180.79s  |
|   o-> comms_gcopy_integer                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> dftd_sedc_calculate_energy_cell               1           1    181.02s  |
|    /                                                                        |
|   o-> dftd_sedc_calculate_cell                    1           1    181.02s  |
|   o-> comms_gcopy_real                            1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_sedc_calculate_cell                    1                         |
|  /                                                                          |
| O-> sedc                                          1           1    180.79s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    11                         |
|  /                                                                          |
| O-> hamiltonian_diagonalise_ks                   11          11      8.34s  |
|    /                                                                        |
|   o-> wave_allocate_slice                        55          55      0.00s  |
|   o-> wave_initialise_slice                     253         253      0.06s  |
|   o-> comms_reduce_gv_logical                    11          11      0.00s  |
|   o-> comms_reduce_bnd_logical                  192         192      0.00s  |
|   o-> hamiltonian_apply_ks                        9           9      1.42s  |
|   o-> wave_diagonalise_H_ks                      11          11      0.32s  |
|   o-> comms_reduce_bnd_real                      11          11      0.00s  |
|   o-> wave_calc_precon                           11          11      0.00s  |
|   o-> nlpot_prepare_precon_ks                    11          11      0.30s  |
|   o-> comms_reduce_bnd_integer                  697         697      0.00s  |
|   o-> wave_copy_wv_slice                        154         154      0.02s  |
|   o-> wave_copy_slice_slice                    3542        3542      0.03s  |
|   o-> wave_dot_all_slice_slice                  258         258      0.16s  |
|   o-> hamiltonian_searchspace_ks                181         181      3.17s  |
|   o-> wave_Sorthog_lower_slice_s                181         181      0.14s  |
|   o-> wave_Sorthonormalise_slice                181         181      0.07s  |
|   o-> hamiltonian_apply_slice                   181         181      1.60s  |
|   o-> comms_reduce_bnd_complex                  181         181      0.03s  |
|   o-> algor_diagonalise_complex                 181         181      0.19s  |
|   o-> wave_rotate_slice                         362         362      0.65s  |
|   o-> wave_copy_slice_wv                        362         362      0.03s  |
|   o-> comms_copy_gv_logical                     543         543      0.00s  |
|   o-> wave_deallocate_slice                      55          55      0.00s  |
|   o-> wave_kinetic_eigenvalues_wv_ks              2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_all_beta_multi_phi_recip              376                         |
|  /o-- ion_beta_add_multi_recip_all              704                         |
|  /o-- local_dot_all_many_many                   632                         |
|  /o-- local_q_dot_all_many_many_c               363                         |
|  /o-- ion_Q_apply_and_add_wsf                    90                         |
| |/                                                                          |
| O-> algor_matmul_cmplx_cmplx                   2165        2165      4.28s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                181                         |
|  /                                                                          |
| O-> hamiltonian_searchspace_ks                  181         181      3.17s  |
|    /                                                                        |
|   o-> wave_initialise_slice                     181         181      0.04s  |
|   o-> nlpot_apply_precon_ES_slice               181         181      1.80s  |
|   o-> wave_Sorthogonalise_wv_slice              181         181      1.33s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_apply_add_slice_r                   523                         |
|  /o-- nlpot_apply_precon_ES_slice               181                         |
| |/                                                                          |
| O-> ion_beta_add_multi_recip_all                704         704      2.35s  |
|    /                                                                        |
|   o-> ion_beta_recip_set                        704         704      0.01s  |
|   o-> algor_matmul_cmplx_cmplx                  704         704      2.33s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> check_forces_stresses                         1           1      1.92s  |
|    /                                                                        |
|   o-> firstd_calculate_forces                     1           1      0.56s  |
|   o-> firstd_symmetrise_forces                    1           1      0.00s  |
|   o-> firstd_output_forces                        1           1      0.01s  |
|   o-> firstd_calculate_stress                     1           1      1.35s  |
|   o-> firstd_symmetrise_stress                    1           1      0.00s  |
|   o-> firstd_output_stress                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_beta_phi_wv_ks                        14                         |
|  /o-- wave_beta_phi_slice                       362                         |
| |/                                                                          |
| O-> ion_all_beta_multi_phi_recip                376         376      1.86s  |
|    /                                                                        |
|   o-> ion_beta_recip_set                        376         376      0.02s  |
|   o-> algor_matmul_cmplx_cmplx                  376         376      1.53s  |
|   o-> comms_reduce_gv_complex                   376         376      0.28s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_searchspace_ks                181                         |
|  /                                                                          |
| O-> nlpot_apply_precon_ES_slice                 181         181      1.80s  |
|    /                                                                        |
|   o-> wave_beta_phi_slice                       362         362      0.83s  |
|   o-> comms_reduce_gv_complex                   181         181      0.05s  |
|   o-> ion_beta_add_multi_recip_all              181         181      0.66s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                      342                         |
|  /o-- hamiltonian_apply_slice                   181                         |
| |/                                                                          |
| O-> nlpot_apply_add_slice                       523         523      1.79s  |
|    /                                                                        |
|   o-> nlpot_apply_add_slice_r                   523         523      1.79s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_apply_add_slice                     523                         |
|  /                                                                          |
| O-> nlpot_apply_add_slice_r                     523         523      1.79s  |
|    /                                                                        |
|   o-> wave_spin_type_slice                      523         523      0.00s  |
|   o-> wave_beta_phi_slice                       523         523      0.01s  |
|   o-> comms_reduce_gv_complex                   523         523      0.07s  |
|   o-> ion_beta_add_multi_recip_all              523         523      1.69s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                181                         |
|  /                                                                          |
| O-> hamiltonian_apply_slice                     181         181      1.60s  |
|    /                                                                        |
|   o-> wave_initialise_slice                     181         181      0.03s  |
|   o-> pot_apply_slice                           181         181      0.77s  |
|   o-> nlpot_apply_add_slice                     181         181      0.77s  |
|   o-> wave_add_kinetic_energy_slice             181         181      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_apply_add_slice_r                   523                         |
|  /o-- nlpot_apply_precon_ES_slice               362                         |
|  /o-- wave_coeffs_Sdot_all_wv_slice             181                         |
|  /o-- wave_Sdot_lower_slice                     362                         |
|  /o-- wave_calc_Soverlap_slice                  181                         |
|  /o-- wave_orthonormalise_over_slice            181                         |
| |/                                                                          |
| O-> wave_beta_phi_slice                        1790        1790      1.55s  |
|    /                                                                        |
|   o-> ion_set_projectors                       1790        1790      0.02s  |
|   o-> wave_calc_ps_q_nonzero                   1790        1790      0.00s  |
|   o-> ion_all_beta_multi_phi_recip              362         362      1.49s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> popn_calculate_mulliken                       1           1      1.51s  |
|    /                                                                        |
|   o-> popn_create_lcao_basis_allkpt               1           1      0.58s  |
|   o-> popn_calculate_matrices                     1           1      0.10s  |
|   o-> wave_deallocate_wv                          1           1      0.00s  |
|   o-> popn_calculate_charge_spilling              1           1      0.07s  |
|   o-> popn_calculate_density_matrix               1           1      0.01s  |
|   o-> popn_calculate_populations                  1           1      0.76s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> write_eigenvalues                             1           1      1.48s  |
|    /                                                                        |
|   o-> comms_reduce_bnd_real                       1           1      0.00s  |
|   o-> comms_gather_kp_integer                     1           1      0.00s  |
|   o-> global_kpoint_index                         4           4      0.00s  |
|   o-> comms_recv_real                             9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                  9                         |
|  /                                                                          |
| O-> hamiltonian_apply_ks                          9           9      1.42s  |
|    /                                                                        |
|   o-> wave_allocate_slice                        18          18      0.00s  |
|   o-> wave_initialise_slice                      18          18      0.00s  |
|   o-> wave_copy_wv_slice                        342         342      0.00s  |
|   o-> pot_apply_slice                           342         342      0.38s  |
|   o-> nlpot_apply_add_slice                     342         342      1.02s  |
|   o-> wave_copy_slice_wv                        342         342      0.00s  |
|   o-> wave_add_kinetic_energy_wv_ks               9           9      0.01s  |
|   o-> wave_deallocate_slice                      18          18      0.00s  |
|   o-> wave_dot_wv_wv_ks                           9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- check_forces_stresses                       1                         |
|  /                                                                          |
| O-> firstd_calculate_stress                       1           1      1.35s  |
|    /                                                                        |
|   o-> ewald_calculate_stress                      1           1      0.11s  |
|   o-> wave_kinetic_stress_wv                      1           1      0.00s  |
|   o-> locpot_calculate_stress                     1           1      0.04s  |
|   o-> hubbard_calculate_stress                    1           1      0.00s  |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> locpot_calculate                            1           1      0.01s  |
|   o-> nlpot_calculate_stress                      1           1      1.19s  |
|   o-> pot_deallocate                              1           1      0.00s  |
|   o-> dftd_sedc_calc_stress_curr_cell             1           1      0.01s  |
|   o-> tddft_calculate_stress                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_searchspace_ks                181                         |
|  /                                                                          |
| O-> wave_Sorthogonalise_wv_slice                181         181      1.33s  |
|    /                                                                        |
|   o-> wave_Sdot_all_wv_slice                    181         181      1.06s  |
|   o-> wave_orthog_over_wv_slice                 181         181      0.26s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> nlpot_calculate_stress                        1           1      1.19s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                           1           1      0.00s  |
|   o-> nlpot_calculate_stress_r                    1           1      1.19s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_stress                      1                         |
|  /                                                                          |
| O-> nlpot_calculate_stress_r                      1           1      1.19s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                           1           1      0.00s  |
|   o-> wave_beta_phi_wv                            1           1      0.00s  |
|   o-> nlpot_allocate_nl_d                         1           1      0.00s  |
|   o-> nlpot_calculate_d                           1           1      0.05s  |
|   o-> wave_dbetadcell_phi_wv_iks                160         160      0.25s  |
|   o-> nlpot_calculate_dddcell_real                1           1      0.84s  |
|   o-> comms_reduce_gv_real                        1           1      0.00s  |
|   o-> comms_reduce_bnd_real                       1           1      0.00s  |
|   o-> comms_reduce_kp_real                        1           1      0.05s  |
|   o-> nlpot_deallocate_nl_d                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                      342                         |
|  /o-- hamiltonian_apply_slice                   181                         |
| |/                                                                          |
| O-> pot_apply_slice                             523         523      1.15s  |
|    /                                                                        |
|   o-> pot_debug_check                           523         523      0.00s  |
|   o-> pot_interpolate                             9           9      0.01s  |
|   o-> wave_spin_type_slice                      523         523      0.00s  |
|   o-> pot_nongamma_apply_slice                  523         523      1.14s  |
+-----------------------------------------------------------------------------+
|   o-- pot_apply_slice                           523                         |
|  /                                                                          |
| O-> pot_nongamma_apply_slice                    523         523      1.14s  |
|    /                                                                        |
|   o-> wave_spin_type_slice                      523         523      0.00s  |
|   o-> wave_recip_to_real_slice                 1151        1151      0.56s  |
|   o-> wave_real_to_recip_slice_slice           1151        1151      0.48s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthogonalise_wv_slice              181                         |
|  /                                                                          |
| O-> wave_Sdot_all_wv_slice                      181         181      1.06s  |
|    /                                                                        |
|   o-> wave_coeffs_Sdot_all_wv_slice             181         181      0.99s  |
|   o-> comms_reduce_gv_complex                   181         181      0.07s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sdot_all_wv_slice                    181                         |
|  /                                                                          |
| O-> wave_coeffs_Sdot_all_wv_slice               181         181      0.99s  |
|    /                                                                        |
|   o-> coeffs_dot_all_many_many                  181         181      0.13s  |
|   o-> wave_beta_phi_wv_ks                       181         181      0.02s  |
|   o-> wave_beta_phi_slice                       181         181      0.70s  |
|   o-> wave_q_dot_all_many_many_c                181         181      0.14s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_stress_r                    1                         |
|  /                                                                          |
| O-> nlpot_calculate_dddcell_real                  1           1      0.84s  |
|    /                                                                        |
|   o-> basis_real_to_recip_grid                    1           1      0.00s  |
|   o-> basis_recip_fine_to_half_grid               1           1      0.00s  |
|   o-> ion_int_dQdcell_at_origin_recip          1536        1536      0.82s  |
|   o-> comms_reduce_gv_real                        1           1      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              2                         |
|  /o-- ion_generate_orbitals                       2                         |
| |/                                                                          |
| O-> ion_atom_pseudo_scf                           4           4      0.83s  |
|    /                                                                        |
|   o-> ion_atom_init_pseudo_basis                  4           4      0.01s  |
|   o-> ion_atom_init_pseudo_atom                   4           4      0.09s  |
|   o-> ion_atom_init_pseudo_H                      4           4      0.00s  |
|   o-> ion_atom_ps_diag                           62          62      0.53s  |
|   o-> ion_atom_set_pseudo_H                      62          62      0.02s  |
|   o-> ion_atom_regin                          14007       14007      0.02s  |
|   o-> ion_atom_basis_pseudo_dealloc               4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_dddcell_real             1536                         |
|  /                                                                          |
| O-> ion_int_dQdcell_at_origin_recip            1536        1536      0.82s  |
|    /                                                                        |
|   o-> ion_set_dQdcell_at_origin_recip          1536        1536      0.78s  |
+-----------------------------------------------------------------------------+
|   o-- basis_recip_grid_to_real                  820                         |
|  /o-- basis_real_to_recip_grid                  192                         |
|  /o-- basis_real_fine_to_std_grid                36                         |
|  /o-- basis_recip_red_real_3d_coeffs           9074                         |
|  /o-- basis_real_recip_red_3d_coeffs           9074                         |
|  /o-- basis_recip_reduced_real_one             2338                         |
|  /o-- basis_real_std_to_fine_grid                36                         |
| |/                                                                          |
| O-> comms_transpose                           21570       21570      0.80s  |
|    /                                                                        |
|   o-> comms_transpose_n                       21570       21570      0.74s  |
+-----------------------------------------------------------------------------+
|   o-- ion_int_dQdcell_at_origin_recip          1536                         |
|  /                                                                          |
| O-> ion_set_dQdcell_at_origin_recip            1536        1536      0.78s  |
|    /                                                                        |
|   o-> ion_atom_radial_transform                1536        1536      0.00s  |
|   o-> ion_Q_recip_interpolation                1824        1824      0.22s  |
|   o-> ion_apply_and_add_ylm                   31680       31680      0.16s  |
|   o-> ion_apply_and_add_ylmp                  31680       31680      0.15s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_mulliken                     1                         |
|  /                                                                          |
| O-> popn_calculate_populations                    1           1      0.76s  |
|    /                                                                        |
|   o-> comms_reduce_kp_complex                     1           1      0.00s  |
|   o-> cell_format_atom_label                  24112       24112      0.11s  |
|   o-> cell_verify_mixture_components          25440       25440      0.11s  |
|   o-> comms_reduce_kp_real                        1           1      0.02s  |
|   o-> popn_sort                                   1           1      0.07s  |
+-----------------------------------------------------------------------------+
|   o-- wave_rotate_wv_ks                          33                         |
|  /o-- wave_orthog_over_wv_slice                 362                         |
|  /o-- wave_orthog_lower_over_slice_s            362                         |
|  /o-- wave_rotate_slice                         543                         |
| |/                                                                          |
| O-> algor_matmul_cmplx_cmplx_3D                1300        1300      0.74s  |
+-----------------------------------------------------------------------------+
|   o-- comms_transpose                         21570                         |
|  /                                                                          |
| O-> comms_transpose_n                         21570       21570      0.74s  |
|    /                                                                        |
|   o-> comms_transpose_exchange                21570       21570      0.67s  |
+-----------------------------------------------------------------------------+
|   o-- comms_reduce_gv_complex                  2243                         |
|  /o-- comms_reduce_bnd_complex                  214                         |
|  /o-- comms_reduce_kp_complex                    10                         |
| |/                                                                          |
| O-> comms_reduce_array_complex                 2467        2467      0.73s  |
+-----------------------------------------------------------------------------+
|   o-- ion_all_beta_multi_phi_recip              376                         |
|  /o-- wave_calc_Soverlap_wv_ks                    2                         |
|  /o-- nlpot_apply_add_slice_r                   523                         |
|  /o-- wave_dot_wv_wv_ks                           9                         |
|  /o-- wave_dot_all_wv_wv_ks                      11                         |
|  /o-- ion_beta_beta_recip_cmplx                   5                         |
|  /o-- nlpot_prepare_precon_ks                     5                         |
|  /o-- wave_dot_all_slice_slice                  258                         |
|  /o-- nlpot_apply_precon_ES_slice               181                         |
|  /o-- wave_Sdot_all_wv_slice                    181                         |
|  /o-- wave_Sdot_lower_slice                     181                         |
|  /o-- wave_calc_Soverlap_slice                  181                         |
|  /o-- ion_augment_charge_nospin_kp                9                         |
|  /o-- dm_mix_density_dot                        165                         |
|  /o-- basis_sum_recip_grid                      492                         |
|  /o-- wave_dbetadR_phi_wv                       160                         |
|  /o-- wave_dbetadcell_phi_wv_iks                160                         |
|  /o-- wave_Sdot_all_wv_wv_ks                      1                         |
| |/                                                                          |
| O-> comms_reduce_gv_complex                    2900        2900      0.69s  |
|    /                                                                        |
|   o-> comms_reduce_array_complex               2243        2243      0.67s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    10                         |
|  /                                                                          |
| O-> electronic_prepare_H                         10          10      0.68s  |
|    /                                                                        |
|   o-> density_symmetrise                         10          10      0.00s  |
|   o-> locpot_calculate                           10          10      0.12s  |
|   o-> nlpot_calculate_d                          10          10      0.56s  |
+-----------------------------------------------------------------------------+
|   o-- comms_transpose_n                       21570                         |
|  /                                                                          |
| O-> comms_transpose_exchange                  21570       21570      0.67s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_prepare_H                       10                         |
|  /o-- nlpot_calculate_forces_r                    1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
| |/                                                                          |
| O-> nlpot_calculate_d                            12          12      0.66s  |
|    /                                                                        |
|   o-> nlpot_calculate_d_real                     12          12      0.66s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_d                          12                         |
|  /                                                                          |
| O-> nlpot_calculate_d_real                       12          12      0.66s  |
|    /                                                                        |
|   o-> basis_real_to_recip_grid                   12          12      0.00s  |
|   o-> basis_recip_fine_to_half_grid              12          12      0.00s  |
|   o-> ion_int_Q_at_origin_recip               18432       18432      0.50s  |
|   o-> comms_reduce_gv_real                       12          12      0.00s  |
|   o-> comms_reduce_bnd_real                      12          12      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_d_real                  18432                         |
|  /o-- nlpot_calculate_dddR                     4608                         |
| |/                                                                          |
| O-> ion_int_Q_at_origin_recip                 23040       23040      0.66s  |
|    /                                                                        |
|   o-> ion_set_Q_at_origin_recip               23040       23040      0.16s  |
|   o-> algor_re_dot_cmplx_cmplx                23040       23040      0.21s  |
|   o-> comms_reduce_gv_real                     4608        4608      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                362                         |
|  /                                                                          |
| O-> wave_rotate_slice                           362         362      0.65s  |
|    /                                                                        |
|   o-> wave_allocate_slice                       362         362      0.00s  |
|   o-> wave_initialise_slice                     362         362      0.27s  |
|   o-> algor_matmul_cmplx_cmplx_3D               543         543      0.33s  |
|   o-> wave_copy_slice_slice                     362         362      0.02s  |
|   o-> wave_deallocate_slice                     362         362      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_mulliken                     1                         |
|  /                                                                          |
| O-> popn_create_lcao_basis_allkpt                 1           1      0.58s  |
|    /                                                                        |
|   o-> popn_check_lcao_basis                       1           1      0.00s  |
|   o-> ion_generate_orbitals                       2           2      0.57s  |
|   o-> wave_allocate_wv                            1           1      0.00s  |
|   o-> wave_initialise_wv                          1           1      0.00s  |
|   o-> basis_radial_to_recip_reduced             256         256      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_create_lcao_basis_allkpt               2                         |
|  /                                                                          |
| O-> ion_generate_orbitals                         2           2      0.57s  |
|    /                                                                        |
|   o-> ion_atom_allocate_pspot                     2           2      0.00s  |
|   o-> ion_set_psp                                 2           2      0.00s  |
|   o-> ion_atom_pseudo_scf                         2           2      0.57s  |
|   o-> ion_atom_deallocate_pspot                   2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_nongamma_apply_slice                 1151                         |
|  /                                                                          |
| O-> wave_recip_to_real_slice                   1151        1151      0.56s  |
|    /                                                                        |
|   o-> basis_recip_red_real_3d_coeffs           1151        1151      0.50s  |
+-----------------------------------------------------------------------------+
|   o-- check_forces_stresses                       1                         |
|  /                                                                          |
| O-> firstd_calculate_forces                       1           1      0.56s  |
|    /                                                                        |
|   o-> ewald_calculate_forces                      1           1      0.08s  |
|   o-> locpot_calculate_forces                     1           1      0.03s  |
|   o-> hubbard_calculate_forces                    1           1      0.00s  |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> locpot_calculate                            1           1      0.01s  |
|   o-> nlpot_calculate_forces                      1           1      0.42s  |
|   o-> pot_deallocate                              1           1      0.00s  |
|   o-> ewald_dipole_corr                           1           1      0.00s  |
|   o-> dftd_sedc_calc_forces_curr_cell             1           1      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_pseudo_scf                        62                         |
|  /                                                                          |
| O-> ion_atom_ps_diag                             62          62      0.53s  |
|    /                                                                        |
|   o-> ion_atom_regin                          54675       54675      0.24s  |
|   o-> ion_atom_rectoreal                        245         245      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_recip_to_real_slice                 1151                         |
|  /                                                                          |
| O-> basis_recip_red_real_3d_coeffs             1151        1151      0.50s  |
|    /                                                                        |
|   o-> comms_transpose                          9074        9074      0.33s  |
+-----------------------------------------------------------------------------+
|   o-- pot_nongamma_apply_slice                 1151                         |
|  /                                                                          |
| O-> wave_real_to_recip_slice_slice             1151        1151      0.48s  |
|    /                                                                        |
|   o-> basis_real_recip_red_3d_coeffs           1151        1151      0.47s  |
+-----------------------------------------------------------------------------+
|   o-- wave_real_to_recip_slice_slice           1151                         |
|  /                                                                          |
| O-> basis_real_recip_red_3d_coeffs             1151        1151      0.47s  |
|    /                                                                        |
|   o-> comms_transpose                          9074        9074      0.31s  |
+-----------------------------------------------------------------------------+
|   o-- check_elec_ground_state                     1                         |
|  /o-- castep                                      1                         |
| |/                                                                          |
| O-> model_write                                   2           2      0.46s  |
|    /                                                                        |
|   o-> model_write_all                             4           4      0.45s  |
|   o-> comms_barrier_farm                          2           2      0.00s  |
|   o-> comms_barrier                               2           2      0.00s  |
|   o-> io_delete_file                              1           1      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_sedc_update_TS                         1                         |
|  /o-- castep                                      1                         |
| |/                                                                          |
| O-> hirshfeld_calculate                           2           2      0.45s  |
|    /                                                                        |
|   o-> hirshfeld_init                              2           2      0.00s  |
|   o-> cell_allocate                               8           8      0.00s  |
|   o-> cell_copy                                   8           8      0.00s  |
|   o-> comms_gcopy_integer                         2           2      0.00s  |
|   o-> comms_gcopy_logical                         4           4      0.00s  |
|   o-> cell_supercell                              2           2      0.02s  |
|   o-> cell_distribute_kpoints_wrapped             2           2      0.00s  |
|   o-> cell_deallocate                             6           6      0.00s  |
|   o-> basis_radial_to_recip_grid                  4           4      0.00s  |
|   o-> basis_recip_grid_to_real                  320         320      0.07s  |
|   o-> comms_reduce_gv_real                     1282        1282      0.00s  |
|   o-> comms_gcopy_real                            4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_write                                 4                         |
|  /                                                                          |
| O-> model_write_all                               4           4      0.45s  |
|    /                                                                        |
|   o-> parameters_dump                             4           4      0.02s  |
|   o-> cell_dump                                   8           8      0.03s  |
|   o-> model_write_occ_eigenvalues                 4           4      0.00s  |
|   o-> density_write                               4           4      0.19s  |
|   o-> wave_write_all                              2           2      0.11s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                253                         |
|  /o-- hamiltonian_apply_ks                       18                         |
|  /o-- wave_rotate_wv_ks                          22                         |
|  /o-- hamiltonian_searchspace_ks                181                         |
|  /o-- hamiltonian_apply_slice                   181                         |
|  /o-- wave_rotate_slice                         362                         |
|  /o-- electronic_permute_eigenstates             44                         |
| |/                                                                          |
| O-> wave_initialise_slice                      1061        1061      0.43s  |
|     +- section "initialisation", branch on value of method:-                |
|       z                                         880         880      0.35s  |
|         \                                                                   |
|          o-> wave_zero_slice                    880         880      0.35s  |
|       Z                                         181         181      0.03s  |
|         \                                                                   |
|   o-> wave_zero_slice                           181         181      0.03s  |
|    /                                                                        |
|   o-> wave_setup_slice                         1061        1061      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_forces                     1                         |
|  /                                                                          |
| O-> nlpot_calculate_forces                        1           1      0.42s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                           1           1      0.00s  |
|   o-> nlpot_calculate_forces_r                    1           1      0.42s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_forces                      1                         |
|  /                                                                          |
| O-> nlpot_calculate_forces_r                      1           1      0.42s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                           1           1      0.00s  |
|   o-> wave_beta_phi_wv                            1           1      0.00s  |
|   o-> nlpot_allocate_nl_d                         1           1      0.00s  |
|   o-> nlpot_calculate_d                           1           1      0.05s  |
|   o-> wave_dbetadR_phi_wv                       160         160      0.15s  |
|   o-> nlpot_calculate_dddR                        1           1      0.19s  |
|   o-> comms_reduce_gv_real                        1           1      0.00s  |
|   o-> comms_reduce_bnd_real                       1           1      0.00s  |
|   o-> comms_reduce_kp_real                      160         160      0.02s  |
|   o-> nlpot_deallocate_nl_d                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_calc_Soverlap_wv_ks                    2                         |
|  /o-- nlpot_calc_eigenvals_nkns_r                12                         |
|  /o-- wave_coeffs_Sdot_all_wv_slice             181                         |
|  /o-- wave_beta_phi_wv                            2                         |
|  /o-- wave_coeffs_Sdot_all_wv_ks                  2                         |
| |/                                                                          |
| O-> wave_beta_phi_wv_ks                         199         199      0.40s  |
|    /                                                                        |
|   o-> ion_set_projectors                        199         199      0.02s  |
|   o-> wave_calc_ps_q_nonzero                     14          14      0.00s  |
|   o-> ion_all_beta_multi_phi_recip               14          14      0.37s  |
+-----------------------------------------------------------------------------+
|   o-- wave_initialise_slice                    1061                         |
|  /                                                                          |
| O-> wave_zero_slice                            1061        1061      0.37s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> basis_initialise                              1           1      0.35s  |
|    /                                                                        |
|   o-> cell_copy                                   1           1      0.00s  |
|   o-> basis_utils_prime_factors                   6           6      0.00s  |
|   o-> basis_distribute_grids                      1           1      0.30s  |
|   o-> basis_map_standard_to_fine                  1           1      0.00s  |
|   o-> basis_map_fine_recip_half_full              1           1      0.00s  |
|   o-> basis_assign_grid_coordinates               1           1      0.00s  |
|   o-> basis_count_plane_waves                     1           1      0.04s  |
|   o-> basis_assign_plane_wave_indexes             1           1      0.00s  |
|   o-> basis_assign_pw_gvectors                    1           1      0.00s  |
|   o-> basis_calculate_cut_off                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                 11                         |
|  /                                                                          |
| O-> wave_diagonalise_H_ks                        11          11      0.32s  |
|    /                                                                        |
|   o-> wave_dot_all_wv_wv_ks                      11          11      0.06s  |
|   o-> algor_diagonalise_complex                  11          11      0.09s  |
|   o-> wave_rotate_wv_ks                          22          22      0.17s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_distribute_grids                        1           1      0.30s  |
|    /                                                                        |
|   o-> basis_utils_sort_columns                    1           1      0.00s  |
|   o-> comms_gather_gv_integer                    12          12      0.00s  |
|   o-> comms_gather_kp_integer                    12          12      0.02s  |
|   o-> comms_map_transpose                         2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                 11                         |
|  /                                                                          |
| O-> nlpot_prepare_precon_ks                      11          11      0.30s  |
|    /                                                                        |
|   o-> ion_set_projectors                         11          11      0.00s  |
|   o-> comms_copy_gv_logical                      11          11      0.00s  |
|   o-> comms_copy_bnd_logical                     11          11      0.00s  |
|   o-> ion_beta_beta_recip_cmplx                   5           5      0.28s  |
|   o-> comms_reduce_gv_complex                     5           5      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- coeffs_dot_all_many_many                  451                         |
|  /o-- coeffs_dot_lower_many_many                181                         |
| |/                                                                          |
| O-> local_dot_all_many_many                     632         632      0.30s  |
|    /                                                                        |
|   o-> algor_matmul_cmplx_cmplx                  632         632      0.28s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    11                         |
|  /                                                                          |
| O-> electronic_apply_H_energy_eigen              11          11      0.30s  |
|    /                                                                        |
|   o-> wave_kinetic_eigenvalues_wv_ks             11          11      0.01s  |
|   o-> nlpot_calc_eigenvalues_nkns                11          11      0.28s  |
|   o-> comms_reduce_bnd_real                      33          33      0.00s  |
|   o-> comms_reduce_kp_real                       33          33      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> ion_initialise                                1           1      0.30s  |
|    /                                                                        |
|   o-> ion_allocate                                1           1      0.00s  |
|   o-> ion_atom_allocate_pspot                     4           4      0.01s  |
|   o-> ion_atom_read_usp                           2           2      0.02s  |
|   o-> ion_set_data                                2           2      0.00s  |
|   o-> ion_atom_deallocate_pspot                   5           5      0.00s  |
|   o-> ion_set_psp                                 2           2      0.00s  |
|   o-> ion_atom_pseudo_scf                         2           2      0.26s  |
|   o-> ion_clebsch_gordan                          1           1      0.00s  |
|   o-> ion_atom_radial_transform                   1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_diagonalise_H_ks                      11                         |
|  /o-- hamiltonian_diagonalise_ks                181                         |
| |/                                                                          |
| O-> algor_diagonalise_complex                   192         192      0.29s  |
|    /                                                                        |
|   o-> algor_diagonalise_hermitian_lapack        192         192      0.28s  |
+-----------------------------------------------------------------------------+
|   o-- algor_diagonalise_complex                 192                         |
|  /                                                                          |
| O-> algor_diagonalise_hermitian_lapack          192         192      0.28s  |
|    /                                                                        |
|   o-> comms_copy_gv_complex                     192         192      0.01s  |
|   o-> comms_copy_gv_real                        192         192      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_apply_H_energy_local             1                         |
|  /o-- electronic_apply_H_energy_eigen            11                         |
| |/                                                                          |
| O-> nlpot_calc_eigenvalues_nkns                  12          12      0.28s  |
|    /                                                                        |
|   o-> nlpot_calc_eigenvals_nkns_r                12          12      0.28s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calc_eigenvalues_nkns                12                         |
|  /                                                                          |
| O-> nlpot_calc_eigenvals_nkns_r                  12          12      0.28s  |
|    /                                                                        |
|   o-> wave_beta_phi_wv_ks                        12          12      0.28s  |
|   o-> comms_reduce_gv_real                       12          12      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_prepare_precon_ks                     5                         |
|  /                                                                          |
| O-> ion_beta_beta_recip_cmplx                     5           5      0.28s  |
|    /                                                                        |
|   o-> ion_beta_recip_set                          5           5      0.00s  |
|   o-> comms_reduce_bnd_complex                    5           5      0.03s  |
|   o-> comms_reduce_gv_complex                     5           5      0.08s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_basis                392                         |
|  /o-- ion_atom_init_pseudo_H                     44                         |
|  /o-- ion_atom_pseudo_hartree                    66                         |
|  /o-- ion_atom_pseudo_xc                         66                         |
|  /o-- ion_atom_ps_diag                        54675                         |
|  /o-- ion_atom_calc_pseudo_rho                   62                         |
|  /o-- ion_atom_set_pseudo_H                     742                         |
|  /o-- ion_atom_pseudo_scf                     14007                         |
| |/                                                                          |
| O-> ion_atom_regin                            70054       70054      0.27s  |
+-----------------------------------------------------------------------------+
|   o-- wave_dot_all_wv_wv_ks                      11                         |
|  /o-- wave_dot_all_slice_slice                  258                         |
|  /o-- wave_coeffs_Sdot_all_wv_slice             181                         |
|  /o-- wave_coeffs_Sdot_all_wv_ks                  1                         |
| |/                                                                          |
| O-> coeffs_dot_all_many_many                    451         451      0.27s  |
|    /                                                                        |
|   o-> local_dot_all_many_many                   451         451      0.27s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthogonalise_wv_slice              181                         |
|  /                                                                          |
| O-> wave_orthog_over_wv_slice                   181         181      0.26s  |
|    /                                                                        |
|   o-> algor_matmul_cmplx_cmplx_3D               362         362      0.25s  |
+-----------------------------------------------------------------------------+
|   o-- basis_calculate_cut_off                     1                         |
|  /o-- castep_report_storage                       2                         |
|  /o-- ewald_calculate_energy                      3                         |
|  /o-- electronic_apply_H_energy_local             2                         |
|  /o-- electronic_find_fermi_free                526                         |
|  /o-- electronic_occupancy_update                33                         |
|  /o-- electronic_entropy_correction              11                         |
|  /o-- electronic_apply_H_energy_eigen            33                         |
|  /o-- density_calc_soft_wvfn_real                 9                         |
|  /o-- ewald_calculate_forces                      3                         |
|  /o-- nlpot_calculate_forces_r                  160                         |
|  /o-- ewald_calculate_stress                      3                         |
|  /o-- wave_kinetic_stress_wv                      1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
|  /o-- popn_calculate_populations                  1                         |
| |/                                                                          |
| O-> comms_reduce_kp_real                        789         789      0.25s  |
|    /                                                                        |
|   o-> comms_reduce_array_real                   181         181      0.13s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_stress_r                  160                         |
|  /                                                                          |
| O-> wave_dbetadcell_phi_wv_iks                  160         160      0.25s  |
|    /                                                                        |
|   o-> ion_dbetadcell_recip_ion                  160         160      0.06s  |
|   o-> comms_reduce_gv_complex                   160         160      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- comms_reduce_gv_real                      250                         |
|  /o-- comms_reduce_kp_real                      181                         |
|  /o-- comms_reduce_bnd_real                      40                         |
| |/                                                                          |
| O-> comms_reduce_array_real                     471         471      0.25s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_sedc_calculate_cell                    1                         |
|  /                                                                          |
| O-> dftd_sedc_update_TS                           1           1      0.22s  |
|    /                                                                        |
|   o-> hirshfeld_calculate                         1           1      0.22s  |
+-----------------------------------------------------------------------------+
|   o-- ion_set_Q_at_origin_recip                  48                         |
|  /o-- ion_set_dQdcell_at_origin_recip          1824                         |
| |/                                                                          |
| O-> ion_Q_recip_interpolation                  1872        1872      0.22s  |
+-----------------------------------------------------------------------------+
|   o-- ion_int_Q_at_origin_recip               23040                         |
|  /                                                                          |
| O-> algor_re_dot_cmplx_cmplx                  23040       23040      0.21s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_forces_r                    1                         |
|  /                                                                          |
| O-> nlpot_calculate_dddR                          1           1      0.19s  |
|    /                                                                        |
|   o-> basis_real_to_recip_grid                    1           1      0.00s  |
|   o-> basis_recip_fine_to_half_grid               1           1      0.00s  |
|   o-> ion_int_Q_at_origin_recip                4608        4608      0.15s  |
+-----------------------------------------------------------------------------+
|   o-- model_write_all                             4                         |
|  /                                                                          |
| O-> density_write                                 4           4      0.19s  |
|    /                                                                        |
|   o-> density_allocate                            4           4      0.00s  |
|   o-> density_copy                                4           4      0.00s  |
|   o-> density_write_parallel                      4           4      0.19s  |
|   o-> density_deallocate                          4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_write                               4                         |
|  /                                                                          |
| O-> density_write_parallel                        4           4      0.19s  |
|    /                                                                        |
|   o-> comms_gather_gv_integer                    12          12      0.00s  |
|   o-> comms_gather_gv_real                        4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_q_dot_all_many_many_c                182                         |
|  /o-- wave_q_dot_lower_many_many_c              181                         |
| |/                                                                          |
| O-> local_q_dot_all_many_many_c                 363         363      0.18s  |
|    /                                                                        |
|   o-> weight_beta_phi_many_cmplx                363         363      0.12s  |
|   o-> algor_matmul_cmplx_cmplx                  363         363      0.04s  |
+-----------------------------------------------------------------------------+
|   o-- wave_diagonalise_H_ks                      22                         |
|  /                                                                          |
| O-> wave_rotate_wv_ks                            22          22      0.17s  |
|    /                                                                        |
|   o-> wave_allocate_slice                        22          22      0.00s  |
|   o-> wave_initialise_slice                      22          22      0.02s  |
|   o-> algor_matmul_cmplx_cmplx_3D                33          33      0.13s  |
|   o-> wave_copy_slice_wv                         22          22      0.02s  |
|   o-> wave_deallocate_slice                      22          22      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                258                         |
|  /                                                                          |
| O-> wave_dot_all_slice_slice                    258         258      0.16s  |
|    /                                                                        |
|   o-> coeffs_dot_all_many_many                  258         258      0.08s  |
|   o-> comms_reduce_gv_complex                   258         258      0.05s  |
+-----------------------------------------------------------------------------+
|   o-- ion_set_Q_at_origin_recip                 156                         |
|  /o-- ion_set_dQdcell_at_origin_recip         31680                         |
| |/                                                                          |
| O-> ion_apply_and_add_ylm                     31836       31836      0.16s  |
+-----------------------------------------------------------------------------+
|   o-- ion_int_Q_at_origin_recip               23040                         |
|  /o-- ion_Q_apply_and_add_wsf                   351                         |
| |/                                                                          |
| O-> ion_set_Q_at_origin_recip                 23391       23391      0.16s  |
|    /                                                                        |
|   o-> ion_generate_QLnm                          48          48      0.01s  |
|   o-> ion_atom_radial_transform                  50          50      0.06s  |
|   o-> ion_Q_recip_interpolation                  48          48      0.00s  |
|   o-> ion_apply_and_add_ylm                     156         156      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_forces_r                  160                         |
|  /                                                                          |
| O-> wave_dbetadR_phi_wv                         160         160      0.15s  |
|    /                                                                        |
|   o-> ion_dbetadR_recip_ion                     160         160      0.01s  |
|   o-> comms_reduce_gv_complex                   160         160      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     9                         |
|  /                                                                          |
| O-> density_calculate_soft_wvfn                   9           9      0.15s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                          18          18      0.00s  |
|   o-> density_allocate                            9           9      0.00s  |
|   o-> density_calc_soft_wvfn_real                 9           9      0.14s  |
+-----------------------------------------------------------------------------+
|   o-- ion_set_dQdcell_at_origin_recip         31680                         |
|  /                                                                          |
| O-> ion_apply_and_add_ylmp                    31680       31680      0.15s  |
+-----------------------------------------------------------------------------+
|   o-- density_calculate_soft_wvfn                 9                         |
|  /                                                                          |
| O-> density_calc_soft_wvfn_real                   9           9      0.14s  |
|    /                                                                        |
|   o-> wave_recip_to_real_wv_bks                1169        1169      0.13s  |
|   o-> comms_reduce_kp_real                        9           9      0.00s  |
|   o-> comms_reduce_bnd_real                       9           9      0.00s  |
|   o-> basis_real_std_to_fine_gamma                9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_prepare_H                       10                         |
|  /o-- firstd_calculate_forces                     1                         |
|  /o-- firstd_calculate_stress                     1                         |
|  /o-- write_local_pot_and_den                     1                         |
| |/                                                                          |
| O-> locpot_calculate                             13          13      0.14s  |
|    /                                                                        |
|   o-> pot_complex_to_real                        13          13      0.00s  |
|   o-> pot_allocate                               13          13      0.00s  |
|   o-> pot_zero                                   26          26      0.00s  |
|   o-> locpot_check_locps_core                    13          13      0.02s  |
|   o-> hartree_calculate_potential                13          13      0.01s  |
|   o-> pot_add                                    26          26      0.00s  |
|   o-> pot_calc_energy_real                       33          33      0.00s  |
|   o-> xc_calculate_potential                     13          13      0.06s  |
|   o-> cell_cart_lattice_to_abc                   13          13      0.00s  |
|   o-> ewald_dipole_corr                          13          13      0.05s  |
|   o-> pot_deallocate                             13          13      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                181                         |
|  /                                                                          |
| O-> wave_Sorthog_lower_slice_s                  181         181      0.14s  |
|    /                                                                        |
|   o-> wave_Sdot_lower_slice                     181         181      0.10s  |
|   o-> wave_orthog_lower_over_slice_s            181         181      0.04s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    11                         |
|  /                                                                          |
| O-> electronic_find_occupancies                  11          11      0.14s  |
|    /                                                                        |
|   o-> electronic_order_eigenvalues               11          11      0.00s  |
|   o-> electronic_permute_eigenstates             22          22      0.01s  |
|   o-> electronic_occupancy_update                11          11      0.13s  |
+-----------------------------------------------------------------------------+
|   o-- wave_coeffs_Sdot_all_wv_slice             181                         |
|  /o-- wave_coeffs_Sdot_all_wv_ks                  1                         |
| |/                                                                          |
| O-> wave_q_dot_all_many_many_c                  182         182      0.14s  |
|    /                                                                        |
|   o-> local_q_dot_all_many_many_c               182         182      0.14s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     9                         |
|  /                                                                          |
| O-> density_augment                               9           9      0.14s  |
|    /                                                                        |
|   o-> density_augment_complex                     9           9      0.14s  |
+-----------------------------------------------------------------------------+
|   o-- density_augment                             9                         |
|  /                                                                          |
| O-> density_augment_complex                       9           9      0.14s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                           9           9      0.00s  |
|   o-> ion_augment_charge_nospin_kp                9           9      0.14s  |
+-----------------------------------------------------------------------------+
|   o-- density_augment_complex                     9                         |
|  /                                                                          |
| O-> ion_augment_charge_nospin_kp                  9           9      0.14s  |
|    /                                                                        |
|   o-> comms_reduce_gv_complex                     9           9      0.00s  |
|   o-> comms_reduce_bnd_complex                    9           9      0.00s  |
|   o-> comms_reduce_kp_complex                     9           9      0.00s  |
|   o-> ion_Q_apply_and_add_wsf                    18          18      0.12s  |
|   o-> basis_recip_half_to_fine_grid               9           9      0.00s  |
|   o-> basis_recip_grid_to_real                    9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_find_occupancies                11                         |
|  /                                                                          |
| O-> electronic_occupancy_update                  11          11      0.13s  |
|    /                                                                        |
|   o-> electronic_find_fermi_energy               11          11      0.13s  |
|   o-> comms_reduce_bnd_real                      33          33      0.00s  |
|   o-> comms_reduce_kp_real                       33          33      0.00s  |
|   o-> electronic_entropy_correction              11          11      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_occupancy_update                11                         |
|  /                                                                          |
| O-> electronic_find_fermi_energy                 11          11      0.13s  |
|    /                                                                        |
|   o-> electronic_find_fermi_free                 11          11      0.13s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_find_fermi_energy               11                         |
|  /                                                                          |
| O-> electronic_find_fermi_free                   11          11      0.13s  |
|    /                                                                        |
|   o-> parameters_band_degeneracy                537         537      0.00s  |
|   o-> algor_broadening                         1672        1672      0.00s  |
|   o-> comms_reduce_bnd_real                     526         526      0.00s  |
|   o-> comms_reduce_kp_real                      526         526      0.10s  |
|   o-> comms_gcopy_logical                       515         515      0.00s  |
|   o-> comms_gcopy_real                          504         504      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_calc_soft_wvfn_real              1169                         |
|  /                                                                          |
| O-> wave_recip_to_real_wv_bks                  1169        1169      0.13s  |
|    /                                                                        |
|   o-> basis_recip_reduced_real_one             1169        1169      0.13s  |
+-----------------------------------------------------------------------------+
|   o-- basis_calculate_cut_off                     1                         |
|  /o-- density_initialise_cell                     1                         |
|  /o-- castep_report_storage                       2                         |
|  /o-- ewald_calculate_energy                      3                         |
|  /o-- hartree_calculate_potential                13                         |
|  /o-- pot_calc_energy_really_real                34                         |
|  /o-- xc_gga                                     42                         |
|  /o-- nlpot_calculate_d_real                     12                         |
|  /o-- wave_kinetic_eigenvalues_wv_ks             14                         |
|  /o-- nlpot_calc_eigenvals_nkns_r                12                         |
|  /o-- wave_add_kinetic_energy_wv_ks               9                         |
|  /o-- wave_add_kinetic_energy_slice             181                         |
|  /o-- hirshfeld_calculate                      1282                         |
|  /o-- ewald_calculate_forces                      3                         |
|  /o-- ion_int_Q_at_origin_recip                4608                         |
|  /o-- nlpot_calculate_forces_r                    1                         |
|  /o-- ewald_calculate_stress                      3                         |
|  /o-- wave_kinetic_stress_wv                      1                         |
|  /o-- hartree_calculate_stress                    1                         |
|  /o-- nlpot_calculate_dddcell_real                1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
| |/                                                                          |
| O-> comms_reduce_gv_real                       6225        6225      0.13s  |
|    /                                                                        |
|   o-> comms_reduce_array_real                   250         250      0.09s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> cell_read_wrapped                             1           1      0.13s  |
|    /                                                                        |
|   o-> cell_setup_keywords                         1           1      0.00s  |
|   o-> comms_gcopy_integer                        63          63      0.00s  |
|   o-> cell_check_keywords                         1           1      0.00s  |
|   o-> cell_read_line_real                       480         480      0.01s  |
|   o-> cell_read_line_char                       160         160      0.00s  |
|   o-> cell_allocate                               1           1      0.00s  |
|   o-> cell_abc_to_cart_lattice                    1           1      0.00s  |
|   o-> cell_calculate_volume                       2           2      0.00s  |
|   o-> cell_recip_lattice                          1           1      0.00s  |
|   o-> cell_analyse_symmetry_wrapped               1           1      0.02s  |
|   o-> cell_kpoints_mp_wrapped                     1           1      0.00s  |
|   o-> cell_read_optics_data                       1           1      0.00s  |
|   o-> cell_read_bs_data                           1           1      0.00s  |
|   o-> cell_read_spectral_data                     1           1      0.00s  |
|   o-> cell_read_phonon_data                       1           1      0.00s  |
|   o-> cell_read_phonon_fine_data                  1           1      0.00s  |
|   o-> cell_read_magres_data                       1           1      0.00s  |
|   o-> cell_read_elnes_data                        1           1      0.00s  |
|   o-> cell_read_supercell_data                    1           1      0.00s  |
|   o-> cell_generate_cell_constraints              1           1      0.00s  |
|   o-> cell_check_cell_constraints                 1           1      0.00s  |
|   o-> cell_detect_primitive                       1           1      0.00s  |
|   o-> cell_check_group                            1           1      0.00s  |
|   o-> cell_count_symmetry_translations_wrap       1           1      0.00s  |
|   o-> cell_symmetry_test                          2           2      0.00s  |
|   o-> cell_generate_ionic_constraints             1           1      0.00s  |
|   o-> comms_gcopy_logical                        11          11      0.00s  |
|   o-> comms_gcopy_real                           67          67      0.00s  |
|   o-> comms_gcopy_character                     164         164      0.00s  |
|   o-> cell_generate_qpoints_local                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_recip_to_real_wv_bks                1169                         |
|  /                                                                          |
| O-> basis_recip_reduced_real_one               1169        1169      0.13s  |
|    /                                                                        |
|   o-> comms_transpose                          2338        2338      0.08s  |
+-----------------------------------------------------------------------------+
|   o-- ion_augment_charge_nospin_kp               18                         |
|  /                                                                          |
| O-> ion_Q_apply_and_add_wsf                      18          18      0.12s  |
|    /                                                                        |
|   o-> algor_matmul_cmplx_cmplx                   90          90      0.10s  |
|   o-> ion_set_Q_at_origin_recip                 351         351      0.00s  |
|   o-> comms_reduce_bnd_complex                   18          18      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- local_q_dot_all_many_many_c               363                         |
|  /                                                                          |
| O-> weight_beta_phi_many_cmplx                  363         363      0.12s  |
+-----------------------------------------------------------------------------+
|   o-- model_write_all                             2                         |
|  /                                                                          |
| O-> wave_write_all                                2           2      0.11s  |
|    /                                                                        |
|   o-> wave_write_all_par                          2           2      0.11s  |
+-----------------------------------------------------------------------------+
|   o-- wave_write_all                              2                         |
|  /                                                                          |
| O-> wave_write_all_par                            2           2      0.11s  |
|    /                                                                        |
|   o-> comms_gather_kp_integer                     2           2      0.00s  |
|   o-> comms_gather_gv_integer                     8           8      0.00s  |
|   o-> comms_gcopy_integer                         2           2      0.00s  |
|   o-> comms_gather_gv_complex                   304         304      0.00s  |
|   o-> comms_send_integer                         12          12      0.00s  |
|   o-> comms_recv_real                             6           6      0.00s  |
|   o-> comms_recv_integer                         24          24      0.00s  |
|   o-> comms_recv_complex                        912         912      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> parameters_read                               1           1      0.11s  |
|    /                                                                        |
|   o-> parameters_keywords_setup                   1           1      0.00s  |
|   o-> comms_gcopy_integer                         3           3      0.00s  |
|   o-> parameters_read_xc_block                    1           1      0.00s  |
|   o-> parameters_validate                         1           1      0.00s  |
|   o-> parameters_bcast                            1           1      0.00s  |
|   o-> algor_set_random_seed                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_populations              25440                         |
|  /                                                                          |
| O-> cell_verify_mixture_components            25440       25440      0.11s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_output_forces                      160                         |
|  /o-- popn_calculate_populations              24112                         |
| |/                                                                          |
| O-> cell_format_atom_label                    24272       24272      0.11s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> ewald_calculate_stress                        1           1      0.11s  |
|    /                                                                        |
|   o-> comms_reduce_kp_logical                     2           2      0.00s  |
|   o-> comms_reduce_bnd_logical                    2           2      0.00s  |
|   o-> comms_reduce_gv_logical                     2           2      0.00s  |
|   o-> ewald_calculate_num_cells                   1           1      0.00s  |
|   o-> comms_reduce_gv_real                        3           3      0.03s  |
|   o-> comms_reduce_bnd_real                       3           3      0.01s  |
|   o-> comms_reduce_kp_real                        3           3      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthonormalise_wv_ks                  1                         |
|  /o-- popn_calculate_matrices                     1                         |
| |/                                                                          |
| O-> wave_calc_Soverlap_wv_ks                      2           2      0.10s  |
|    /                                                                        |
|   o-> coeffs_dot_all_self                         2           2      0.01s  |
|   o-> wave_beta_phi_wv_ks                         2           2      0.08s  |
|   o-> wave_q_dot_all_self_c                       2           2      0.00s  |
|   o-> comms_reduce_gv_complex                     2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthog_lower_slice_s                181                         |
|  /                                                                          |
| O-> wave_Sdot_lower_slice                       181         181      0.10s  |
|    /                                                                        |
|   o-> coeffs_dot_lower_many_many                181         181      0.03s  |
|   o-> wave_beta_phi_slice                       362         362      0.01s  |
|   o-> wave_q_dot_lower_many_many_c              181         181      0.04s  |
|   o-> comms_reduce_gv_complex                   181         181      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_mulliken                     1                         |
|  /                                                                          |
| O-> popn_calculate_matrices                       1           1      0.10s  |
|    /                                                                        |
|   o-> wave_calc_Soverlap_wv_ks                    1           1      0.06s  |
|   o-> wave_Sdot_all_wv_wv_ks                      1           1      0.04s  |
|   o-> algor_invert_complex                        1           1      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_pseudo_scf                         4                         |
|  /                                                                          |
| O-> ion_atom_init_pseudo_atom                     4           4      0.09s  |
|    /                                                                        |
|   o-> ion_atom_resolve_pseudo_cfg                 4           4      0.00s  |
|   o-> ion_atom_locate                          1336        1336      0.01s  |
|   o-> ion_atom_interpolate                     1336        1336      0.01s  |
|   o-> ion_atom_set_pseudo_occ                     4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_initialise_cell                     1                         |
|  /o-- locps_calculate_potential                   2                         |
|  /o-- hartree_calculate_potential                13                         |
|  /o-- density_gradient_cmplx                     42                         |
|  /o-- xc_calc_vector_divergence_cmplx            14                         |
|  /o-- ion_augment_charge_nospin_kp                9                         |
|  /o-- dm_mix_density_to_density                   9                         |
|  /o-- hirshfeld_calculate                       320                         |
| |/                                                                          |
| O-> basis_recip_grid_to_real                    410         410      0.09s  |
|    /                                                                        |
|   o-> comms_transpose                           820         820      0.05s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> model_initialise                              1           1      0.09s  |
|    /                                                                        |
|   o-> model_reset                                 1           1      0.00s  |
|   o-> cell_allocate                               2           2      0.00s  |
|   o-> cell_copy                                   2           2      0.00s  |
|   o-> wave_allocate_wv                            1           1      0.00s  |
|   o-> wave_initialise_wv                          1           1      0.05s  |
|   o-> density_allocate                            1           1      0.00s  |
|   o-> density_initialise_cell                     1           1      0.03s  |
|   o-> model_init_occ_eigenvalues_ef               1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ewald_calc_storage                          1                         |
|  /o-- ewald_calculate_energy                      1                         |
|  /o-- ewald_dipole_corr                          15                         |
|  /o-- ewald_calculate_forces                      1                         |
|  /o-- ewald_calculate_stress                      1                         |
| |/                                                                          |
| O-> ewald_calculate_num_cells                    19          19      0.08s  |
|    /                                                                        |
|   o-> comms_gcopy_real                           95          95      0.01s  |
|   o-> comms_gcopy_logical                        57          57      0.00s  |
|   o-> comms_gcopy_character                      19          19      0.00s  |
|   o-> comms_gcopy_integer                       114         114      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_forces                     1                         |
|  /                                                                          |
| O-> ewald_calculate_forces                        1           1      0.08s  |
|    /                                                                        |
|   o-> comms_reduce_kp_logical                     2           2      0.00s  |
|   o-> comms_reduce_bnd_logical                    2           2      0.00s  |
|   o-> comms_reduce_gv_logical                     2           2      0.00s  |
|   o-> ewald_calculate_num_cells                   1           1      0.00s  |
|   o-> comms_reduce_gv_real                        3           3      0.02s  |
|   o-> comms_reduce_bnd_real                       3           3      0.01s  |
|   o-> comms_reduce_kp_real                        3           3      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- wave_initialise_wv                          3                         |
|  /o-- wave_beta_phi_wv_ks                       199                         |
|  /o-- ion_beta_recip_set                       1245                         |
|  /o-- wave_calc_storage_wv                        3                         |
|  /o-- wave_calc_storage_slice                     3                         |
|  /o-- wave_calc_storage_bnd                       1                         |
|  /o-- wave_setup_slice                         1061                         |
|  /o-- wave_beta_phi_slice                      1790                         |
|  /o-- nlpot_prepare_precon_ks                    11                         |
| |/                                                                          |
| O-> ion_set_projectors                         4316        4316      0.08s  |
|    /                                                                        |
|   o-> comms_reduce_gv_logical                  4316        4316      0.06s  |
|   o-> comms_reduce_gv_integer                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                181                         |
|  /                                                                          |
| O-> wave_Sorthonormalise_slice                  181         181      0.07s  |
|    /                                                                        |
|   o-> wave_calc_Soverlap_slice                  181         181      0.05s  |
|   o-> wave_orthonormalise_over_slice            181         181      0.02s  |
|   o-> comms_reduce_bnd_integer                  181         181      0.00s  |
|   o-> comms_reduce_gv_integer                   181         181      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_populations                  1                         |
|  /                                                                          |
| O-> popn_sort                                     1           1      0.07s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              1                         |
|  /o-- ion_set_Q_at_origin_recip                  50                         |
|  /o-- ion_set_dQdcell_at_origin_recip          1536                         |
| |/                                                                          |
| O-> ion_atom_radial_transform                  1587        1587      0.07s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           13                         |
|  /o-- electronic_minimisation                     1                         |
|  /o-- firstd_calculate_forces                     1                         |
| |/                                                                          |
| O-> ewald_dipole_corr                            15          15      0.07s  |
|    /                                                                        |
|   o-> ewald_calculate_num_cells                  15          15      0.07s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_mulliken                     1                         |
|  /                                                                          |
| O-> popn_calculate_charge_spilling                1           1      0.07s  |
|    /                                                                        |
|   o-> popn_nbands_occupied                        1           1      0.00s  |
|   o-> comms_reduce_kp_complex                     1           1      0.06s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           13                         |
|  /o-- xc_calculate_stress_density                 1                         |
| |/                                                                          |
| O-> xc_calculate_potential                       14          14      0.06s  |
|    /                                                                        |
|   o-> pot_allocate                               14          14      0.00s  |
|   o-> pot_zero                                   14          14      0.00s  |
|   o-> xc_gga                                     14          14      0.06s  |
|   o-> pot_add                                    14          14      0.00s  |
|   o-> pot_copy                                   14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- xc_calculate_potential                     14                         |
|  /                                                                          |
| O-> xc_gga                                       14          14      0.06s  |
|    /                                                                        |
|   o-> density_gradient_cmplx                     14          14      0.01s  |
|   o-> comms_reduce_gv_real                       42          42      0.00s  |
|   o-> xc_calc_vector_divergence_cmplx            14          14      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- wave_dbetadcell_phi_wv_iks                160                         |
|  /                                                                          |
| O-> ion_dbetadcell_recip_ion                    160         160      0.06s  |
|    /                                                                        |
|   o-> ion_cc_structure_factor                   160         160      0.00s  |
|   o-> ion_beta_recip_interpolation              512         512      0.01s  |
|   o-> ion_apply_ylmp                           3072        3072      0.01s  |
|   o-> basis_multiply_recip_reduced             3072        3072      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_set_projectors                       4316                         |
|  /o-- ewald_calculate_energy                      2                         |
|  /o-- locpot_check_locps_core                    26                         |
|  /o-- electronic_store_energy                    24                         |
|  /o-- hamiltonian_diagonalise_ks                 11                         |
|  /o-- electronic_check_occupancies               11                         |
|  /o-- ewald_calculate_forces                      2                         |
|  /o-- ewald_calculate_stress                      2                         |
|  /o-- bib_output                                  1                         |
| |/                                                                          |
| O-> comms_reduce_gv_logical                    4395        4395      0.06s  |
|    /                                                                        |
|   o-> comms_lcopy                                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> ion_read                                      1           1      0.06s  |
|    /                                                                        |
|   o-> ion_atom_inquire_usp                        2           2      0.06s  |
+-----------------------------------------------------------------------------+
|   o-- ion_augment_charge_nospin_kp                9                         |
|  /o-- popn_calculate_charge_spilling              1                         |
|  /o-- popn_calculate_populations                  1                         |
| |/                                                                          |
| O-> comms_reduce_kp_complex                      11          11      0.06s  |
|    /                                                                        |
|   o-> comms_reduce_array_complex                 10          10      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_read                                    2                         |
|  /                                                                          |
| O-> ion_atom_inquire_usp                          2           2      0.06s  |
|    /                                                                        |
|   o-> comms_gcopy_integer                        10          10      0.00s  |
|   o-> comms_gcopy_real                            2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_diagonalise_H_ks                      11                         |
|  /                                                                          |
| O-> wave_dot_all_wv_wv_ks                        11          11      0.06s  |
|    /                                                                        |
|   o-> coeffs_dot_all_many_many                   11          11      0.05s  |
|   o-> comms_reduce_gv_complex                    11          11      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- model_initialise                            1                         |
|  /o-- electronic_initialise                       1                         |
|  /o-- popn_create_lcao_basis_allkpt               1                         |
| |/                                                                          |
| O-> wave_initialise_wv                            3           3      0.06s  |
|     +- section "initialisation", branch on value of method:-                |
|       R                                           1           1      0.05s  |
|         \                                                                   |
|          o-> wave_spin_type_wv                    1           1      0.00s  |
|          o-> wave_prepare_init_wvfn               1           1      0.00s  |
|          o-> algor_set_random_seed                1           1      0.00s  |
|          o-> algor_uniform_random_array         152         152      0.00s  |
|          o-> wave_Sorthonormalise_wv              1           1      0.05s  |
|       Z                                           2           2      0.00s  |
|    /                                                                        |
|   o-> ion_set_projectors                          3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_orthonormalise_over_wv_ks              1                         |
|  /o-- ion_beta_beta_recip_cmplx                   5                         |
|  /o-- hamiltonian_diagonalise_ks                181                         |
|  /o-- ion_augment_charge_nospin_kp                9                         |
|  /o-- ion_Q_apply_and_add_wsf                    18                         |
| |/                                                                          |
| O-> comms_reduce_bnd_complex                    214         214      0.05s  |
|    /                                                                        |
|   o-> comms_reduce_array_complex                214         214      0.05s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks               3542                         |
|  /o-- wave_rotate_slice                         362                         |
|  /o-- electronic_permute_eigenstates            422                         |
| |/                                                                          |
| O-> wave_copy_slice_slice                      4326        4326      0.05s  |
+-----------------------------------------------------------------------------+
|   o-- wave_initialise_wv                          1                         |
|  /                                                                          |
| O-> wave_Sorthonormalise_wv                       1           1      0.05s  |
|    /                                                                        |
|   o-> wave_Sorthonormalise_wv_ks                  1           1      0.05s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthonormalise_wv                     1                         |
|  /                                                                          |
| O-> wave_Sorthonormalise_wv_ks                    1           1      0.05s  |
|    /                                                                        |
|   o-> wave_calc_Soverlap_wv_ks                    1           1      0.05s  |
|   o-> wave_orthonormalise_over_wv_ks              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locps_calculate_potential                   4                         |
|  /o-- locps_calculate_stress                      2                         |
| |/                                                                          |
| O-> locps_structure_factor                        6           6      0.05s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                      342                         |
|  /o-- wave_rotate_wv_ks                          22                         |
|  /o-- hamiltonian_diagonalise_ks                362                         |
|  /o-- electronic_permute_eigenstates            422                         |
| |/                                                                          |
| O-> wave_copy_slice_wv                         1148        1148      0.05s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthonormalise_slice                181                         |
|  /                                                                          |
| O-> wave_calc_Soverlap_slice                    181         181      0.05s  |
|    /                                                                        |
|   o-> coeffs_dot_all_self                       181         181      0.02s  |
|   o-> wave_beta_phi_slice                       181         181      0.00s  |
|   o-> wave_q_dot_all_self_c                     181         181      0.02s  |
|   o-> comms_reduce_gv_complex                   181         181      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_all_beta_multi_phi_recip              376                         |
|  /o-- ion_beta_add_multi_recip_all              704                         |
|  /o-- ion_beta_beta_recip_cmplx                   5                         |
|  /o-- ion_dbetadR_recip_ion                     160                         |
| |/                                                                          |
| O-> ion_beta_recip_set                         1245        1245      0.04s  |
|    /                                                                        |
|   o-> ion_set_projectors                       1245        1245      0.02s  |
|   o-> ion_cc_structure_factor                   160         160      0.00s  |
|   o-> ion_beta_recip_interpolation              512         512      0.00s  |
|   o-> basis_multiply_recip_reduced              512         512      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_count_plane_waves                     3                         |
|  /o-- ion_set_projectors                          1                         |
|  /o-- wave_prepare_init_wvfn                      2                         |
|  /o-- dm_count_plane_waves                        6                         |
|  /o-- wave_Sorthonormalise_slice                181                         |
| |/                                                                          |
| O-> comms_reduce_gv_integer                     193         193      0.04s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_count_plane_waves                       1           1      0.04s  |
|    /                                                                        |
|   o-> comms_reduce_gv_integer                     3           3      0.04s  |
|   o-> comms_reduce_kp_integer                     3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> locpot_calculate_stress                       1           1      0.04s  |
|    /                                                                        |
|   o-> hartree_calculate_stress                    1           1      0.00s  |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> pot_zero                                    1           1      0.00s  |
|   o-> locps_calculate_stress                      1           1      0.04s  |
|   o-> xc_calculate_stress_density                 1           1      0.00s  |
|   o-> pot_deallocate                              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sdot_lower_slice                     181                         |
|  /                                                                          |
| O-> wave_q_dot_lower_many_many_c                181         181      0.04s  |
|    /                                                                        |
|   o-> local_q_dot_all_many_many_c               181         181      0.04s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> electronic_initialise                         1           1      0.04s  |
|    /                                                                        |
|   o-> density_allocate                            1           1      0.00s  |
|   o-> density_zero                                1           1      0.00s  |
|   o-> wave_allocate_wv                            1           1      0.00s  |
|   o-> wave_initialise_wv                          1           1      0.00s  |
|   o-> ion_real_initialise                         1           1      0.00s  |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> hubbard_initialise                          1           1      0.00s  |
|   o-> nlpot_allocate_nl_d                         1           1      0.00s  |
|   o-> electronic_restore                          1           1      0.00s  |
|   o-> dm_flush_history                            1           1      0.00s  |
|   o-> dm_mix_density                              1           1      0.00s  |
|   o-> ewald_calculate_energy                      1           1      0.04s  |
|   o-> locps_calculate_non_coulomb                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthog_lower_slice_s                181                         |
|  /                                                                          |
| O-> wave_orthog_lower_over_slice_s              181         181      0.04s  |
|    /                                                                        |
|   o-> algor_matmul_cmplx_cmplx_3D               362         362      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /                                                                          |
| O-> ewald_calculate_energy                        1           1      0.04s  |
|    /                                                                        |
|   o-> comms_reduce_kp_logical                     2           2      0.00s  |
|   o-> comms_reduce_bnd_logical                    2           2      0.00s  |
|   o-> comms_reduce_gv_logical                     2           2      0.00s  |
|   o-> ewald_calculate_num_cells                   1           1      0.00s  |
|   o-> comms_reduce_gv_real                        3           3      0.01s  |
|   o-> comms_reduce_bnd_real                       3           3      0.00s  |
|   o-> comms_reduce_kp_real                        3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_check_locps_core                     1                         |
|  /o-- locps_calculate_stress                      1                         |
| |/                                                                          |
| O-> locps_calculate_potential                     2           2      0.04s  |
|    /                                                                        |
|   o-> pot_zero                                    2           2      0.00s  |
|   o-> basis_radial_to_recip_grid                  4           4      0.00s  |
|   o-> basis_scale_recip_grid                      4           4      0.00s  |
|   o-> locps_structure_factor                      4           4      0.03s  |
|   o-> basis_multiply_recip_grid                   4           4      0.00s  |
|   o-> basis_recip_grid_to_real                    2           2      0.00s  |
|   o-> basis_scale_real_grid                       2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_matrices                     1                         |
|  /                                                                          |
| O-> wave_Sdot_all_wv_wv_ks                        1           1      0.04s  |
|    /                                                                        |
|   o-> wave_coeffs_Sdot_all_wv_ks                  1           1      0.03s  |
|   o-> comms_reduce_gv_complex                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate_stress                     1                         |
|  /                                                                          |
| O-> locps_calculate_stress                        1           1      0.04s  |
|    /                                                                        |
|   o-> locps_calculate_non_coulomb                 1           1      0.00s  |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> locps_calculate_potential                   1           1      0.02s  |
|   o-> pot_calc_energy_real                        1           1      0.00s  |
|   o-> pot_deallocate                              1           1      0.00s  |
|   o-> density_to_recip                            1           1      0.00s  |
|   o-> basis_radial_to_recip_grid                  2           2      0.00s  |
|   o-> basis_scale_recip_grid                      2           2      0.00s  |
|   o-> locps_structure_factor                      2           2      0.02s  |
|   o-> basis_multiply_recip_grid                   2           2      0.00s  |
|   o-> basis_sum_recip_grid                       12          12      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sdot_all_wv_wv_ks                      1                         |
|  /                                                                          |
| O-> wave_coeffs_Sdot_all_wv_ks                    1           1      0.03s  |
|    /                                                                        |
|   o-> coeffs_dot_all_many_many                    1           1      0.01s  |
|   o-> wave_beta_phi_wv_ks                         2           2      0.02s  |
|   o-> wave_q_dot_all_many_many_c                  1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> write_local_pot_and_den                       1           1      0.03s  |
|    /                                                                        |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> locpot_calculate                            1           1      0.01s  |
|   o-> pot_write                                   1           1      0.02s  |
|   o-> pot_deallocate                              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_forces                     1                         |
|  /                                                                          |
| O-> locpot_calculate_forces                       1           1      0.03s  |
|    /                                                                        |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> pot_zero                                    1           1      0.00s  |
|   o-> locps_calculate_forces                      1           1      0.03s  |
|   o-> pot_deallocate                              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate_forces                     1                         |
|  /                                                                          |
| O-> locps_calculate_forces                        1           1      0.03s  |
|    /                                                                        |
|   o-> density_to_recip                            1           1      0.00s  |
|   o-> basis_radial_to_recip_grid                  2           2      0.00s  |
|   o-> basis_scale_recip_grid                      2           2      0.00s  |
|   o-> basis_multiply_recip_grid                   2           2      0.00s  |
|   o-> basis_sum_recip_grid                      480         480      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- model_initialise                            1                         |
|  /                                                                          |
| O-> density_initialise_cell                       1           1      0.03s  |
|    /                                                                        |
|   o-> density_zero                                1           1      0.00s  |
|   o-> basis_radial_to_recip_grid                  2           2      0.00s  |
|   o-> basis_recip_grid_to_real                    1           1      0.00s  |
|   o-> comms_reduce_gv_real                        1           1      0.00s  |
|   o-> density_complex_to_real                     1           1      0.00s  |
|   o-> comms_copy_bnd_real                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sdot_lower_slice                     181                         |
|  /                                                                          |
| O-> coeffs_dot_lower_many_many                  181         181      0.03s  |
|    /                                                                        |
|   o-> local_dot_all_many_many                   181         181      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- model_write_all                             8                         |
|  /                                                                          |
| O-> cell_dump                                     8           8      0.03s  |
|    /                                                                        |
|   o-> cell_dump_cell                              8           8      0.01s  |
|   o-> cell_dump_global                            8           8      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_slice                   181                         |
|  /                                                                          |
| O-> wave_add_kinetic_energy_slice               181         181      0.03s  |
|    /                                                                        |
|   o-> local_kinetic_energy_many                 181         181      0.01s  |
|   o-> comms_reduce_gv_real                      181         181      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- wave_initialise_slice                    1061                         |
|  /                                                                          |
| O-> wave_setup_slice                           1061        1061      0.03s  |
|    /                                                                        |
|   o-> ion_set_projectors                       1061        1061      0.02s  |
|   o-> wave_band_basis_initialise                346         346      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_calc_Soverlap_wv_ks                    2                         |
|  /o-- wave_calc_Soverlap_slice                  181                         |
| |/                                                                          |
| O-> coeffs_dot_all_self                         183         183      0.03s  |
|    /                                                                        |
|   o-> local_dot_all_self                        183         183      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- coeffs_dot_all_self                       183                         |
|  /                                                                          |
| O-> local_dot_all_self                          183         183      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- ewald_calculate_energy                      3                         |
|  /o-- nlpot_calculate_d_real                     12                         |
|  /o-- electronic_apply_H_energy_local             2                         |
|  /o-- hamiltonian_diagonalise_ks                 11                         |
|  /o-- electronic_find_fermi_free                526                         |
|  /o-- electronic_occupancy_update                33                         |
|  /o-- electronic_entropy_correction              11                         |
|  /o-- electronic_apply_H_energy_eigen            33                         |
|  /o-- density_calc_soft_wvfn_real                 9                         |
|  /o-- model_write_occ_eigenvalues                 8                         |
|  /o-- ewald_calculate_forces                      3                         |
|  /o-- nlpot_calculate_forces_r                    1                         |
|  /o-- ewald_calculate_stress                      3                         |
|  /o-- wave_kinetic_stress_wv                      1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
|  /o-- write_eigenvalues                           1                         |
| |/                                                                          |
| O-> comms_reduce_bnd_real                       658         658      0.03s  |
|    /                                                                        |
|   o-> comms_reduce_array_real                    40          40      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- dm_density_to_mix_density                  10                         |
|  /o-- density_to_recip                           30                         |
|  /o-- xc_calc_vector_divergence_cmplx            42                         |
|  /o-- nlpot_calculate_d_real                     12                         |
|  /o-- nlpot_calculate_dddR                        1                         |
|  /o-- nlpot_calculate_dddcell_real                1                         |
| |/                                                                          |
| O-> basis_real_to_recip_grid                     96          96      0.02s  |
|    /                                                                        |
|   o-> comms_transpose                           192         192      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- model_write_all                             4                         |
|  /                                                                          |
| O-> parameters_dump                               4           4      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- wave_calc_Soverlap_wv_ks                    2                         |
|  /o-- wave_calc_Soverlap_slice                  181                         |
| |/                                                                          |
| O-> wave_q_dot_all_self_c                       183         183      0.02s  |
|    /                                                                        |
|   o-> local_q_dot_all_self_c                    183         183      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- wave_q_dot_all_self_c                     183                         |
|  /                                                                          |
| O-> local_q_dot_all_self_c                      183         183      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthonormalise_slice                181                         |
|  /                                                                          |
| O-> wave_orthonormalise_over_slice              181         181      0.02s  |
|    /                                                                        |
|   o-> wave_beta_phi_slice                       181         181      0.00s  |
|   o-> algor_invert_complex                      181         181      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- write_local_pot_and_den                     1                         |
|  /                                                                          |
| O-> pot_write                                     1           1      0.02s  |
|    /                                                                        |
|   o-> pot_write_parallel                          1           1      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- pot_write                                   1                         |
|  /                                                                          |
| O-> pot_write_parallel                            1           1      0.02s  |
|    /                                                                        |
|   o-> comms_gather_gv_integer                     3           3      0.00s  |
|   o-> comms_gather_gv_real                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_dump                                   8                         |
|  /                                                                          |
| O-> cell_dump_global                              8           8      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                      342                         |
|  /o-- hamiltonian_diagonalise_ks                154                         |
|  /o-- electronic_permute_eigenstates            480                         |
| |/                                                                          |
| O-> wave_copy_wv_slice                          976         976      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              2                         |
|  /                                                                          |
| O-> ion_atom_read_usp                             2           2      0.02s  |
|    /                                                                        |
|   o-> ion_atom_derivative                         2           2      0.00s  |
|   o-> ion_atom_radin                             20          20      0.00s  |
|   o-> comms_gcopy_real                           36          36      0.00s  |
|   o-> comms_gcopy_integer                        14          14      0.00s  |
|   o-> comms_gcopy_logical                         4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           13                         |
|  /                                                                          |
| O-> locpot_check_locps_core                      13          13      0.02s  |
|    /                                                                        |
|   o-> comms_reduce_gv_logical                    26          26      0.00s  |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> pot_zero                                    1           1      0.00s  |
|   o-> locps_calculate_potential                   1           1      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- basis_distribute_grids                     12                         |
|  /o-- model_write_occ_eigenvalues                 4                         |
|  /o-- wave_write_all_par                          2                         |
|  /o-- write_eigenvalues                           1                         |
| |/                                                                          |
| O-> comms_gather_kp_integer                      19          19      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /o-- cell_supercell                              2                         |
| |/                                                                          |
| O-> cell_generate_ionic_constraints               3           3      0.02s  |
|    /                                                                        |
|   o-> algor_uniform_random                     2304        2304      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_beta_recip_set                        512                         |
|  /o-- ion_dbetadcell_recip_ion                  512                         |
| |/                                                                          |
| O-> ion_beta_recip_interpolation               1024        1024      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_pseudo_scf                        62                         |
|  /                                                                          |
| O-> ion_atom_set_pseudo_H                        62          62      0.02s  |
|    /                                                                        |
|   o-> ion_atom_calc_pseudo_rho                   62          62      0.00s  |
|   o-> ion_atom_pseudo_hartree                    62          62      0.00s  |
|   o-> ion_atom_pseudo_xc                         62          62      0.01s  |
|   o-> ion_atom_regin                            742         742      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_analyse_symmetry_wrapped                 1           1      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- hirshfeld_calculate                         2                         |
|  /                                                                          |
| O-> cell_supercell                                2           2      0.02s  |
|    /                                                                        |
|   o-> cell_deallocate                             2           2      0.00s  |
|   o-> cell_num_supercells                         2           2      0.00s  |
|   o-> algor_invert_real                           2           2      0.00s  |
|   o-> cell_allocate                               2           2      0.00s  |
|   o-> cell_recip_lattice                          2           2      0.00s  |
|   o-> cell_calculate_volume                       2           2      0.00s  |
|   o-> cell_generate_supercell_origins             2           2      0.00s  |
|   o-> cell_copy_kpoints                           2           2      0.00s  |
|   o-> cell_supercell_reduce_kpoints_wrapped       2           2      0.00s  |
|   o-> cell_set_supercell_symmetry                 2           2      0.00s  |
|   o-> cell_generate_ionic_constraints             2           2      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_sedc_calculate_cell                    3                         |
|  /                                                                          |
| O-> dftd_sedc_update_cell                         3           3      0.02s  |
|    /                                                                        |
|   o-> cell_frac_to_cart_vector_wrapped          480         480      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /o-- cell_read_wrapped                          67                         |
|  /o-- cell_generate_qpoints_local                 4                         |
|  /o-- ion_atom_inquire_usp                        2                         |
|  /o-- parameters_bcast                           94                         |
|  /o-- ion_atom_read_usp                          36                         |
|  /o-- ewald_calculate_num_cells                  95                         |
|  /o-- electronic_find_fermi_free                504                         |
|  /o-- hirshfeld_calculate                         4                         |
|  /o-- dftd_sedc_calculate_energy_cell             1                         |
|  /o-- dftd_sedc_calculate_forces_cell             1                         |
|  /o-- dftd_sedc_calculate_stress_cell             1                         |
| |/                                                                          |
| O-> comms_gcopy_real                            810         810      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- xc_gga                                     14                         |
|  /                                                                          |
| O-> density_gradient_cmplx                       14          14      0.01s  |
|    /                                                                        |
|   o-> density_to_recip                           14          14      0.00s  |
|   o-> basis_recip_grid_to_real                   42          42      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- xc_gga                                     14                         |
|  /                                                                          |
| O-> xc_calc_vector_divergence_cmplx              14          14      0.01s  |
|    /                                                                        |
|   o-> basis_real_to_recip_grid                   42          42      0.01s  |
|   o-> basis_recip_grid_to_real                   14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_orthonormalise_over_wv_ks              1                         |
|  /o-- wave_orthonormalise_over_slice            181                         |
|  /o-- popn_calculate_matrices                     1                         |
| |/                                                                          |
| O-> algor_invert_complex                        183         183      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_mulliken                     1                         |
|  /                                                                          |
| O-> popn_calculate_density_matrix                 1           1      0.01s  |
|    /                                                                        |
|   o-> popn_nbands_occupied                        1           1      0.00s  |
|   o-> popn_invert_complex                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_add_kinetic_energy_wv_ks               9                         |
|  /o-- wave_add_kinetic_energy_slice             181                         |
| |/                                                                          |
| O-> local_kinetic_energy_many                   190         190      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> parameters_output                             1           1      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- model_write                                 1                         |
|  /                                                                          |
| O-> io_delete_file                                1           1      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_pseudo_scf                         4                         |
|  /                                                                          |
| O-> ion_atom_init_pseudo_basis                    4           4      0.01s  |
|    /                                                                        |
|   o-> ion_atom_find_root                        792         792      0.00s  |
|   o-> ion_atom_regin                            392         392      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_dump                                   8                         |
|  /                                                                          |
| O-> cell_dump_cell                                8           8      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /o-- electronic_minimisation                     9                         |
| |/                                                                          |
| O-> dm_mix_density                               10          10      0.01s  |
|    /                                                                        |
|   o-> dm_mix_density_pulay                       10          10      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density                             10                         |
|  /                                                                          |
| O-> dm_mix_density_pulay                         10          10      0.01s  |
|    /                                                                        |
|   o-> dm_initialise                               1           1      0.00s  |
|   o-> dm_density_to_mix_density                   9           9      0.00s  |
|   o-> dm_mix_density_kerker                       1           1      0.00s  |
|   o-> dm_mix_density_copy                        40          40      0.00s  |
|   o-> dm_mix_density_add                        104         104      0.00s  |
|   o-> dm_mix_density_dot                        164         164      0.00s  |
|   o-> dm_apply_kerker                             8           8      0.00s  |
|   o-> dm_mix_density_to_density                   8           8      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_forces                     1                         |
|  /                                                                          |
| O-> dftd_sedc_calc_forces_curr_cell               1           1      0.01s  |
|    /                                                                        |
|   o-> dftd_sedc_calculate_forces_cell             1           1      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_sedc_calc_forces_curr_cell             1                         |
|  /                                                                          |
| O-> dftd_sedc_calculate_forces_cell               1           1      0.01s  |
|    /                                                                        |
|   o-> dftd_sedc_calculate_cell                    1           1      0.01s  |
|   o-> comms_gcopy_real                            1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                         480                         |
|  /                                                                          |
| O-> cell_read_line_real                         480         480      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- wave_dbetadR_phi_wv                       160                         |
|  /                                                                          |
| O-> ion_dbetadR_recip_ion                       160         160      0.01s  |
|    /                                                                        |
|   o-> ion_beta_recip_set                        160         160      0.01s  |
|   o-> ion_proj_index                            512         512      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_write_all_par                        912                         |
|  /                                                                          |
| O-> comms_recv_complex                          912         912      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- check_forces_stresses                       1                         |
|  /                                                                          |
| O-> firstd_output_forces                          1           1      0.01s  |
|    /                                                                        |
|   o-> cell_format_atom_label                    160         160      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locps_calculate_forces                    480                         |
|  /o-- locps_calculate_stress                     12                         |
| |/                                                                          |
| O-> basis_sum_recip_grid                        492         492      0.01s  |
|    /                                                                        |
|   o-> comms_reduce_gv_complex                   492         492      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_apply_H_energy_local             1                         |
|  /o-- electronic_apply_H_energy_eigen            11                         |
|  /o-- hamiltonian_diagonalise_ks                  2                         |
| |/                                                                          |
| O-> wave_kinetic_eigenvalues_wv_ks               14          14      0.01s  |
|    /                                                                        |
|   o-> comms_reduce_gv_real                       14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_generate_ionic_constraints          2304                         |
|  /                                                                          |
| O-> algor_uniform_random                       2304        2304      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- hartree_calculate_potential                13                         |
|  /o-- density_gradient_cmplx                     14                         |
|  /o-- locps_calculate_forces                      1                         |
|  /o-- hartree_calculate_stress                    1                         |
|  /o-- locps_calculate_stress                      1                         |
| |/                                                                          |
| O-> density_to_recip                             30          30      0.01s  |
|    /                                                                        |
|   o-> basis_real_to_recip_grid                   30          30      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- wave_orthonormalise_over_wv_ks              1                         |
|  /o-- algor_diagonalise_hermitian_lapack        192                         |
| |/                                                                          |
| O-> comms_copy_gv_complex                       193         193      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_beta_recip_set                        512                         |
|  /o-- ion_dbetadcell_recip_ion                 3072                         |
| |/                                                                          |
| O-> basis_multiply_recip_reduced               3584        3584      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                        9                         |
|  /                                                                          |
| O-> wave_add_kinetic_energy_wv_ks                 9           9      0.01s  |
|    /                                                                        |
|   o-> local_kinetic_energy_many                   9           9      0.01s  |
|   o-> comms_reduce_gv_real                        9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_H                      4                         |
|  /o-- ion_atom_set_pseudo_H                      62                         |
| |/                                                                          |
| O-> ion_atom_pseudo_xc                           66          66      0.01s  |
|    /                                                                        |
|   o-> ion_atom_derivative                       132         132      0.00s  |
|   o-> ion_atom_regin                             66          66      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           13                         |
|  /                                                                          |
| O-> hartree_calculate_potential                  13          13      0.01s  |
|    /                                                                        |
|   o-> hartree_check_inv_gsqr                     13          13      0.00s  |
|   o-> density_to_recip                           13          13      0.00s  |
|   o-> comms_reduce_gv_real                       13          13      0.00s  |
|   o-> basis_recip_grid_to_real                   13          13      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              4                         |
|  /o-- ion_generate_orbitals                       2                         |
| |/                                                                          |
| O-> ion_atom_allocate_pspot                       6           6      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_dbetadcell_recip_ion                 3072                         |
|  /                                                                          |
| O-> ion_apply_ylmp                             3072        3072      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- cell_kpoints_mp_wrapped                     1                         |
|  /o-- cell_generate_MP_set_wrapped                1                         |
|  /o-- castep                                      5                         |
|  /o-- dftd_sedc_initialize                        1                         |
| |/                                                                          |
| O-> bib_add                                       8           8      0.01s  |
|    /                                                                        |
|   o-> bib_setup                                   1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> cell_output_wrapped                           1           1      0.01s  |
|    /                                                                        |
|   o-> cell_cart_lattice_to_abc                    1           1      0.00s  |
|   o-> cell_check_group                            1           1      0.00s  |
|   o-> cell_factor_group_symmetry_wrapped          1           1      0.00s  |
|   o-> cell_symmetry_symbol_wrapped_wrapped        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_apply_slice                             9                         |
|  /                                                                          |
| O-> pot_interpolate                               9           9      0.01s  |
|    /                                                                        |
|   o-> pot_debug_check                             9           9      0.00s  |
|   o-> basis_real_fine_to_std_grid                 9           9      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_find_occupancies                22                         |
|  /                                                                          |
| O-> electronic_permute_eigenstates               22          22      0.01s  |
|    /                                                                        |
|   o-> wave_allocate_slice                        44          44      0.00s  |
|   o-> wave_initialise_slice                      44          44      0.00s  |
|   o-> wave_copy_wv_slice                        480         480      0.00s  |
|   o-> wave_copy_slice_wv                        422         422      0.00s  |
|   o-> wave_copy_slice_slice                     422         422      0.00s  |
|   o-> wave_deallocate_slice                      44          44      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_set_Q_at_origin_recip                  48                         |
|  /                                                                          |
| O-> ion_generate_QLnm                            48          48      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> dftd_sedc_calc_stress_curr_cell               1           1      0.01s  |
|    /                                                                        |
|   o-> dftd_sedc_calculate_stress_cell             1           1      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_sedc_calc_stress_curr_cell             1                         |
|  /                                                                          |
| O-> dftd_sedc_calculate_stress_cell               1           1      0.01s  |
|    /                                                                        |
|   o-> dftd_sedc_calculate_cell                    1           1      0.01s  |
|   o-> comms_gcopy_real                            1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_interpolate                             9                         |
|  /                                                                          |
| O-> basis_real_fine_to_std_grid                   9           9      0.01s  |
|    /                                                                        |
|   o-> comms_transpose                            36          36      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              2                         |
|  /o-- ion_generate_orbitals                       2                         |
| |/                                                                          |
| O-> ion_set_psp                                   4           4      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- check_elec_ground_state                     1                         |
|  /                                                                          |
| O-> castep_calc_storage                           1           1      0.01s  |
|    /                                                                        |
|   o-> cell_calc_storage_global_wrapped            1           1      0.00s  |
|   o-> basis_calc_storage                          1           1      0.00s  |
|   o-> ion_calc_storage                            1           1      0.00s  |
|   o-> model_calc_storage                          1           1      0.00s  |
|   o-> ewald_calc_storage                          1           1      0.00s  |
|   o-> density_calc_storage                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_atom                1336                         |
|  /                                                                          |
| O-> ion_atom_interpolate                       1336        1336      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_atom                1336                         |
|  /                                                                          |
| O-> ion_atom_locate                            1336        1336      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                          11                         |
|  /o-- cell_generate_qpoints_local                 2                         |
|  /o-- parameters_bcast                           63                         |
|  /o-- ion_atom_read_usp                           4                         |
|  /o-- ewald_calculate_num_cells                  57                         |
|  /o-- electronic_minimisation                     1                         |
|  /o-- hubbard_initialise                          3                         |
|  /o-- electronic_find_fermi_free                515                         |
|  /o-- dftd_sedc_calculate_cell                    4                         |
|  /o-- hirshfeld_calculate                         4                         |
| |/                                                                          |
| O-> comms_gcopy_logical                         664         664      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              2                         |
|  /                                                                          |
| O-> ion_set_data                                  2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_calc_soft_wvfn_real                 9                         |
|  /                                                                          |
| O-> basis_real_std_to_fine_gamma                  9           9      0.00s  |
|    /                                                                        |
|   o-> basis_real_std_to_fine_grid                 9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> electronic_apply_H_energy_local               1           1      0.00s  |
|    /                                                                        |
|   o-> wave_kinetic_eigenvalues_wv_ks              1           1      0.00s  |
|   o-> nlpot_calc_eigenvalues_nkns                 1           1      0.00s  |
|   o-> comms_reduce_bnd_real                       2           2      0.00s  |
|   o-> comms_reduce_kp_real                        2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep_calc_storage                         1                         |
|  /                                                                          |
| O-> ewald_calc_storage                            1           1      0.00s  |
|    /                                                                        |
|   o-> ewald_calculate_num_cells                   1           1      0.00s  |
|   o-> algor_sizeof_real                          12          12      0.00s  |
|   o-> algor_sizeof_cmplx                          9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> comms_parallel_strategy                       1           1      0.00s  |
|    /                                                                        |
|   o-> find_strategy                               1           1      0.00s  |
|   o-> assign_nodes                                1           1      0.00s  |
|   o-> reassign_nodes                              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate_stress                     1                         |
|  /                                                                          |
| O-> xc_calculate_stress_density                   1           1      0.00s  |
|    /                                                                        |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> xc_calculate_potential                      1           1      0.00s  |
|   o-> pot_deallocate                              1           1      0.00s  |
|   o-> xc_calculate_stress_correction              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_real_std_to_fine_gamma                9                         |
|  /                                                                          |
| O-> basis_real_std_to_fine_grid                   9           9      0.00s  |
|    /                                                                        |
|   o-> comms_transpose                            36          36      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_beta_phi_wv_ks                        14                         |
|  /o-- wave_beta_phi_slice                      1790                         |
| |/                                                                          |
| O-> wave_calc_ps_q_nonzero                     1804        1804      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthonormalise_wv_ks                  1                         |
|  /                                                                          |
| O-> wave_orthonormalise_over_wv_ks                1           1      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_bnd_complex                    1           1      0.00s  |
|   o-> algor_invert_complex                        1           1      0.00s  |
|   o-> comms_copy_gv_complex                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- parameters_read                             1                         |
|  /                                                                          |
| O-> parameters_bcast                              1           1      0.00s  |
|    /                                                                        |
|   o-> comms_gcopy_logical                        63          63      0.00s  |
|   o-> comms_gcopy_character                      93          93      0.00s  |
|   o-> comms_gcopy_integer                        77          77      0.00s  |
|   o-> parameters_reallocate_xc                    7           7      0.00s  |
|   o-> comms_gcopy_real                           94          94      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_create_lcao_basis_allkpt             256                         |
|  /                                                                          |
| O-> basis_radial_to_recip_reduced               256         256      0.00s  |
|    /                                                                        |
|   o-> basis_utils_interpolation                 256         256      0.00s  |
|   o-> basis_utils_apply_s_harmonic              160         160      0.00s  |
|   o-> basis_utils_apply_p_harmonic               96          96      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_basis                792                         |
|  /                                                                          |
| O-> ion_atom_find_root                          792         792      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_write_all_par                        304                         |
|  /                                                                          |
| O-> comms_gather_gv_complex                     304         304      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_distribute_grids                      2                         |
|  /                                                                          |
| O-> comms_map_transpose                           2           2      0.00s  |
|    /                                                                        |
|   o-> comms_map_transpose_n                       2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- comms_map_transpose                         2                         |
|  /                                                                          |
| O-> comms_map_transpose_n                         2           2      0.00s  |
|    /                                                                        |
|   o-> comms_local_map                             2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                         160                         |
|  /                                                                          |
| O-> cell_read_line_char                         160         160      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- comms_map_transpose_n                       2                         |
|  /                                                                          |
| O-> comms_local_map                               2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                697                         |
|  /o-- wave_Sorthonormalise_slice                181                         |
| |/                                                                          |
| O-> comms_reduce_bnd_integer                    878         878      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_radial_to_recip_grid                 14                         |
|  /o-- basis_radial_to_recip_reduced             256                         |
| |/                                                                          |
| O-> basis_utils_interpolation                   270         270      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> hirshfeld_output                              1           1      0.00s  |
|    /                                                                        |
|   o-> parameters_nspins                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_pulay                        9                         |
|  /o-- dm_mix_density_kerker                       1                         |
| |/                                                                          |
| O-> dm_density_to_mix_density                    10          10      0.00s  |
|    /                                                                        |
|   o-> basis_real_to_recip_grid                   10          10      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_beta_recip_set                        160                         |
|  /o-- ion_dbetadcell_recip_ion                  160                         |
| |/                                                                          |
| O-> ion_cc_structure_factor                     320         320      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_apply_slice                           523                         |
|  /o-- pot_nongamma_apply_slice                  523                         |
|  /o-- nlpot_apply_add_slice_r                   523                         |
| |/                                                                          |
| O-> wave_spin_type_slice                       1569        1569      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_find_fermi_free                537                         |
|  /                                                                          |
| O-> parameters_band_degeneracy                  537         537      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> dftd_sedc_initialize                          1           1      0.00s  |
|    /                                                                        |
|   o-> bib_add                                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_distribute_grids                     12                         |
|  /o-- density_write_parallel                     12                         |
|  /o-- wave_write_all_par                          8                         |
|  /o-- pot_write_parallel                          3                         |
| |/                                                                          |
| O-> comms_gather_gv_integer                      35          35      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_kerker                       1                         |
|  /o-- dm_mix_density_pulay                      164                         |
| |/                                                                          |
| O-> dm_mix_density_dot                          165         165      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_gv_complex                   165         165      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_find_fermi_free               1672                         |
|  /                                                                          |
| O-> algor_broadening                           1672        1672      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                        9                         |
|  /                                                                          |
| O-> wave_dot_wv_wv_ks                             9           9      0.00s  |
|    /                                                                        |
|   o-> local_dot_many                              9           9      0.00s  |
|   o-> comms_reduce_gv_complex                     9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                       18                         |
|  /o-- wave_rotate_wv_ks                          22                         |
|  /o-- wave_rotate_slice                         362                         |
|  /o-- hamiltonian_diagonalise_ks                 55                         |
|  /o-- electronic_permute_eigenstates             44                         |
| |/                                                                          |
| O-> wave_deallocate_slice                       501         501      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_kerker                       1                         |
|  /o-- dm_mix_density_pulay                        8                         |
| |/                                                                          |
| O-> dm_mix_density_to_density                     9           9      0.00s  |
|    /                                                                        |
|   o-> basis_recip_grid_to_real                    9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                 55                         |
|  /o-- hamiltonian_apply_ks                       18                         |
|  /o-- wave_rotate_wv_ks                          22                         |
|  /o-- wave_rotate_slice                         362                         |
|  /o-- electronic_permute_eigenstates             44                         |
| |/                                                                          |
| O-> wave_allocate_slice                         501         501      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                         164                         |
|  /o-- parameters_bcast                           93                         |
|  /o-- ewald_calculate_num_cells                  19                         |
|  /o-- electronic_minimisation                     1                         |
| |/                                                                          |
| O-> comms_gcopy_character                       277         277      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_magres_data                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_prepare_precon_ks                    11                         |
|  /o-- hamiltonian_diagonalise_ks                543                         |
| |/                                                                          |
| O-> comms_copy_gv_logical                       554         554      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    12                         |
|  /                                                                          |
| O-> electronic_write_scf_energies                12          12      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_ps_diag                          245                         |
|  /                                                                          |
| O-> ion_atom_rectoreal                          245         245      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                          63                         |
|  /o-- cell_generate_qpoints_local                 4                         |
|  /o-- ion_atom_inquire_usp                       10                         |
|  /o-- parameters_read                             3                         |
|  /o-- parameters_bcast                           77                         |
|  /o-- parameters_reallocate_xc                    7                         |
|  /o-- ion_atom_read_usp                          14                         |
|  /o-- ewald_calculate_num_cells                 114                         |
|  /o-- hirshfeld_calculate                         2                         |
|  /o-- dftd_sedc_calculate_cell                    1                         |
|  /o-- wave_write_all_par                          2                         |
| |/                                                                          |
| O-> comms_gcopy_integer                         297         297      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           33                         |
|  /o-- locps_calculate_stress                      1                         |
| |/                                                                          |
| O-> pot_calc_energy_real                         34          34      0.00s  |
|    /                                                                        |
|   o-> pot_debug_check                            34          34      0.00s  |
|   o-> pot_calc_energy_really_real                34          34      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> bib_output                                    1           1      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_farm_logical                   1           1      0.00s  |
|   o-> comms_reduce_bnd_logical                    1           1      0.00s  |
|   o-> comms_reduce_gv_logical                     1           1      0.00s  |
|   o-> comms_reduce_kp_logical                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_dot_wv_wv_ks                           9                         |
|  /                                                                          |
| O-> local_dot_many                                9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_zero                                   45                         |
|  /o-- pot_add                                    80                         |
|  /o-- pot_calc_energy_real                       34                         |
|  /o-- pot_calc_energy_really_real                34                         |
|  /o-- pot_copy                                   28                         |
|  /o-- pot_apply_slice                           523                         |
|  /o-- pot_interpolate                             9                         |
| |/                                                                          |
| O-> pot_debug_check                             753         753      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_elnes_data                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_density_matrix               1                         |
|  /                                                                          |
| O-> popn_invert_complex                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_initialise_cell                     2                         |
|  /o-- locps_calculate_potential                   4                         |
|  /o-- hirshfeld_calculate                         4                         |
|  /o-- locps_calculate_forces                      2                         |
|  /o-- locps_calculate_stress                      2                         |
| |/                                                                          |
| O-> basis_radial_to_recip_grid                   14          14      0.00s  |
|    /                                                                        |
|   o-> basis_utils_interpolation                  14          14      0.00s  |
|   o-> basis_utils_apply_s_harmonic               14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_calc_energy_real                       34                         |
|  /                                                                          |
| O-> pot_calc_energy_really_real                  34          34      0.00s  |
|    /                                                                        |
|   o-> pot_debug_check                            34          34      0.00s  |
|   o-> comms_reduce_gv_real                       34          34      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_prepare_H                       10                         |
|  /o-- electronic_minimisation                     9                         |
| |/                                                                          |
| O-> density_symmetrise                           19          19      0.00s  |
|    /                                                                        |
|   o-> density_real_to_complex                    19          19      0.00s  |
|   o-> density_symmetrise_gv_parallel             19          19      0.00s  |
|   o-> density_complex_to_real                    19          19      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_set_pseudo_H                      62                         |
|  /                                                                          |
| O-> ion_atom_calc_pseudo_rho                     62          62      0.00s  |
|    /                                                                        |
|   o-> ion_atom_regin                             62          62      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_write_all                             4                         |
|  /                                                                          |
| O-> model_write_occ_eigenvalues                   4           4      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_bnd_real                       8           8      0.00s  |
|   o-> comms_gather_kp_integer                     4           4      0.00s  |
|   o-> comms_recv_real                            36          36      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- check_elec_ground_state                     1                         |
|  /                                                                          |
| O-> electronic_calc_storage                       1           1      0.00s  |
|    /                                                                        |
|   o-> wave_calc_storage_wv                        2           2      0.00s  |
|   o-> wave_calc_storage_slice                     1           1      0.00s  |
|   o-> dm_calc_storage                             1           1      0.00s  |
|   o-> density_calc_storage                        1           1      0.00s  |
|   o-> pot_calc_storage                            1           1      0.00s  |
|   o-> algor_sizeof_real                           2           2      0.00s  |
|   o-> algor_sizeof_cmplx                          1           1      0.00s  |
|   o-> hamiltonian_calc_storage                    1           1      0.00s  |
|   o-> wave_calc_storage_bnd                       1           1      0.00s  |
|   o-> locpot_calc_storage                         1           1      0.00s  |
|   o-> density_symmetrise_calc_storage             1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_phonon_fine_data                    1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_pulay                        1                         |
|  /                                                                          |
| O-> dm_initialise                                 1           1      0.00s  |
|    /                                                                        |
|   o-> dm_mix_density_initialise                   1           1      0.00s  |
|   o-> dm_mix_density_allocate                    46          46      0.00s  |
|   o-> dm_mix_density_zero                        46          46      0.00s  |
|   o-> density_allocate                            1           1      0.00s  |
|   o-> density_zero                                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_pseudo_scf                         4                         |
|  /                                                                          |
| O-> ion_atom_init_pseudo_H                        4           4      0.00s  |
|    /                                                                        |
|   o-> ion_atom_regin                             44          44      0.00s  |
|   o-> ion_atom_pseudo_hartree                     4           4      0.00s  |
|   o-> ion_atom_pseudo_xc                          4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_H                      4                         |
|  /o-- ion_atom_set_pseudo_H                      62                         |
| |/                                                                          |
| O-> ion_atom_pseudo_hartree                      66          66      0.00s  |
|    /                                                                        |
|   o-> ion_atom_regin                             66          66      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_read_usp                           2                         |
|  /o-- ion_atom_pseudo_xc                        132                         |
| |/                                                                          |
| O-> ion_atom_derivative                         134         134      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           26                         |
|  /o-- xc_calculate_potential                     14                         |
| |/                                                                          |
| O-> pot_add                                      40          40      0.00s  |
|    /                                                                        |
|   o-> pot_debug_check                            80          80      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_generate_cell_constraints              1                         |
|  /o-- wave_initialise_wv                        152                         |
| |/                                                                          |
| O-> algor_uniform_random_array                  153         153      0.00s  |
|    /                                                                        |
|   o-> algor_set_random_seed                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- algor_diagonalise_hermitian_lapack        192                         |
|  /o-- electronic_order_eigenvalues               11                         |
| |/                                                                          |
| O-> comms_copy_gv_real                          203         203      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_dbetadR_recip_ion                     512                         |
|  /                                                                          |
| O-> ion_proj_index                              512         512      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           26                         |
|  /o-- locpot_check_locps_core                     1                         |
|  /o-- locps_calculate_potential                   2                         |
|  /o-- xc_calculate_potential                     14                         |
|  /o-- locpot_calculate_forces                     1                         |
|  /o-- locpot_calculate_stress                     1                         |
| |/                                                                          |
| O-> pot_zero                                     45          45      0.00s  |
|    /                                                                        |
|   o-> pot_debug_check                            45          45      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_initialise                            1                         |
|  /o-- nlxc_initialise                             1                         |
|  /o-- electronic_initialise                       1                         |
|  /o-- dm_initialise                               1                         |
|  /o-- density_real_to_complex                    19                         |
|  /o-- density_complex_to_real                    19                         |
|  /o-- density_calculate_soft_wvfn                 9                         |
|  /o-- density_write                               4                         |
| |/                                                                          |
| O-> density_allocate                             55          55      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /o-- model_initialise                            2                         |
|  /o-- hirshfeld_calculate                         8                         |
|  /o-- cell_supercell                              2                         |
| |/                                                                          |
| O-> cell_allocate                                13          13      0.00s  |
|    /                                                                        |
|   o-> cell_deallocate                            13          13      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_symmetrise                         19                         |
|  /                                                                          |
| O-> density_real_to_complex                      19          19      0.00s  |
|    /                                                                        |
|   o-> density_allocate                           19          19      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_allocate_wv                            3                         |
|  /o-- wave_setup_slice                          346                         |
| |/                                                                          |
| O-> wave_band_basis_initialise                  349         349      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /o-- model_initialise                            2                         |
|  /o-- hirshfeld_calculate                         8                         |
| |/                                                                          |
| O-> cell_copy                                    11          11      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /o-- nlxc_initialise                             1                         |
| |/                                                                          |
| O-> cell_generate_qpoints_local                   2           2      0.00s  |
|    /                                                                        |
|   o-> cell_detect_MP                              2           2      0.00s  |
|   o-> comms_gcopy_real                            4           4      0.00s  |
|   o-> comms_gcopy_integer                         4           4      0.00s  |
|   o-> comms_gcopy_logical                         2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_kpoints_mp_wrapped                       1           1      0.00s  |
|    /                                                                        |
|   o-> bib_add                                     1           1      0.00s  |
|   o-> cell_generate_MP_set_wrapped                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> nlxc_initialise                               1           1      0.00s  |
|    /                                                                        |
|   o-> density_allocate                            1           1      0.00s  |
|   o-> density_copy                                1           1      0.00s  |
|   o-> cell_generate_qpoints_local                 1           1      0.00s  |
|   o-> density_deallocate                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    12                         |
|  /                                                                          |
| O-> electronic_store_energy                      12          12      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_gv_logical                    24          24      0.00s  |
|   o-> comms_reduce_bnd_logical                   24          24      0.00s  |
|   o-> comms_reduce_kp_logical                    24          24      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- check_forces_stresses                       1                         |
|  /                                                                          |
| O-> firstd_output_stress                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> wave_kinetic_stress_wv                        1           1      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_gv_real                        1           1      0.00s  |
|   o-> comms_reduce_bnd_real                       1           1      0.00s  |
|   o-> comms_reduce_kp_real                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              1                         |
|  /                                                                          |
| O-> ion_clebsch_gordan                            1           1      0.00s  |
|    /                                                                        |
|   o-> init_factorial                              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_phonon_data                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              5                         |
|  /o-- ion_generate_orbitals                       2                         |
| |/                                                                          |
| O-> ion_atom_deallocate_pspot                     7           7      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /o-- locpot_calculate                           13                         |
|  /o-- locpot_check_locps_core                     1                         |
|  /o-- xc_calculate_potential                     14                         |
|  /o-- locpot_calculate_forces                     1                         |
|  /o-- firstd_calculate_forces                     1                         |
|  /o-- locpot_calculate_stress                     1                         |
|  /o-- locps_calculate_stress                      1                         |
|  /o-- xc_calculate_stress_density                 1                         |
|  /o-- firstd_calculate_stress                     1                         |
|  /o-- write_local_pot_and_den                     1                         |
| |/                                                                          |
| O-> pot_allocate                                 36          36      0.00s  |
|    /                                                                        |
|   o-> pot_deallocate_real                         3           3      0.00s  |
|   o-> pot_deallocate_nc                          36          36      0.00s  |
|   o-> pot_deallocate_cmplx                       33          33      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_pulay                        1                         |
|  /                                                                          |
| O-> dm_mix_density_kerker                         1           1      0.00s  |
|    /                                                                        |
|   o-> dm_density_to_mix_density                   1           1      0.00s  |
|   o-> dm_mix_density_copy                         1           1      0.00s  |
|   o-> dm_mix_density_add                          2           2      0.00s  |
|   o-> dm_mix_density_dot                          1           1      0.00s  |
|   o-> dm_apply_kerker                             1           1      0.00s  |
|   o-> dm_mix_density_to_density                   1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- parameters_read                             1                         |
|  /                                                                          |
| O-> parameters_keywords_setup                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_calc_storage                            3                         |
|  /o-- wave_calc_storage_wv                        6                         |
|  /o-- density_calc_storage                        6                         |
|  /o-- ewald_calc_storage                          9                         |
|  /o-- phonon_store_mem_estimate                   1                         |
|  /o-- wave_calc_storage_slice                     6                         |
|  /o-- dm_calc_storage                             1                         |
|  /o-- pot_calc_storage                           14                         |
|  /o-- electronic_calc_storage                     1                         |
|  /o-- hamiltonian_calc_storage                    4                         |
|  /o-- wave_calc_storage_bnd                       2                         |
|  /o-- density_symmetrise_calc_storage             1                         |
|  /o-- locps_calc_storage                          3                         |
| |/                                                                          |
| O-> algor_sizeof_cmplx                           57          57      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_real                          57          57      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    11                         |
|  /                                                                          |
| O-> electronic_check_occupancies                 11          11      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_bnd_logical                   11          11      0.00s  |
|   o-> comms_reduce_gv_logical                    11          11      0.00s  |
|   o-> comms_reduce_kp_logical                    11          11      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ewald_calculate_energy                      2                         |
|  /o-- electronic_store_energy                    24                         |
|  /o-- hamiltonian_diagonalise_ks                192                         |
|  /o-- electronic_check_occupancies               11                         |
|  /o-- ewald_calculate_forces                      2                         |
|  /o-- ewald_calculate_stress                      2                         |
|  /o-- bib_output                                  1                         |
| |/                                                                          |
| O-> comms_reduce_bnd_logical                    234         234      0.00s  |
|    /                                                                        |
|   o-> comms_lcopy                                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate_stress                     1                         |
|  /                                                                          |
| O-> hartree_calculate_stress                      1           1      0.00s  |
|    /                                                                        |
|   o-> hartree_check_inv_gsqr                      1           1      0.00s  |
|   o-> density_to_recip                            1           1      0.00s  |
|   o-> comms_reduce_gv_real                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_prepare_init_wvfn                      1                         |
|  /o-- ewald_calculate_energy                      2                         |
|  /o-- electronic_store_energy                    24                         |
|  /o-- electronic_check_occupancies               11                         |
|  /o-- ewald_calculate_forces                      2                         |
|  /o-- ewald_calculate_stress                      2                         |
|  /o-- bib_output                                  1                         |
| |/                                                                          |
| O-> comms_reduce_kp_logical                      43          43      0.00s  |
|    /                                                                        |
|   o-> comms_lcopy                                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_kerker                       2                         |
|  /o-- dm_mix_density_pulay                      104                         |
| |/                                                                          |
| O-> dm_mix_density_add                          106         106      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_d_real                     12                         |
|  /o-- nlpot_calculate_dddR                        1                         |
|  /o-- nlpot_calculate_dddcell_real                1                         |
| |/                                                                          |
| O-> basis_recip_fine_to_half_grid                14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_find_occupancies                11                         |
|  /                                                                          |
| O-> electronic_order_eigenvalues                 11          11      0.00s  |
|    /                                                                        |
|   o-> algor_sort                                 11          11      0.00s  |
|   o-> comms_copy_gv_integer                      11          11      0.00s  |
|   o-> comms_copy_gv_real                         11          11      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           13                         |
|  /o-- electronic_finalise                         1                         |
|  /o-- locpot_calculate_forces                     1                         |
|  /o-- firstd_calculate_forces                     1                         |
|  /o-- locps_calculate_stress                      1                         |
|  /o-- xc_calculate_stress_density                 1                         |
|  /o-- locpot_calculate_stress                     1                         |
|  /o-- firstd_calculate_stress                     1                         |
|  /o-- write_local_pot_and_den                     1                         |
| |/                                                                          |
| O-> pot_deallocate                               21          21      0.00s  |
|    /                                                                        |
|   o-> pot_deallocate_real                        21          21      0.00s  |
|   o-> pot_deallocate_cmplx                       21          21      0.00s  |
|   o-> pot_deallocate_nc                          21          21      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_initialise_cell                     1                         |
|  /o-- density_symmetrise                         19                         |
| |/                                                                          |
| O-> density_complex_to_real                      20          20      0.00s  |
|    /                                                                        |
|   o-> density_allocate                           19          19      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_initialise                            1                         |
|  /o-- electronic_initialise                       1                         |
|  /o-- popn_create_lcao_basis_allkpt               1                         |
| |/                                                                          |
| O-> wave_allocate_wv                              3           3      0.00s  |
|    /                                                                        |
|   o-> wave_band_basis_initialise                  3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_radial_to_recip_grid                 14                         |
|  /o-- basis_radial_to_recip_reduced             160                         |
| |/                                                                          |
| O-> basis_utils_apply_s_harmonic                174         174      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- check_elec_ground_state                     1                         |
|  /                                                                          |
| O-> castep_report_storage                         1           1      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_gv_real                        2           2      0.00s  |
|   o-> comms_reduce_kp_real                        2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_assign_grid_coordinates                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_supercell                              2                         |
|  /                                                                          |
| O-> cell_supercell_reduce_kpoints_wrapped         2           2      0.00s  |
|    /                                                                        |
|   o-> algor_invert_real                           2           2      0.00s  |
|   o-> cell_reduce_kpoints_arg_wrapped             2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_occupancy_update                11                         |
|  /                                                                          |
| O-> electronic_entropy_correction                11          11      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_bnd_real                      11          11      0.00s  |
|   o-> comms_reduce_kp_real                       11          11      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_sedc_update_cell                     480                         |
|  /                                                                          |
| O-> cell_frac_to_cart_vector_wrapped            480         480      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_generate_MP_set_wrapped                1                         |
|  /o-- cell_reduce_kpoints_arg_wrapped             2                         |
| |/                                                                          |
| O-> cell_reduce_kpoints_internal                  3           3      0.00s  |
|    /                                                                        |
|   o-> cell_sort_kpoints_with_recip                6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_write_parallel                      4                         |
|  /o-- pot_write_parallel                          1                         |
| |/                                                                          |
| O-> comms_gather_gv_real                          5           5      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_calc_storage                     1                         |
|  /o-- firstd_calc_storage                         2                         |
| |/                                                                          |
| O-> locpot_calc_storage                           3           3      0.00s  |
|    /                                                                        |
|   o-> pot_calc_storage                            3           3      0.00s  |
|   o-> density_calc_storage                        3           3      0.00s  |
|   o-> hartree_calc_storage                        3           3      0.00s  |
|   o-> locps_calc_storage                          3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_supercell_reduce_kpoints_wrapped       2                         |
|  /                                                                          |
| O-> cell_reduce_kpoints_arg_wrapped               2           2      0.00s  |
|    /                                                                        |
|   o-> cell_reduce_kpoints_internal                2           2      0.00s  |
|   o-> cell_set_kpoints_array_wrapped              2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> electronic_finalise                           1           1      0.00s  |
|    /                                                                        |
|   o-> density_deallocate                          1           1      0.00s  |
|   o-> pot_deallocate                              1           1      0.00s  |
|   o-> wave_deallocate_wv                          4           4      0.00s  |
|   o-> hubbard_finalise                            1           1      0.00s  |
|   o-> nlpot_deallocate_nl_d                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- check_elec_ground_state                     2                         |
|  /                                                                          |
| O-> firstd_calc_storage                           2           2      0.00s  |
|    /                                                                        |
|   o-> pot_calc_storage                            2           2      0.00s  |
|   o-> locpot_calc_storage                         2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_calc_storage_cell                     20                         |
|  /o-- cell_calc_storage_global_wrapped           34                         |
|  /o-- basis_calc_storage                          9                         |
|  /o-- algor_sizeof_cmplx                         57                         |
|  /o-- ion_calc_storage                           24                         |
|  /o-- model_calc_storage                          2                         |
|  /o-- ewald_calc_storage                         12                         |
|  /o-- phonon_store_mem_estimate                   2                         |
|  /o-- dm_calc_storage                             3                         |
|  /o-- electronic_calc_storage                     2                         |
|  /o-- hamiltonian_calc_storage                    6                         |
|  /o-- hartree_calc_storage                        3                         |
| |/                                                                          |
| O-> algor_sizeof_real                           174         174      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> dftd_sedc_print_corr_energies                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_kpoints_mp_wrapped                     1                         |
|  /                                                                          |
| O-> cell_generate_MP_set_wrapped                  1           1      0.00s  |
|    /                                                                        |
|   o-> bib_add                                     1           1      0.00s  |
|   o-> cell_reduce_kpoints_internal                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_reduce_kpoints_internal                6                         |
|  /                                                                          |
| O-> cell_sort_kpoints_with_recip                  6           6      0.00s  |
|    /                                                                        |
|   o-> cell_find_reduced_cell                      6           6      0.00s  |
|   o-> algor_invert_real                          12          12      0.00s  |
|   o-> algor_sort                                  6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_augment_charge_nospin_kp                9                         |
|  /                                                                          |
| O-> basis_recip_half_to_fine_grid                 9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_calc_storage                     1                         |
|  /o-- locpot_calc_storage                         3                         |
|  /o-- firstd_calc_storage                         2                         |
|  /o-- locps_calc_storage                          1                         |
| |/                                                                          |
| O-> pot_calc_storage                              7           7      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_cmplx                         14          14      0.00s  |
|   o-> algor_sizeof_logical                        7           7      0.00s  |
|   o-> algor_sizeof_int                            7           7      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> electronic_scf_footer                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_calc_storage                          1                         |
|  /o-- electronic_calc_storage                     2                         |
| |/                                                                          |
| O-> wave_calc_storage_wv                          3           3      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_cmplx                          6           6      0.00s  |
|   o-> algor_sizeof_int                            6           6      0.00s  |
|   o-> ion_set_projectors                          3           3      0.00s  |
|   o-> algor_sizeof_logical                        6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_radial_to_recip_reduced              96                         |
|  /                                                                          |
| O-> basis_utils_apply_p_harmonic                 96          96      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_reset                                 2                         |
|  /o-- electronic_finalise                         4                         |
|  /o-- popn_calculate_mulliken                     1                         |
| |/                                                                          |
| O-> wave_deallocate_wv                            7           7      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_write_occ_eigenvalues                36                         |
|  /o-- wave_write_all_par                          6                         |
|  /o-- write_eigenvalues                           9                         |
| |/                                                                          |
| O-> comms_recv_real                              51          51      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- parameters_read                             1                         |
|  /                                                                          |
| O-> parameters_read_xc_block                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_calc_storage_cell                      2                         |
|  /o-- basis_calc_storage                          3                         |
|  /o-- ion_calc_storage                            2                         |
|  /o-- wave_calc_storage_wv                        6                         |
|  /o-- density_calc_storage                        6                         |
|  /o-- wave_calc_storage_slice                     3                         |
|  /o-- pot_calc_storage                            7                         |
|  /o-- hamiltonian_calc_storage                    1                         |
|  /o-- wave_calc_storage_bnd                       1                         |
| |/                                                                          |
| O-> algor_sizeof_logical                         31          31      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_int                           31          31      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_initialise                              46                         |
|  /                                                                          |
| O-> dm_mix_density_zero                          46          46      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_calc_storage                     1                         |
|  /                                                                          |
| O-> hamiltonian_calc_storage                      1           1      0.00s  |
|    /                                                                        |
|   o-> wave_calc_storage_slice                     2           2      0.00s  |
|   o-> algor_sizeof_cmplx                          4           4      0.00s  |
|   o-> algor_sizeof_real                           6           6      0.00s  |
|   o-> algor_sizeof_logical                        1           1      0.00s  |
|   o-> algor_sizeof_int                            1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /o-- nlpot_calculate_forces_r                    1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
| |/                                                                          |
| O-> nlpot_allocate_nl_d                           3           3      0.00s  |
|    /                                                                        |
|   o-> nlpot_deallocate_nl_d                       3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_write_all_par                         24                         |
|  /                                                                          |
| O-> comms_recv_integer                           24          24      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_initialise                              46                         |
|  /                                                                          |
| O-> dm_mix_density_allocate                      46          46      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_setup_keywords                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_calc_storage                     1                         |
|  /o-- hamiltonian_calc_storage                    2                         |
| |/                                                                          |
| O-> wave_calc_storage_slice                       3           3      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_cmplx                          6           6      0.00s  |
|   o-> ion_set_projectors                          3           3      0.00s  |
|   o-> algor_sizeof_int                            3           3      0.00s  |
|   o-> algor_sizeof_logical                        3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_detect_primitive                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_calc_storage_cell                     26                         |
|  /o-- algor_sizeof_logical                       31                         |
|  /o-- cell_calc_storage_global_wrapped            7                         |
|  /o-- basis_calc_storage                         25                         |
|  /o-- ion_calc_storage                           14                         |
|  /o-- wave_calc_storage_wv                        6                         |
|  /o-- density_calc_storage                        6                         |
|  /o-- wave_calc_storage_slice                     3                         |
|  /o-- dm_calc_storage                             2                         |
|  /o-- pot_calc_storage                            7                         |
|  /o-- hamiltonian_calc_storage                    1                         |
| |/                                                                          |
| O-> algor_sizeof_int                            128         128      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- xc_calculate_potential                     14                         |
|  /                                                                          |
| O-> pot_copy                                     14          14      0.00s  |
|    /                                                                        |
|   o-> pot_debug_check                            28          28      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_calc_storage                          1                         |
|  /o-- castep_calc_storage                         1                         |
|  /o-- electronic_calc_storage                     1                         |
|  /o-- locpot_calc_storage                         3                         |
| |/                                                                          |
| O-> density_calc_storage                          6           6      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_cmplx                          6           6      0.00s  |
|   o-> algor_sizeof_int                            6           6      0.00s  |
|   o-> algor_sizeof_logical                        6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> electronic_scf_banner                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_allocate                              13                         |
|  /o-- model_reset                                 8                         |
|  /o-- cell_supercell                              2                         |
|  /o-- hirshfeld_calculate                         6                         |
| |/                                                                          |
| O-> cell_deallocate                              29          29      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_calc_storage                     1                         |
|  /                                                                          |
| O-> dm_calc_storage                               1           1      0.00s  |
|    /                                                                        |
|   o-> dm_count_plane_waves                        1           1      0.00s  |
|   o-> algor_sizeof_int                            2           2      0.00s  |
|   o-> algor_sizeof_real                           3           3      0.00s  |
|   o-> algor_sizeof_cmplx                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_count_plane_waves                     3                         |
|  /o-- electronic_analyse_occ_all                  2                         |
| |/                                                                          |
| O-> comms_reduce_kp_integer                       5           5      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_supercell                              2                         |
|  /                                                                          |
| O-> cell_set_supercell_symmetry                   2           2      0.00s  |
|    /                                                                        |
|   o-> cell_num_supercells                         2           2      0.00s  |
|   o-> algor_invert_real                           2           2      0.00s  |
|   o-> cell_reduce_symmetry_supercell              2           2      0.00s  |
|   o-> cell_find_related_atoms_supercell           2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_complex_to_real                         3                         |
|  /o-- pot_allocate                               33                         |
|  /o-- pot_deallocate                             21                         |
| |/                                                                          |
| O-> pot_deallocate_cmplx                         57          57      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_calc_storage                             1                         |
|  /o-- dm_mix_density_initialise                   1                         |
| |/                                                                          |
| O-> dm_count_plane_waves                          2           2      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_gv_integer                     6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_kerker                       1                         |
|  /o-- dm_mix_density_pulay                       40                         |
| |/                                                                          |
| O-> dm_mix_density_copy                          41          41      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_initialise                               1                         |
|  /                                                                          |
| O-> dm_mix_density_initialise                     1           1      0.00s  |
|    /                                                                        |
|   o-> dm_count_plane_waves                        1           1      0.00s  |
|   o-> dm_assign_plane_wave_indices                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    12                         |
|  /                                                                          |
| O-> electronic_write_energies                    12          12      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep_calc_storage                         1                         |
|  /                                                                          |
| O-> cell_calc_storage_global_wrapped              1           1      0.00s  |
|    /                                                                        |
|   o-> cell_calc_storage_cell                      1           1      0.00s  |
|   o-> algor_sizeof_real                          34          34      0.00s  |
|   o-> algor_sizeof_int                            7           7      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- bib_add                                     1                         |
|  /                                                                          |
| O-> bib_setup                                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- parameters_read                             1                         |
|  /                                                                          |
| O-> parameters_validate                           1           1      0.00s  |
|    /                                                                        |
|   o-> parameters_xc_allowed                       7           7      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_initialise                            1                         |
|  /o-- model_deallocate                            1                         |
| |/                                                                          |
| O-> model_reset                                   2           2      0.00s  |
|    /                                                                        |
|   o-> cell_deallocate                             8           8      0.00s  |
|   o-> wave_deallocate_wv                          2           2      0.00s  |
|   o-> density_deallocate                          2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_allocate                               36                         |
|  /o-- pot_deallocate                             21                         |
| |/                                                                          |
| O-> pot_deallocate_nc                            57          57      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep_calc_storage                         1                         |
|  /                                                                          |
| O-> model_calc_storage                            1           1      0.00s  |
|    /                                                                        |
|   o-> cell_calc_storage_cell                      1           1      0.00s  |
|   o-> wave_calc_storage_wv                        1           1      0.00s  |
|   o-> density_calc_storage                        1           1      0.00s  |
|   o-> algor_sizeof_real                           2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           13                         |
|  /                                                                          |
| O-> pot_complex_to_real                          13          13      0.00s  |
|    /                                                                        |
|   o-> pot_deallocate_cmplx                        3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_sort_kpoints_with_recip                6                         |
|  /o-- electronic_order_eigenvalues               11                         |
| |/                                                                          |
| O-> algor_sort                                   17          17      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                 11                         |
|  /                                                                          |
| O-> wave_calc_precon                             11          11      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locps_calculate_potential                   4                         |
|  /o-- locps_calculate_forces                      2                         |
|  /o-- locps_calculate_stress                      2                         |
| |/                                                                          |
| O-> basis_multiply_recip_grid                     8           8      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_calc_storage_global_wrapped            1                         |
|  /o-- model_calc_storage                          1                         |
| |/                                                                          |
| O-> cell_calc_storage_cell                        2           2      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_int                           26          26      0.00s  |
|   o-> algor_sizeof_real                          20          20      0.00s  |
|   o-> algor_sizeof_logical                        2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- parameters_bcast                            7                         |
|  /                                                                          |
| O-> parameters_reallocate_xc                      7           7      0.00s  |
|    /                                                                        |
|   o-> comms_gcopy_integer                         7           7      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep_calc_storage                         1                         |
|  /                                                                          |
| O-> ion_calc_storage                              1           1      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_int                           14          14      0.00s  |
|   o-> algor_sizeof_cmplx                          3           3      0.00s  |
|   o-> algor_sizeof_real                          24          24      0.00s  |
|   o-> algor_sizeof_logical                        2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_forces_r                    1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
| |/                                                                          |
| O-> wave_beta_phi_wv                              2           2      0.00s  |
|    /                                                                        |
|   o-> wave_beta_phi_wv_ks                         2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_atom                   4                         |
|  /                                                                          |
| O-> ion_atom_set_pseudo_occ                       4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locps_calculate_potential                   4                         |
|  /o-- locps_calculate_forces                      2                         |
|  /o-- locps_calculate_stress                      2                         |
| |/                                                                          |
| O-> basis_scale_recip_grid                        8           8      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_allocate                                3                         |
|  /o-- pot_deallocate                             21                         |
| |/                                                                          |
| O-> pot_deallocate_real                          24          24      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /                                                                          |
| O-> dm_flush_history                              1           1      0.00s  |
|    /                                                                        |
|   o-> dm_finalise                                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_find_reduced_cell                      6                         |
|  /o-- cell_sort_kpoints_with_recip               12                         |
|  /o-- cell_check_group                            2                         |
|  /o-- cell_symmetry_symbol_wrapped_wrapped        1                         |
|  /o-- cell_supercell                              2                         |
|  /o-- cell_generate_supercell_origins             2                         |
|  /o-- cell_supercell_reduce_kpoints_wrapped       2                         |
|  /o-- cell_set_supercell_symmetry                 2                         |
| |/                                                                          |
| O-> algor_invert_real                            29          29      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_initialise_wv                          1                         |
|  /o-- model_init_occ_eigenvalues_ef               1                         |
|  /o-- density_calculate_soft_wvfn                18                         |
|  /o-- density_augment_complex                     9                         |
|  /o-- nlpot_calculate_forces                      1                         |
|  /o-- nlpot_calculate_forces_r                    1                         |
|  /o-- nlpot_calculate_stress                      1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
| |/                                                                          |
| O-> wave_spin_type_wv                            33          33      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_order_eigenvalues               11                         |
|  /                                                                          |
| O-> comms_copy_gv_integer                        11          11      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- check_elec_ground_state                     3                         |
|  /                                                                          |
| O-> phonon_store_mem_estimate                     3           3      0.00s  |
|    /                                                                        |
|   o-> secondd_store_mem_estimate                  3           3      0.00s  |
|   o-> algor_sizeof_cmplx                          1           1      0.00s  |
|   o-> algor_sizeof_real                           2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_initialise_wv                          1                         |
|  /                                                                          |
| O-> wave_prepare_init_wvfn                        1           1      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_gv_integer                     2           2      0.00s  |
|   o-> comms_reduce_kp_logical                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_flush_history                            1                         |
|  /                                                                          |
| O-> dm_finalise                                   1           1      0.00s  |
|    /                                                                        |
|   o-> density_deallocate                          1           1      0.00s  |
|   o-> dm_mix_density_deallocate                   6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep_calc_storage                         1                         |
|  /                                                                          |
| O-> basis_calc_storage                            1           1      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_int                           25          25      0.00s  |
|   o-> algor_sizeof_real                           9           9      0.00s  |
|   o-> algor_sizeof_logical                        3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_sort_kpoints_with_recip                6                         |
|  /                                                                          |
| O-> cell_find_reduced_cell                        6           6      0.00s  |
|    /                                                                        |
|   o-> algor_invert_real                           6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_calc_storage                     1                         |
|  /                                                                          |
| O-> wave_calc_storage_bnd                         1           1      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_cmplx                          2           2      0.00s  |
|   o-> ion_set_projectors                          1           1      0.00s  |
|   o-> algor_sizeof_logical                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> model_deallocate                              1           1      0.00s  |
|    /                                                                        |
|   o-> model_reset                                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              1                         |
|  /                                                                          |
| O-> ion_allocate                                  1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_assign_plane_wave_indexes               1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_symmetrise                         19                         |
|  /                                                                          |
| O-> density_symmetrise_gv_parallel               19          19      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /o-- hubbard_calculate_forces                    1                         |
|  /o-- hubbard_calculate_stress                    1                         |
| |/                                                                          |
| O-> hubbard_initialise                            3           3      0.00s  |
|    /                                                                        |
|   o-> hubbard_is_on                               3           3      0.00s  |
|   o-> comms_gcopy_logical                         3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calc_storage                         3                         |
|  /                                                                          |
| O-> locps_calc_storage                            3           3      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_cmplx                          3           3      0.00s  |
|   o-> pot_calc_storage                            1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hartree_calculate_potential                13                         |
|  /o-- hartree_calculate_stress                    1                         |
| |/                                                                          |
| O-> hartree_check_inv_gsqr                       14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    11                         |
|  /                                                                          |
| O-> electronic_write_spin_density                11          11      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_initialise                   1                         |
|  /                                                                          |
| O-> dm_assign_plane_wave_indices                  1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_set_supercell_symmetry                 2                         |
|  /                                                                          |
| O-> cell_find_related_atoms_supercell             2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     9                         |
|  /                                                                          |
| O-> hubbard_mix_occ_matrix                        9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     9                         |
|  /                                                                          |
| O-> hubbard_calculate_occ_matrix                  9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlxc_initialise                             1                         |
|  /o-- density_write                               4                         |
| |/                                                                          |
| O-> density_copy                                  5           5      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_generate_qpoints_local                 2                         |
|  /o-- cell_set_kpoints_array_wrapped              2                         |
| |/                                                                          |
| O-> cell_detect_MP                                4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_supercell                              2                         |
|  /                                                                          |
| O-> cell_generate_supercell_origins               2           2      0.00s  |
|    /                                                                        |
|   o-> cell_num_supercells                         2           2      0.00s  |
|   o-> algor_invert_real                           2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_assign_pw_gvectors                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_read_usp                          20                         |
|  /                                                                          |
| O-> ion_atom_radin                               20          20      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_kerker                       1                         |
|  /o-- dm_mix_density_pulay                        8                         |
| |/                                                                          |
| O-> dm_apply_kerker                               9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_map_standard_to_fine                    1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_reduce_kpoints_arg_wrapped             2                         |
|  /                                                                          |
| O-> cell_set_kpoints_array_wrapped                2           2      0.00s  |
|    /                                                                        |
|   o-> cell_detect_MP                              2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_pseudo_scf                         4                         |
|  /                                                                          |
| O-> ion_atom_basis_pseudo_dealloc                 4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> electronic_write_occupancies                  1           1      0.00s  |
|    /                                                                        |
|   o-> electronic_analyse_occ_all                  1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> hubbard_calculate_stress                      1           1      0.00s  |
|    /                                                                        |
|   o-> hubbard_initialise                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_reset                                 2                         |
|  /o-- nlxc_initialise                             1                         |
|  /o-- dm_finalise                                 1                         |
|  /o-- electronic_finalise                         1                         |
|  /o-- density_write                               4                         |
| |/                                                                          |
| O-> density_deallocate                            9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_generate_cell_constraints                1           1      0.00s  |
|    /                                                                        |
|   o-> algor_uniform_random_array                  1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- algor_uniform_random_array                  1                         |
|  /o-- parameters_read                             1                         |
|  /o-- wave_initialise_wv                          1                         |
| |/                                                                          |
| O-> algor_set_random_seed                         3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_check_cell_constraints                 1                         |
|  /o-- cell_output_wrapped                         1                         |
|  /o-- locpot_calculate                           13                         |
| |/                                                                          |
| O-> cell_cart_lattice_to_abc                     15          15      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_calculate_cut_off                       1           1      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_gv_real                        1           1      0.00s  |
|   o-> comms_reduce_kp_real                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_write_all_par                         12                         |
|  /                                                                          |
| O-> comms_send_integer                           12          12      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calc_storage                         3                         |
|  /                                                                          |
| O-> hartree_calc_storage                          3           3      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_real                           3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_finalise                                 6                         |
|  /                                                                          |
| O-> dm_mix_density_deallocate                     6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_map_fine_recip_half_full                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_write_occupancies                1                         |
|  /                                                                          |
| O-> electronic_analyse_occ_all                    1           1      0.00s  |
|    /                                                                        |
|   o-> comms_copy_bnd_integer                      2           2      0.00s  |
|   o-> comms_reduce_kp_integer                     2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    11                         |
|  /                                                                          |
| O-> electronic_dump                              11          11      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_supercell                              2                         |
|  /o-- cell_generate_supercell_origins             2                         |
|  /o-- cell_set_supercell_symmetry                 2                         |
| |/                                                                          |
| O-> cell_num_supercells                           6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_forces                     1                         |
|  /                                                                          |
| O-> hubbard_calculate_forces                      1           1      0.00s  |
|    /                                                                        |
|   o-> hubbard_initialise                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_supercell                              2                         |
|  /                                                                          |
| O-> cell_copy_kpoints                             2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           2                         |
|  /                                                                          |
| O-> cell_symmetry_test                            2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /o-- cell_output_wrapped                         1                         |
| |/                                                                          |
| O-> cell_check_group                              2           2      0.00s  |
|    /                                                                        |
|   o-> algor_invert_real                           2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_set_supercell_symmetry                 2                         |
|  /                                                                          |
| O-> cell_reduce_symmetry_supercell                2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locps_calculate_potential                   2                         |
|  /                                                                          |
| O-> basis_scale_real_grid                         2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_initialise_cell                     1                         |
|  /o-- electronic_initialise                       1                         |
|  /o-- dm_initialise                               1                         |
| |/                                                                          |
| O-> density_zero                                  3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hirshfeld_calculate                         2                         |
|  /                                                                          |
| O-> hirshfeld_init                                2           2      0.00s  |
|    /                                                                        |
|   o-> parameters_nspins                           4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_calc_storage                     1                         |
|  /                                                                          |
| O-> density_symmetrise_calc_storage               1           1      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_cmplx                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_prepare_precon_ks                    11                         |
|  /                                                                          |
| O-> comms_copy_bnd_logical                       11          11      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_abc_to_cart_lattice                      1           1      0.00s  |
|    /                                                                        |
|   o-> cell_bravais_lattice                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_spectral_data                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_allocate_nl_d                         3                         |
|  /o-- electronic_finalise                         1                         |
|  /o-- nlpot_calculate_forces_r                    1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
| |/                                                                          |
| O-> nlpot_deallocate_nl_d                         6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hubbard_initialise                          3                         |
|  /                                                                          |
| O-> hubbard_is_on                                 3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_bs_data                             1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- phonon_store_mem_estimate                   3                         |
|  /                                                                          |
| O-> secondd_store_mem_estimate                    3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           2                         |
|  /o-- cell_supercell                              2                         |
| |/                                                                          |
| O-> cell_calculate_volume                         4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- write_eigenvalues                           4                         |
|  /                                                                          |
| O-> global_kpoint_index                           4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_initialise                            1                         |
|  /                                                                          |
| O-> model_init_occ_eigenvalues_ef                 1           1      0.00s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- parameters_validate                         7                         |
|  /                                                                          |
| O-> parameters_xc_allowed                         7           7      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- bib_output                                  1                         |
|  /                                                                          |
| O-> comms_reduce_farm_logical                     1           1      0.00s  |
|    /                                                                        |
|   o-> comms_lcopy                                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- comms_parallel_strategy                     1                         |
|  /                                                                          |
| O-> find_strategy                                 1           1      0.00s  |
|    /                                                                        |
|   o-> best_mixed_strategy                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /o-- hirshfeld_calculate                         2                         |
| |/                                                                          |
| O-> cell_distribute_kpoints_wrapped               3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            6                         |
|  /                                                                          |
| O-> basis_utils_prime_factors                     6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_count_symmetry_translations_wrappe       1           1      0.00s  |
|    /                                                                        |
|   o-> cell_factor_group_symmetry_wrapped          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_distribute_grids                      1                         |
|  /                                                                          |
| O-> basis_utils_sort_columns                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_check_cell_constraints                   1           1      0.00s  |
|    /                                                                        |
|   o-> cell_cart_lattice_to_abc                    1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> tddft_calculate_stress                        1           1      0.00s  |
|    /                                                                        |
|   o-> tddft_set_tddft_on                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_atom                   4                         |
|  /                                                                          |
| O-> ion_atom_resolve_pseudo_cfg                   4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_write                                 2                         |
|  /                                                                          |
| O-> comms_barrier                                 2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_output_wrapped                         1                         |
|  /                                                                          |
| O-> cell_symmetry_symbol_wrapped_wrapped          1           1      0.00s  |
|    /                                                                        |
|   o-> algor_invert_real                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /o-- cell_supercell                              2                         |
| |/                                                                          |
| O-> cell_recip_lattice                            3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- check_forces_stresses                       1                         |
|  /                                                                          |
| O-> firstd_symmetrise_forces                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /o-- locps_calculate_stress                      1                         |
| |/                                                                          |
| O-> locps_calculate_non_coulomb                   2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_count_symmetry_translations_wrap       1                         |
|  /o-- cell_output_wrapped                         1                         |
| |/                                                                          |
| O-> cell_factor_group_symmetry_wrapped            2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hirshfeld_init                              4                         |
|  /o-- hirshfeld_output                            1                         |
| |/                                                                          |
| O-> parameters_nspins                             5           5      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- find_strategy                               1                         |
|  /                                                                          |
| O-> best_mixed_strategy                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_analyse_occ_all                  2                         |
|  /                                                                          |
| O-> comms_copy_bnd_integer                        2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_write                                 2                         |
|  /                                                                          |
| O-> comms_barrier_farm                            2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- comms_reduce_farm_logical                   1                         |
|  /o-- comms_reduce_bnd_logical                    1                         |
|  /o-- comms_reduce_gv_logical                     1                         |
|  /o-- comms_reduce_kp_logical                     1                         |
| |/                                                                          |
| O-> comms_lcopy                                   4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_check_keywords                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /o-- electronic_initialise                       1                         |
| |/                                                                          |
| O-> ion_real_initialise                           2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_supercell_data                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_optics_data                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /o-- tddft_calculate_stress                      1                         |
| |/                                                                          |
| O-> tddft_set_tddft_on                            2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_charge_spilling              1                         |
|  /o-- popn_calculate_density_matrix               1                         |
| |/                                                                          |
| O-> popn_nbands_occupied                          2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_abc_to_cart_lattice                    1                         |
|  /                                                                          |
| O-> cell_bravais_lattice                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- comms_parallel_strategy                     1                         |
|  /                                                                          |
| O-> reassign_nodes                                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_clebsch_gordan                          1                         |
|  /                                                                          |
| O-> init_factorial                                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- xc_calculate_stress_density                 1                         |
|  /                                                                          |
| O-> xc_calculate_stress_correction                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_initialise_cell                     1                         |
|  /                                                                          |
| O-> comms_copy_bnd_real                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- comms_parallel_strategy                     1                         |
|  /                                                                          |
| O-> assign_nodes                                  1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /                                                                          |
| O-> electronic_restore                            1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_create_lcao_basis_allkpt               1                         |
|  /                                                                          |
| O-> popn_check_lcao_basis                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_finalise                         1                         |
|  /                                                                          |
| O-> hubbard_finalise                              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> memory_system_initialise                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- check_forces_stresses                       1                         |
|  /                                                                          |
| O-> firstd_symmetrise_stress                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
Class of operation                  Time spent
COMMS                                   2.19s
COMMS_GV                                1.73s
COMMS_KP                                0.33s
COMMS_BND                               0.08s
COMMS_FARM                              0.00s
     481 different subroutines and functions were traced
Hash collisions:         30
  
