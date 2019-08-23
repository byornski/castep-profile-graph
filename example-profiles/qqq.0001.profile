+-----------------------------------------------------------------------------+
|     Subroutine                              Total      Profiled   Time      |
|                                             calls      calls      (incl.)   |
+-----------------------------------------------------------------------------+
|   o-- <parent(s) not traced>                      1                         |
|  /                                                                          |
| O-> castep                                        1           1   1286.07s  |
|    /                                                                        |
|   o-> memory_system_initialise                    1           1      0.00s  |
|   o-> comms_gcopy_real                            1           1      0.00s  |
|   o-> cell_read_wrapped                           1           1      0.06s  |
|   o-> ion_read                                    1           1      0.07s  |
|   o-> parameters_read                             1           1      0.02s  |
|   o-> bib_add                                     5           5      0.00s  |
|   o-> tddft_set_tddft_on                          1           1      0.00s  |
|   o-> comms_parallel_strategy                     1           1      0.01s  |
|   o-> cell_distribute_kpoints_wrapped             1           1      0.00s  |
|   o-> ion_initialise                              1           1      0.59s  |
|   o-> basis_initialise                            1           1      2.06s  |
|   o-> ion_real_initialise                         1           1      0.00s  |
|   o-> model_initialise                            1           1      2.58s  |
|   o-> nlxc_initialise                             1           1      0.01s  |
|   o-> parameters_output                           1           1      0.01s  |
|   o-> cell_output_wrapped                         1           1      0.00s  |
|   o-> dftd_sedc_initialize                        1           1      0.00s  |
|   o-> check_elec_ground_state                     1           1    708.90s  |
|   o-> check_forces_stresses                       1           1    542.36s  |
|   o-> popn_calculate_mulliken                     1           1      3.60s  |
|   o-> write_eigenvalues                           1           1      0.11s  |
|   o-> write_local_pot_and_den                     1           1      0.98s  |
|   o-> hirshfeld_calculate                         1           1     18.14s  |
|   o-> hirshfeld_output                            1           1      0.00s  |
|   o-> model_write                                 1           1      6.55s  |
|   o-> bib_output                                  1           1      0.00s  |
|   o-> model_deallocate                            1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> check_elec_ground_state                       1           1    708.90s  |
|    /                                                                        |
|   o-> castep_calc_storage                         1           1      0.01s  |
|   o-> phonon_store_mem_estimate                   3           3      0.00s  |
|   o-> electronic_calc_storage                     1           1      0.00s  |
|   o-> firstd_calc_storage                         2           2      0.00s  |
|   o-> castep_report_storage                       1           1      0.00s  |
|   o-> electronic_minimisation                     1           1    702.90s  |
|   o-> model_write                                 1           1      5.99s  |
+-----------------------------------------------------------------------------+
|   o-- check_elec_ground_state                     1                         |
|  /                                                                          |
| O-> electronic_minimisation                       1           1    702.90s  |
|    /                                                                        |
|   o-> comms_gcopy_logical                         1           1      0.00s  |
|   o-> comms_gcopy_character                       1           1      0.00s  |
|   o-> electronic_initialise                       1           1      0.17s  |
|   o-> electronic_prepare_H                       13          13     87.10s  |
|   o-> electronic_apply_H_energy_local             1           1      0.02s  |
|   o-> electronic_write_energies                  15          15      0.00s  |
|   o-> electronic_store_energy                    15          15      0.00s  |
|   o-> electronic_scf_banner                       1           1      0.00s  |
|   o-> electronic_write_scf_energies              15          15      0.02s  |
|   o-> hamiltonian_diagonalise_ks                 14          14    513.67s  |
|   o-> electronic_find_occupancies                14          14      1.62s  |
|   o-> electronic_check_occupancies               14          14      0.00s  |
|   o-> electronic_apply_H_energy_eigen            14          14      8.02s  |
|   o-> electronic_write_spin_density              14          14      0.00s  |
|   o-> electronic_dump                            14          14      0.00s  |
|   o-> hubbard_calculate_occ_matrix               12          12      0.00s  |
|   o-> hubbard_mix_occ_matrix                     12          12      0.00s  |
|   o-> density_calculate_soft_wvfn                12          12     27.66s  |
|   o-> density_symmetrise                         12          12      2.63s  |
|   o-> density_augment                            12          12     26.86s  |
|   o-> dm_mix_density                             12          12      1.78s  |
|   o-> ewald_dipole_corr                           1           1      0.02s  |
|   o-> electronic_scf_footer                       1           1      0.00s  |
|   o-> electronic_write_occupancies                1           1      0.00s  |
|   o-> sedc_energy                                 1           1     33.31s  |
|   o-> dftd_sedc_print_corr_energies               1           1      0.00s  |
|   o-> electronic_finalise                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> check_forces_stresses                         1           1    542.36s  |
|    /                                                                        |
|   o-> firstd_calculate_forces                     1           1     34.38s  |
|   o-> firstd_symmetrise_forces                    1           1      0.00s  |
|   o-> firstd_output_forces                        1           1      0.01s  |
|   o-> firstd_calculate_stress                     1           1    507.96s  |
|   o-> firstd_symmetrise_stress                    1           1      0.00s  |
|   o-> firstd_output_stress                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    14                         |
|  /                                                                          |
| O-> hamiltonian_diagonalise_ks                   14          14    513.67s  |
|    /                                                                        |
|   o-> wave_allocate_slice                        70          70      0.00s  |
|   o-> wave_initialise_slice                     364         364      1.75s  |
|   o-> comms_reduce_gv_logical                    14          14      0.00s  |
|   o-> comms_reduce_bnd_logical                  290         290      0.00s  |
|   o-> hamiltonian_apply_ks                       12          12    131.71s  |
|   o-> wave_diagonalise_H_ks                      14          14     10.20s  |
|   o-> comms_reduce_bnd_real                      14          14      0.00s  |
|   o-> wave_calc_precon                           14          14      0.00s  |
|   o-> nlpot_prepare_precon_ks                    14          14      4.22s  |
|   o-> comms_reduce_bnd_integer                 1052        1052      0.01s  |
|   o-> wave_copy_wv_slice                        224         224      0.40s  |
|   o-> wave_copy_slice_slice                    8444        8444      1.51s  |
|   o-> wave_dot_all_slice_slice                  388         388      3.34s  |
|   o-> hamiltonian_searchspace_ks                276         276     98.09s  |
|   o-> wave_Sorthog_lower_slice_s                276         276      4.81s  |
|   o-> wave_Sorthonormalise_slice                276         276      2.22s  |
|   o-> hamiltonian_apply_slice                   276         276    222.61s  |
|   o-> comms_reduce_bnd_complex                  276         276      0.06s  |
|   o-> algor_diagonalise_complex                 276         276      0.60s  |
|   o-> wave_rotate_slice                         552         552     30.95s  |
|   o-> wave_copy_slice_wv                        552         552      0.89s  |
|   o-> comms_copy_gv_logical                     828         828      0.01s  |
|   o-> wave_deallocate_slice                      70          70      0.00s  |
|   o-> wave_kinetic_eigenvalues_wv_ks              2           2      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- check_forces_stresses                       1                         |
|  /                                                                          |
| O-> firstd_calculate_stress                       1           1    507.96s  |
|    /                                                                        |
|   o-> ewald_calculate_stress                      1           1      0.16s  |
|   o-> wave_kinetic_stress_wv                      1           1      0.01s  |
|   o-> locpot_calculate_stress                     1           1      3.72s  |
|   o-> hubbard_calculate_stress                    1           1      0.00s  |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> locpot_calculate                            1           1      1.02s  |
|   o-> nlpot_calculate_stress                      1           1    366.01s  |
|   o-> pot_deallocate                              1           1      0.00s  |
|   o-> sedc_stress                                 1           1    137.04s  |
|   o-> tddft_calculate_stress                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> nlpot_calculate_stress                        1           1    366.01s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                           1           1      0.00s  |
|   o-> nlpot_calculate_stress_r                    1           1    366.01s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_stress                      1                         |
|  /                                                                          |
| O-> nlpot_calculate_stress_r                      1           1    366.01s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                           1           1      0.00s  |
|   o-> wave_beta_phi_wv                            1           1      0.00s  |
|   o-> nlpot_allocate_nl_d                         1           1      0.00s  |
|   o-> nlpot_calculate_d                           1           1      5.25s  |
|   o-> wave_dbetadcell_phi_wv_iks                108         108      8.61s  |
|   o-> nlpot_calculate_dddcell_real                1           1    352.08s  |
|   o-> comms_reduce_gv_real                        1           1      0.00s  |
|   o-> comms_reduce_bnd_real                       1           1      0.00s  |
|   o-> comms_reduce_kp_real                        1           1      0.06s  |
|   o-> nlpot_deallocate_nl_d                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_stress_r                    1                         |
|  /                                                                          |
| O-> nlpot_calculate_dddcell_real                  1           1    352.08s  |
|    /                                                                        |
|   o-> basis_real_to_recip_grid                    1           1      0.05s  |
|   o-> basis_recip_fine_to_half_grid               1           1      0.00s  |
|   o-> ion_int_dQdcell_at_origin_recip          2436        2436    351.48s  |
|   o-> comms_reduce_gv_real                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_dddcell_real             2436                         |
|  /                                                                          |
| O-> ion_int_dQdcell_at_origin_recip            2436        2436    351.48s  |
|    /                                                                        |
|   o-> ion_set_dQdcell_at_origin_recip          2436        2436    326.10s  |
+-----------------------------------------------------------------------------+
|   o-- ion_int_dQdcell_at_origin_recip          2436                         |
|  /                                                                          |
| O-> ion_set_dQdcell_at_origin_recip            2436        2436    326.10s  |
|    /                                                                        |
|   o-> ion_atom_radial_transform                2436        2436      0.01s  |
|   o-> ion_Q_recip_interpolation                3012        3012    196.54s  |
|   o-> ion_apply_and_add_ylm                   59544       59544     44.43s  |
|   o-> ion_apply_and_add_ylmp                  59544       59544     38.64s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                      708                         |
|  /o-- hamiltonian_apply_slice                   276                         |
| |/                                                                          |
| O-> pot_apply_slice                             984         984    258.04s  |
|    /                                                                        |
|   o-> pot_debug_check                           984         984      0.00s  |
|   o-> pot_interpolate                            12          12      1.23s  |
|   o-> wave_spin_type_slice                      984         984      0.00s  |
|   o-> pot_nongamma_apply_slice                  984         984    256.79s  |
+-----------------------------------------------------------------------------+
|   o-- pot_apply_slice                           984                         |
|  /                                                                          |
| O-> pot_nongamma_apply_slice                    984         984    256.79s  |
|    /                                                                        |
|   o-> wave_spin_type_slice                      984         984      0.00s  |
|   o-> wave_recip_to_real_slice                 2591        2591    132.89s  |
|   o-> wave_real_to_recip_slice_slice           2591        2591    110.93s  |
+-----------------------------------------------------------------------------+
|   o-- basis_recip_grid_to_real                  704                         |
|  /o-- basis_real_to_recip_grid                  284                         |
|  /o-- basis_real_fine_to_std_grid                48                         |
|  /o-- basis_recip_red_real_3d_coeffs          19686                         |
|  /o-- basis_real_recip_red_3d_coeffs          19686                         |
|  /o-- basis_recip_reduced_real_one             3840                         |
|  /o-- basis_real_std_to_fine_grid                48                         |
| |/                                                                          |
| O-> comms_transpose                           44296       44296    222.74s  |
|    /                                                                        |
|   o-> comms_transpose_n                       44296       44296    222.37s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                276                         |
|  /                                                                          |
| O-> hamiltonian_apply_slice                     276         276    222.61s  |
|    /                                                                        |
|   o-> wave_initialise_slice                     276         276      0.88s  |
|   o-> pot_apply_slice                           276         276    193.93s  |
|   o-> nlpot_apply_add_slice                     276         276     27.15s  |
|   o-> wave_add_kinetic_energy_slice             276         276      0.65s  |
+-----------------------------------------------------------------------------+
|   o-- comms_transpose                         44296                         |
|  /                                                                          |
| O-> comms_transpose_n                         44296       44296    222.37s  |
|    /                                                                        |
|   o-> comms_transpose_exchange                44296       44296    222.00s  |
+-----------------------------------------------------------------------------+
|   o-- comms_transpose_n                       44296                         |
|  /                                                                          |
| O-> comms_transpose_exchange                  44296       44296    222.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_all_beta_multi_phi_recip              569                         |
|  /o-- ion_beta_add_multi_recip_all             1260                         |
|  /o-- local_dot_all_many_many                   955                         |
|  /o-- local_q_dot_all_many_many_c               553                         |
|  /o-- ion_Q_apply_and_add_wsf                   312                         |
| |/                                                                          |
| O-> algor_matmul_cmplx_cmplx                   3649        3649    209.07s  |
+-----------------------------------------------------------------------------+
|   o-- ion_set_Q_at_origin_recip                 138                         |
|  /o-- ion_set_dQdcell_at_origin_recip          3012                         |
| |/                                                                          |
| O-> ion_Q_recip_interpolation                  3150        3150    197.13s  |
+-----------------------------------------------------------------------------+
|   o-- sedc_energy                                 1                         |
|  /o-- sedc_forces                                 1                         |
|  /o-- sedc_stress                                 1                         |
| |/                                                                          |
| O-> sedc_calculate_properties                     3           3    170.34s  |
|    /                                                                        |
|   o-> comms_reduce_gv_logical                     3           3      0.00s  |
|   o-> comms_reduce_bnd_logical                    3           3      0.00s  |
|   o-> comms_reduce_kp_logical                     3           3      0.00s  |
|   o-> dftd_new_calc                               3           3      0.00s  |
|   o-> dftd_TS_ratios                              1           1     18.31s  |
|   o-> sedc_calculate                              2           2    152.04s  |
+-----------------------------------------------------------------------------+
|   o-- sedc_calculate_properties                   2                         |
|  /                                                                          |
| O-> sedc_calculate                                2           2    152.04s  |
|    /                                                                        |
|   o-> mbd2_fd                                     2           2    152.04s  |
+-----------------------------------------------------------------------------+
|   o-- sedc_calculate                              2                         |
|  /                                                                          |
| O-> mbd2_fd                                       2           2    152.04s  |
|    /                                                                        |
|   o-> cell_copy                                  14          14      0.00s  |
|   o-> mbd2_calculate                             14          14    152.02s  |
|   o-> dftd_store_results                          3           3      0.00s  |
|   o-> cell_recip_lattice                         12          12      0.00s  |
|   o-> cell_calculate_volume                      12          12      0.00s  |
|   o-> cell_reduce_symmetry_perturb               12          12      0.00s  |
|   o-> dftd_new_calc                              12          12      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- mbd2_fd                                    14                         |
|  /                                                                          |
| O-> mbd2_calculate                               14          14    152.02s  |
|    /                                                                        |
|   o-> cell_copy                                  14          14      0.00s  |
|   o-> comms_save_strategy                        14          14      0.00s  |
|   o-> sedc_get_default_params                    14          14      0.00s  |
|   o-> dftd_read_custom_params                    14          14      0.00s  |
|   o-> dftd_check_params                          14          14      0.00s  |
|   o-> scs_calculate                               2           2      5.22s  |
|   o-> cell_kpoints_spacing_wrapped               14          14      0.01s  |
|   o-> cell_unfold_kpoints_inplace_wrapped        14          14      0.00s  |
|   o-> comms_reassign_strategy                    14          14      0.00s  |
|   o-> cell_set_current_kpoints_wrapped           28          28      0.00s  |
|   o-> cell_distribute_kpoints_wrapped            28          28      0.00s  |
|   o-> calc_cfdm_energy                           14          14    146.77s  |
|   o-> comms_restore_strategy                     14          14      0.00s  |
|   o-> comms_gcopy_real                           16          16      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- mbd2_calculate                             14                         |
|  /                                                                          |
| O-> calc_cfdm_energy                             14          14    146.77s  |
|    /                                                                        |
|   o-> calc_cfdm_energy_k_tensor                  14          14    143.67s  |
|   o-> calc_cfdm_energy_freq                      28          28      1.94s  |
|   o-> comms_reduce_kp_real                       16          16      1.16s  |
+-----------------------------------------------------------------------------+
|   o-- calc_cfdm_energy                           14                         |
|  /                                                                          |
| O-> calc_cfdm_energy_k_tensor                    14          14    143.67s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> sedc_stress                                   1           1    137.04s  |
|    /                                                                        |
|   o-> sedc_calculate_properties                   1           1    137.04s  |
+-----------------------------------------------------------------------------+
|   o-- pot_nongamma_apply_slice                 2591                         |
|  /                                                                          |
| O-> wave_recip_to_real_slice                   2591        2591    132.89s  |
|    /                                                                        |
|   o-> basis_recip_red_real_3d_coeffs           2591        2591    120.81s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                 12                         |
|  /                                                                          |
| O-> hamiltonian_apply_ks                         12          12    131.71s  |
|    /                                                                        |
|   o-> wave_allocate_slice                        24          24      0.00s  |
|   o-> wave_initialise_slice                      24          24      0.01s  |
|   o-> wave_copy_wv_slice                        708         708      0.20s  |
|   o-> pot_apply_slice                           708         708     64.11s  |
|   o-> nlpot_apply_add_slice                     708         708     66.95s  |
|   o-> wave_copy_slice_wv                        708         708      0.18s  |
|   o-> wave_add_kinetic_energy_wv_ks              12          12      0.15s  |
|   o-> wave_deallocate_slice                      24          24      0.01s  |
|   o-> wave_dot_wv_wv_ks                          12          12      0.08s  |
+-----------------------------------------------------------------------------+
|   o-- wave_recip_to_real_slice                 2591                         |
|  /                                                                          |
| O-> basis_recip_red_real_3d_coeffs             2591        2591    120.81s  |
|    /                                                                        |
|   o-> comms_transpose                         19686       19686     91.20s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_apply_add_slice_r                   984                         |
|  /o-- nlpot_apply_precon_ES_slice               276                         |
| |/                                                                          |
| O-> ion_beta_add_multi_recip_all               1260        1260    120.11s  |
|    /                                                                        |
|   o-> ion_beta_recip_set                       1260        1260      0.08s  |
|   o-> algor_matmul_cmplx_cmplx                 1260        1260    120.02s  |
+-----------------------------------------------------------------------------+
|   o-- pot_nongamma_apply_slice                 2591                         |
|  /                                                                          |
| O-> wave_real_to_recip_slice_slice             2591        2591    110.93s  |
|    /                                                                        |
|   o-> basis_real_recip_red_3d_coeffs           2591        2591    110.91s  |
+-----------------------------------------------------------------------------+
|   o-- wave_real_to_recip_slice_slice           2591                         |
|  /                                                                          |
| O-> basis_real_recip_red_3d_coeffs             2591        2591    110.91s  |
|    /                                                                        |
|   o-> comms_transpose                         19686       19686     84.78s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                276                         |
|  /                                                                          |
| O-> hamiltonian_searchspace_ks                  276         276     98.09s  |
|    /                                                                        |
|   o-> wave_initialise_slice                     276         276      1.01s  |
|   o-> nlpot_apply_precon_ES_slice               276         276     52.91s  |
|   o-> wave_Sorthogonalise_wv_slice              276         276     44.16s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                      708                         |
|  /o-- hamiltonian_apply_slice                   276                         |
| |/                                                                          |
| O-> nlpot_apply_add_slice                       984         984     94.10s  |
|    /                                                                        |
|   o-> nlpot_apply_add_slice_r                   984         984     94.09s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_apply_add_slice                     984                         |
|  /                                                                          |
| O-> nlpot_apply_add_slice_r                     984         984     94.09s  |
|    /                                                                        |
|   o-> wave_spin_type_slice                      984         984      0.00s  |
|   o-> wave_beta_phi_slice                       984         984      0.07s  |
|   o-> comms_reduce_gv_complex                   984         984      0.10s  |
|   o-> ion_beta_add_multi_recip_all              984         984     93.84s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    13                         |
|  /                                                                          |
| O-> electronic_prepare_H                         13          13     87.10s  |
|    /                                                                        |
|   o-> density_symmetrise                         13          13      2.95s  |
|   o-> locpot_calculate                           13          13     14.47s  |
|   o-> nlpot_calculate_d                          13          13     69.68s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_d_real                  36540                         |
|  /o-- nlpot_calculate_dddR                     7308                         |
| |/                                                                          |
| O-> ion_int_Q_at_origin_recip                 43848       43848     85.22s  |
|    /                                                                        |
|   o-> ion_set_Q_at_origin_recip               43848       43848      1.69s  |
|   o-> algor_re_dot_cmplx_cmplx                43848       43848     82.51s  |
|   o-> comms_reduce_gv_real                     7308        7308      0.23s  |
+-----------------------------------------------------------------------------+
|   o-- ion_int_Q_at_origin_recip               43848                         |
|  /                                                                          |
| O-> algor_re_dot_cmplx_cmplx                  43848       43848     82.51s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_prepare_H                       13                         |
|  /o-- nlpot_calculate_forces_r                    1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
| |/                                                                          |
| O-> nlpot_calculate_d                            15          15     80.18s  |
|    /                                                                        |
|   o-> nlpot_calculate_d_real                     15          15     80.18s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_d                          15                         |
|  /                                                                          |
| O-> nlpot_calculate_d_real                       15          15     80.18s  |
|    /                                                                        |
|   o-> basis_real_to_recip_grid                   15          15      1.07s  |
|   o-> basis_recip_fine_to_half_grid              15          15      0.05s  |
|   o-> ion_int_Q_at_origin_recip               36540       36540     71.35s  |
|   o-> comms_reduce_gv_real                       15          15      0.01s  |
|   o-> comms_reduce_bnd_real                      15          15      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_beta_phi_wv_ks                        17                         |
|  /o-- wave_beta_phi_slice                       552                         |
| |/                                                                          |
| O-> ion_all_beta_multi_phi_recip                569         569     55.06s  |
|    /                                                                        |
|   o-> ion_beta_recip_set                        569         569      0.28s  |
|   o-> algor_matmul_cmplx_cmplx                  569         569     54.09s  |
|   o-> comms_reduce_gv_complex                   569         569      0.57s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_searchspace_ks                276                         |
|  /                                                                          |
| O-> nlpot_apply_precon_ES_slice                 276         276     52.91s  |
|    /                                                                        |
|   o-> wave_beta_phi_slice                       552         552     23.65s  |
|   o-> comms_reduce_gv_complex                   276         276      0.10s  |
|   o-> ion_beta_add_multi_recip_all              276         276     26.28s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_apply_add_slice_r                   984                         |
|  /o-- nlpot_apply_precon_ES_slice               552                         |
|  /o-- wave_coeffs_Sdot_all_wv_slice             276                         |
|  /o-- wave_Sdot_lower_slice                     552                         |
|  /o-- wave_calc_Soverlap_slice                  276                         |
|  /o-- wave_orthonormalise_over_slice            276                         |
| |/                                                                          |
| O-> wave_beta_phi_slice                        2916        2916     47.34s  |
|    /                                                                        |
|   o-> ion_set_projectors                       2916        2916      0.19s  |
|   o-> wave_calc_ps_q_nonzero                   2916        2916      0.01s  |
|   o-> ion_all_beta_multi_phi_recip              552         552     46.20s  |
+-----------------------------------------------------------------------------+
|   o-- ion_set_Q_at_origin_recip                 462                         |
|  /o-- ion_set_dQdcell_at_origin_recip         59544                         |
| |/                                                                          |
| O-> ion_apply_and_add_ylm                     60006       60006     44.67s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_searchspace_ks                276                         |
|  /                                                                          |
| O-> wave_Sorthogonalise_wv_slice                276         276     44.16s  |
|    /                                                                        |
|   o-> wave_Sdot_all_wv_slice                    276         276     33.14s  |
|   o-> wave_orthog_over_wv_slice                 276         276     11.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_set_dQdcell_at_origin_recip         59544                         |
|  /                                                                          |
| O-> ion_apply_and_add_ylmp                    59544       59544     38.64s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_TS_ratios                              1                         |
|  /o-- castep                                      1                         |
| |/                                                                          |
| O-> hirshfeld_calculate                           2           2     36.45s  |
|    /                                                                        |
|   o-> hirshfeld_init                              2           2      0.00s  |
|   o-> cell_allocate                               8           8      0.00s  |
|   o-> cell_copy                                   8           8      0.00s  |
|   o-> comms_gcopy_integer                         2           2      0.00s  |
|   o-> comms_gcopy_logical                         4           4      0.00s  |
|   o-> cell_supercell                              2           2      0.01s  |
|   o-> cell_distribute_kpoints_wrapped             2           2      0.00s  |
|   o-> cell_deallocate                             6           6      0.00s  |
|   o-> basis_radial_to_recip_grid                  8           8      0.11s  |
|   o-> basis_recip_grid_to_real                  216         216     15.25s  |
|   o-> comms_reduce_gv_real                      866         866      0.01s  |
|   o-> comms_gcopy_real                            4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_rotate_wv_ks                          42                         |
|  /o-- wave_orthog_over_wv_slice                 552                         |
|  /o-- wave_orthog_lower_over_slice_s            552                         |
|  /o-- wave_rotate_slice                         828                         |
| |/                                                                          |
| O-> algor_matmul_cmplx_cmplx_3D                1974        1974     35.86s  |
+-----------------------------------------------------------------------------+
|   o-- check_forces_stresses                       1                         |
|  /                                                                          |
| O-> firstd_calculate_forces                       1           1     34.38s  |
|    /                                                                        |
|   o-> ewald_calculate_forces                      1           1      0.14s  |
|   o-> locpot_calculate_forces                     1           1      5.59s  |
|   o-> hubbard_calculate_forces                    1           1      0.00s  |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> locpot_calculate                            1           1      1.02s  |
|   o-> nlpot_calculate_forces                      1           1     27.64s  |
|   o-> pot_deallocate                              1           1      0.00s  |
|   o-> ewald_dipole_corr                           1           1      0.01s  |
|   o-> sedc_forces                                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> sedc_energy                                   1           1     33.31s  |
|    /                                                                        |
|   o-> sedc_calculate_properties                   1           1     33.31s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthogonalise_wv_slice              276                         |
|  /                                                                          |
| O-> wave_Sdot_all_wv_slice                      276         276     33.14s  |
|    /                                                                        |
|   o-> wave_coeffs_Sdot_all_wv_slice             276         276     32.98s  |
|   o-> comms_reduce_gv_complex                   276         276      0.15s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sdot_all_wv_slice                    276                         |
|  /                                                                          |
| O-> wave_coeffs_Sdot_all_wv_slice               276         276     32.98s  |
|    /                                                                        |
|   o-> coeffs_dot_all_many_many                  276         276      7.75s  |
|   o-> wave_beta_phi_wv_ks                       276         276      0.31s  |
|   o-> wave_beta_phi_slice                       276         276     23.50s  |
|   o-> wave_q_dot_all_many_many_c                276         276      1.40s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                552                         |
|  /                                                                          |
| O-> wave_rotate_slice                           552         552     30.95s  |
|    /                                                                        |
|   o-> wave_allocate_slice                       552         552      0.01s  |
|   o-> wave_initialise_slice                     552         552     11.54s  |
|   o-> algor_matmul_cmplx_cmplx_3D               828         828     16.73s  |
|   o-> wave_copy_slice_slice                     552         552      2.31s  |
|   o-> wave_deallocate_slice                     552         552      0.31s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    12                         |
|  /                                                                          |
| O-> density_calculate_soft_wvfn                  12          12     27.66s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                          24          24      0.00s  |
|   o-> density_allocate                           12          12      0.03s  |
|   o-> density_calc_soft_wvfn_real                12          12     27.63s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_forces                     1                         |
|  /                                                                          |
| O-> nlpot_calculate_forces                        1           1     27.64s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                           1           1      0.00s  |
|   o-> nlpot_calculate_forces_r                    1           1     27.64s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_forces                      1                         |
|  /                                                                          |
| O-> nlpot_calculate_forces_r                      1           1     27.64s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                           1           1      0.00s  |
|   o-> wave_beta_phi_wv                            1           1      0.00s  |
|   o-> nlpot_allocate_nl_d                         1           1      0.00s  |
|   o-> nlpot_calculate_d                           1           1      5.24s  |
|   o-> wave_dbetadR_phi_wv                       108         108      6.80s  |
|   o-> nlpot_calculate_dddR                        1           1     15.54s  |
|   o-> comms_reduce_gv_real                        1           1      0.00s  |
|   o-> comms_reduce_bnd_real                       1           1      0.00s  |
|   o-> comms_reduce_kp_real                      108         108      0.05s  |
|   o-> nlpot_deallocate_nl_d                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_calculate_soft_wvfn                12                         |
|  /                                                                          |
| O-> density_calc_soft_wvfn_real                  12          12     27.63s  |
|    /                                                                        |
|   o-> wave_recip_to_real_wv_bks                1920        1920     24.20s  |
|   o-> comms_reduce_kp_real                       12          12      0.18s  |
|   o-> comms_reduce_bnd_real                      12          12      0.01s  |
|   o-> basis_real_std_to_fine_gamma               12          12      1.24s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    12                         |
|  /                                                                          |
| O-> density_augment                              12          12     26.86s  |
|    /                                                                        |
|   o-> density_augment_complex                    12          12     26.86s  |
+-----------------------------------------------------------------------------+
|   o-- density_augment                            12                         |
|  /                                                                          |
| O-> density_augment_complex                      12          12     26.86s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                          12          12      0.00s  |
|   o-> ion_augment_charge_nospin_kp               12          12     26.80s  |
+-----------------------------------------------------------------------------+
|   o-- density_augment_complex                    12                         |
|  /                                                                          |
| O-> ion_augment_charge_nospin_kp                 12          12     26.80s  |
|    /                                                                        |
|   o-> comms_reduce_gv_complex                    12          12      0.00s  |
|   o-> comms_reduce_bnd_complex                   12          12      0.00s  |
|   o-> comms_reduce_kp_complex                    12          12      0.00s  |
|   o-> ion_Q_apply_and_add_wsf                    48          48     25.67s  |
|   o-> basis_recip_half_to_fine_grid              12          12      0.18s  |
|   o-> basis_recip_grid_to_real                   12          12      0.88s  |
+-----------------------------------------------------------------------------+
|   o-- ion_augment_charge_nospin_kp               48                         |
|  /                                                                          |
| O-> ion_Q_apply_and_add_wsf                      48          48     25.67s  |
|    /                                                                        |
|   o-> algor_matmul_cmplx_cmplx                  312         312     18.85s  |
|   o-> ion_set_Q_at_origin_recip                1332        1332      0.01s  |
|   o-> comms_reduce_bnd_complex                   48          48      0.31s  |
+-----------------------------------------------------------------------------+
|   o-- density_initialise_cell                     1                         |
|  /o-- density_symmetrise_gv_parallel             25                         |
|  /o-- locps_calculate_potential                   2                         |
|  /o-- hartree_calculate_potential                16                         |
|  /o-- density_gradient_cmplx                     51                         |
|  /o-- xc_calc_vector_divergence_cmplx            17                         |
|  /o-- ion_augment_charge_nospin_kp               12                         |
|  /o-- dm_mix_density_to_density                  12                         |
|  /o-- hirshfeld_calculate                       216                         |
| |/                                                                          |
| O-> basis_recip_grid_to_real                    352         352     24.71s  |
|    /                                                                        |
|   o-> comms_transpose                           704         704     19.16s  |
+-----------------------------------------------------------------------------+
|   o-- density_calc_soft_wvfn_real              1920                         |
|  /                                                                          |
| O-> wave_recip_to_real_wv_bks                  1920        1920     24.20s  |
|    /                                                                        |
|   o-> basis_recip_reduced_real_one             1920        1920     24.18s  |
+-----------------------------------------------------------------------------+
|   o-- wave_recip_to_real_wv_bks                1920                         |
|  /                                                                          |
| O-> basis_recip_reduced_real_one               1920        1920     24.18s  |
|    /                                                                        |
|   o-> comms_transpose                          3840        3840     18.65s  |
+-----------------------------------------------------------------------------+
|   o-- sedc_calculate_properties                   1                         |
|  /                                                                          |
| O-> dftd_TS_ratios                                1           1     18.31s  |
|    /                                                                        |
|   o-> hirshfeld_calculate                         1           1     18.31s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_prepare_H                       13                         |
|  /o-- firstd_calculate_forces                     1                         |
|  /o-- firstd_calculate_stress                     1                         |
|  /o-- write_local_pot_and_den                     1                         |
| |/                                                                          |
| O-> locpot_calculate                             16          16     17.20s  |
|    /                                                                        |
|   o-> pot_complex_to_real                        16          16      0.01s  |
|   o-> pot_allocate                               16          16      0.00s  |
|   o-> pot_zero                                   32          32      0.08s  |
|   o-> locpot_check_locps_core                    16          16      1.25s  |
|   o-> hartree_calculate_potential                16          16      2.51s  |
|   o-> pot_add                                    32          32      0.05s  |
|   o-> pot_calc_energy_real                       42          42      0.09s  |
|   o-> xc_calculate_potential                     16          16     13.11s  |
|   o-> cell_cart_lattice_to_abc                   16          16      0.00s  |
|   o-> ewald_dipole_corr                          16          16      0.10s  |
|   o-> pot_deallocate                             16          16      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                364                         |
|  /o-- hamiltonian_apply_ks                       24                         |
|  /o-- wave_rotate_wv_ks                          28                         |
|  /o-- hamiltonian_searchspace_ks                276                         |
|  /o-- hamiltonian_apply_slice                   276                         |
|  /o-- wave_rotate_slice                         552                         |
|  /o-- electronic_permute_eigenstates             56                         |
| |/                                                                          |
| O-> wave_initialise_slice                      1576        1576     16.15s  |
|     +- section "initialisation", branch on value of method:-                |
|       z                                        1300        1300     14.60s  |
|         \                                                                   |
|          o-> wave_zero_slice                   1300        1300     14.59s  |
|       Z                                         276         276      0.86s  |
|         \                                                                   |
|   o-> wave_zero_slice                           276         276      0.86s  |
|    /                                                                        |
|   o-> wave_setup_slice                         1576        1576      0.64s  |
+-----------------------------------------------------------------------------+
|   o-- coeffs_dot_all_many_many                  679                         |
|  /o-- coeffs_dot_lower_many_many                276                         |
| |/                                                                          |
| O-> local_dot_all_many_many                     955         955     15.81s  |
|    /                                                                        |
|   o-> algor_matmul_cmplx_cmplx                  955         955     15.75s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_forces_r                    1                         |
|  /                                                                          |
| O-> nlpot_calculate_dddR                          1           1     15.54s  |
|    /                                                                        |
|   o-> basis_real_to_recip_grid                    1           1      0.06s  |
|   o-> basis_recip_fine_to_half_grid               1           1      0.00s  |
|   o-> ion_int_Q_at_origin_recip                7308        7308     13.87s  |
+-----------------------------------------------------------------------------+
|   o-- wave_initialise_slice                    1576                         |
|  /                                                                          |
| O-> wave_zero_slice                            1576        1576     15.45s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           16                         |
|  /o-- xc_calculate_stress_density                 1                         |
| |/                                                                          |
| O-> xc_calculate_potential                       17          17     13.94s  |
|    /                                                                        |
|   o-> pot_allocate                               17          17      0.00s  |
|   o-> pot_zero                                   17          17      0.04s  |
|   o-> xc_gga                                     17          17     13.87s  |
|   o-> pot_add                                    17          17      0.01s  |
|   o-> pot_copy                                   17          17      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- wave_dot_all_wv_wv_ks                      14                         |
|  /o-- wave_dot_all_slice_slice                  388                         |
|  /o-- wave_coeffs_Sdot_all_wv_slice             276                         |
|  /o-- wave_coeffs_Sdot_all_wv_ks                  1                         |
| |/                                                                          |
| O-> coeffs_dot_all_many_many                    679         679     13.93s  |
|    /                                                                        |
|   o-> local_dot_all_many_many                   679         679     13.93s  |
+-----------------------------------------------------------------------------+
|   o-- xc_calculate_potential                     17                         |
|  /                                                                          |
| O-> xc_gga                                       17          17     13.87s  |
|    /                                                                        |
|   o-> density_gradient_cmplx                     17          17      5.28s  |
|   o-> comms_reduce_gv_real                       51          51      0.01s  |
|   o-> xc_calc_vector_divergence_cmplx            17          17      5.52s  |
+-----------------------------------------------------------------------------+
|   o-- check_elec_ground_state                     1                         |
|  /o-- castep                                      1                         |
| |/                                                                          |
| O-> model_write                                   2           2     12.54s  |
|    /                                                                        |
|   o-> model_write_all                             4           4     12.14s  |
|   o-> comms_barrier_farm                          2           2      0.00s  |
|   o-> comms_barrier                               2           2      0.00s  |
|   o-> io_delete_file                              1           1      0.39s  |
+-----------------------------------------------------------------------------+
|   o-- model_write                                 4                         |
|  /                                                                          |
| O-> model_write_all                               4           4     12.14s  |
|    /                                                                        |
|   o-> parameters_dump                             4           4      0.02s  |
|   o-> cell_dump                                   8           8      0.03s  |
|   o-> model_write_occ_eigenvalues                 4           4      0.00s  |
|   o-> density_write                               4           4      0.85s  |
|   o-> wave_write_all                              2           2     11.19s  |
+-----------------------------------------------------------------------------+
|   o-- model_write_all                             2                         |
|  /                                                                          |
| O-> wave_write_all                                2           2     11.19s  |
|    /                                                                        |
|   o-> wave_write_all_par                          2           2     11.19s  |
+-----------------------------------------------------------------------------+
|   o-- wave_write_all                              2                         |
|  /                                                                          |
| O-> wave_write_all_par                            2           2     11.19s  |
|    /                                                                        |
|   o-> comms_gather_kp_integer                     2           2      0.00s  |
|   o-> comms_gather_gv_integer                     8           8      0.00s  |
|   o-> comms_gcopy_integer                         2           2      0.00s  |
|   o-> comms_gather_gv_complex                   466         466      0.02s  |
|   o-> comms_send_integer                         32          32      0.00s  |
|   o-> comms_recv_real                            16          16      0.00s  |
|   o-> comms_recv_integer                         64          64      0.00s  |
|   o-> comms_recv_complex                       3728        3728      0.15s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthogonalise_wv_slice              276                         |
|  /                                                                          |
| O-> wave_orthog_over_wv_slice                   276         276     11.00s  |
|    /                                                                        |
|   o-> algor_matmul_cmplx_cmplx_3D               552         552     10.98s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                 14                         |
|  /                                                                          |
| O-> wave_diagonalise_H_ks                        14          14     10.20s  |
|    /                                                                        |
|   o-> wave_dot_all_wv_wv_ks                      14          14      2.74s  |
|   o-> algor_diagonalise_complex                  14          14      0.30s  |
|   o-> wave_rotate_wv_ks                          28          28      7.16s  |
+-----------------------------------------------------------------------------+
|   o-- wave_calc_Soverlap_wv_ks                    2                         |
|  /o-- nlpot_calc_eigenvals_nkns_r                15                         |
|  /o-- wave_coeffs_Sdot_all_wv_slice             276                         |
|  /o-- wave_beta_phi_wv                            2                         |
|  /o-- wave_coeffs_Sdot_all_wv_ks                  2                         |
| |/                                                                          |
| O-> wave_beta_phi_wv_ks                         297         297      9.93s  |
|    /                                                                        |
|   o-> ion_set_projectors                        297         297      0.31s  |
|   o-> wave_calc_ps_q_nonzero                     17          17      0.00s  |
|   o-> ion_all_beta_multi_phi_recip               17          17      8.87s  |
+-----------------------------------------------------------------------------+
|   o-- dm_density_to_mix_density                  13                         |
|  /o-- density_to_recip                           61                         |
|  /o-- xc_calc_vector_divergence_cmplx            51                         |
|  /o-- nlpot_calculate_d_real                     15                         |
|  /o-- nlpot_calculate_dddR                        1                         |
|  /o-- nlpot_calculate_dddcell_real                1                         |
| |/                                                                          |
| O-> basis_real_to_recip_grid                    142         142      9.81s  |
|    /                                                                        |
|   o-> comms_transpose                           284         284      7.34s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_stress_r                  108                         |
|  /                                                                          |
| O-> wave_dbetadcell_phi_wv_iks                  108         108      8.61s  |
|    /                                                                        |
|   o-> ion_dbetadcell_recip_ion                  108         108      1.49s  |
|   o-> comms_reduce_gv_complex                   108         108      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    14                         |
|  /                                                                          |
| O-> electronic_apply_H_energy_eigen              14          14      8.02s  |
|    /                                                                        |
|   o-> wave_kinetic_eigenvalues_wv_ks             14          14      0.22s  |
|   o-> nlpot_calc_eigenvalues_nkns                14          14      7.69s  |
|   o-> comms_reduce_bnd_real                      42          42      0.00s  |
|   o-> comms_reduce_kp_real                       42          42      0.11s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_apply_H_energy_local             1                         |
|  /o-- electronic_apply_H_energy_eigen            14                         |
| |/                                                                          |
| O-> nlpot_calc_eigenvalues_nkns                  15          15      7.69s  |
|    /                                                                        |
|   o-> nlpot_calc_eigenvals_nkns_r                15          15      7.69s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calc_eigenvalues_nkns                15                         |
|  /                                                                          |
| O-> nlpot_calc_eigenvals_nkns_r                  15          15      7.69s  |
|    /                                                                        |
|   o-> wave_beta_phi_wv_ks                        15          15      7.68s  |
|   o-> comms_reduce_gv_real                       15          15      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_diagonalise_H_ks                      28                         |
|  /                                                                          |
| O-> wave_rotate_wv_ks                            28          28      7.16s  |
|    /                                                                        |
|   o-> wave_allocate_slice                        28          28      0.03s  |
|   o-> wave_initialise_slice                      28          28      0.93s  |
|   o-> algor_matmul_cmplx_cmplx_3D                42          42      5.79s  |
|   o-> wave_copy_slice_wv                         28          28      0.36s  |
|   o-> wave_deallocate_slice                      28          28      0.04s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_forces_r                  108                         |
|  /                                                                          |
| O-> wave_dbetadR_phi_wv                         108         108      6.80s  |
|    /                                                                        |
|   o-> ion_dbetadR_recip_ion                     108         108      0.10s  |
|   o-> comms_reduce_gv_complex                   108         108      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_forces                     1                         |
|  /                                                                          |
| O-> locpot_calculate_forces                       1           1      5.59s  |
|    /                                                                        |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> pot_zero                                    1           1      0.00s  |
|   o-> locps_calculate_forces                      1           1      5.59s  |
|   o-> pot_deallocate                              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate_forces                     1                         |
|  /                                                                          |
| O-> locps_calculate_forces                        1           1      5.59s  |
|    /                                                                        |
|   o-> density_to_recip                            1           1      0.08s  |
|   o-> basis_radial_to_recip_grid                  4           4      0.05s  |
|   o-> basis_scale_recip_grid                      4           4      0.01s  |
|   o-> basis_multiply_recip_grid                   4           4      0.02s  |
|   o-> basis_sum_recip_grid                      324         324      1.22s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_prepare_H                       13                         |
|  /o-- electronic_minimisation                    12                         |
| |/                                                                          |
| O-> density_symmetrise                           25          25      5.57s  |
|    /                                                                        |
|   o-> density_real_to_complex                    25          25      0.23s  |
|   o-> density_symmetrise_gv_parallel             25          25      5.18s  |
|   o-> density_complex_to_real                    25          25      0.16s  |
+-----------------------------------------------------------------------------+
|   o-- xc_gga                                     17                         |
|  /                                                                          |
| O-> xc_calc_vector_divergence_cmplx              17          17      5.52s  |
|    /                                                                        |
|   o-> basis_real_to_recip_grid                   51          51      3.61s  |
|   o-> basis_recip_grid_to_real                   17          17      1.06s  |
+-----------------------------------------------------------------------------+
|   o-- xc_gga                                     17                         |
|  /                                                                          |
| O-> density_gradient_cmplx                       17          17      5.28s  |
|    /                                                                        |
|   o-> density_to_recip                           17          17      1.25s  |
|   o-> basis_recip_grid_to_real                   51          51      3.44s  |
+-----------------------------------------------------------------------------+
|   o-- mbd2_calculate                              2                         |
|  /                                                                          |
| O-> scs_calculate                                 2           2      5.22s  |
|    /                                                                        |
|   o-> algor_invert_real                           4           4      0.03s  |
|   o-> comms_reduce_gv_real                        2           2      0.00s  |
|   o-> comms_reduce_bnd_real                       2           2      0.00s  |
|   o-> comms_reduce_kp_real                        2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_symmetrise                         25                         |
|  /                                                                          |
| O-> density_symmetrise_gv_parallel               25          25      5.18s  |
|    /                                                                        |
|   o-> density_set_grid_lookup                    25          25      0.00s  |
|   o-> comms_gather_all_gv_integer               100         100      0.01s  |
|   o-> cell_factor_group_symmetry_wrapped         25          25      0.00s  |
|   o-> density_to_recip                           25          25      1.82s  |
|   o-> comms_gather_all_gv_complex                25          25      0.30s  |
|   o-> density_symmetrise_grid                    25          25      0.55s  |
|   o-> comms_copy_bnd_complex                     25          25      0.00s  |
|   o-> comms_copy_kp_complex                      25          25      0.34s  |
|   o-> basis_recip_grid_to_real                   25          25      1.89s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                276                         |
|  /                                                                          |
| O-> wave_Sorthog_lower_slice_s                  276         276      4.81s  |
|    /                                                                        |
|   o-> wave_Sdot_lower_slice                     276         276      2.42s  |
|   o-> wave_orthog_lower_over_slice_s            276         276      2.38s  |
+-----------------------------------------------------------------------------+
|   o-- density_symmetrise_gv_parallel             25                         |
|  /o-- hartree_calculate_potential                16                         |
|  /o-- density_gradient_cmplx                     17                         |
|  /o-- locps_calculate_forces                      1                         |
|  /o-- hartree_calculate_stress                    1                         |
|  /o-- locps_calculate_stress                      1                         |
| |/                                                                          |
| O-> density_to_recip                             61          61      4.56s  |
|    /                                                                        |
|   o-> basis_real_to_recip_grid                   61          61      4.12s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                 14                         |
|  /                                                                          |
| O-> nlpot_prepare_precon_ks                      14          14      4.22s  |
|    /                                                                        |
|   o-> ion_set_projectors                         14          14      0.00s  |
|   o-> comms_copy_gv_logical                      14          14      0.00s  |
|   o-> comms_copy_bnd_logical                     14          14      0.00s  |
|   o-> ion_beta_beta_recip_cmplx                   5           5      4.15s  |
|   o-> comms_reduce_gv_complex                     5           5      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_prepare_precon_ks                     5                         |
|  /                                                                          |
| O-> ion_beta_beta_recip_cmplx                     5           5      4.15s  |
|    /                                                                        |
|   o-> ion_beta_recip_set                          5           5      0.00s  |
|   o-> comms_reduce_bnd_complex                    5           5      0.02s  |
|   o-> comms_reduce_gv_complex                     5           5      0.05s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks               8444                         |
|  /o-- wave_rotate_slice                         552                         |
|  /o-- electronic_permute_eigenstates            676                         |
| |/                                                                          |
| O-> wave_copy_slice_slice                      9672        9672      3.83s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> locpot_calculate_stress                       1           1      3.72s  |
|    /                                                                        |
|   o-> hartree_calculate_stress                    1           1      0.09s  |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> pot_zero                                    1           1      0.00s  |
|   o-> locps_calculate_stress                      1           1      2.81s  |
|   o-> xc_calculate_stress_density                 1           1      0.82s  |
|   o-> pot_deallocate                              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> popn_calculate_mulliken                       1           1      3.60s  |
|    /                                                                        |
|   o-> popn_create_lcao_basis_allkpt               1           1      1.27s  |
|   o-> popn_calculate_matrices                     1           1      1.68s  |
|   o-> wave_deallocate_wv                          1           1      0.00s  |
|   o-> popn_calculate_charge_spilling              1           1      0.04s  |
|   o-> popn_calculate_density_matrix               1           1      0.02s  |
|   o-> popn_calculate_populations                  1           1      0.58s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                388                         |
|  /                                                                          |
| O-> wave_dot_all_slice_slice                    388         388      3.34s  |
|    /                                                                        |
|   o-> coeffs_dot_all_many_many                  388         388      3.21s  |
|   o-> comms_reduce_gv_complex                   388         388      0.09s  |
+-----------------------------------------------------------------------------+
|   o-- basis_calculate_cut_off                     1                         |
|  /o-- castep_report_storage                       2                         |
|  /o-- ewald_calculate_energy                      3                         |
|  /o-- electronic_apply_H_energy_local             2                         |
|  /o-- electronic_find_fermi_free                651                         |
|  /o-- electronic_occupancy_update                42                         |
|  /o-- electronic_entropy_correction              14                         |
|  /o-- electronic_apply_H_energy_eigen            42                         |
|  /o-- density_calc_soft_wvfn_real                12                         |
|  /o-- scs_calculate                               2                         |
|  /o-- calc_cfdm_energy                           16                         |
|  /o-- ewald_calculate_forces                      3                         |
|  /o-- nlpot_calculate_forces_r                  108                         |
|  /o-- ewald_calculate_stress                      3                         |
|  /o-- wave_kinetic_stress_wv                      1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
|  /o-- popn_calculate_populations                  1                         |
| |/                                                                          |
| O-> comms_reduce_kp_real                        904         904      3.17s  |
|    /                                                                        |
|   o-> comms_reduce_array_real                   136         136      0.38s  |
+-----------------------------------------------------------------------------+
|   o-- locps_calculate_potential                   8                         |
|  /o-- locps_calculate_stress                      4                         |
| |/                                                                          |
| O-> locps_structure_factor                       12          12      3.02s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate_stress                     1                         |
|  /                                                                          |
| O-> locps_calculate_stress                        1           1      2.81s  |
|    /                                                                        |
|   o-> locps_calculate_non_coulomb                 1           1      0.00s  |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> locps_calculate_potential                   1           1      1.25s  |
|   o-> pot_calc_energy_real                        1           1      0.00s  |
|   o-> pot_deallocate                              1           1      0.00s  |
|   o-> density_to_recip                            1           1      0.08s  |
|   o-> basis_radial_to_recip_grid                  4           4      0.05s  |
|   o-> basis_scale_recip_grid                      4           4      0.01s  |
|   o-> locps_structure_factor                      4           4      1.01s  |
|   o-> basis_multiply_recip_grid                   4           4      0.02s  |
|   o-> basis_sum_recip_grid                       24          24      0.10s  |
+-----------------------------------------------------------------------------+
|   o-- wave_diagonalise_H_ks                      14                         |
|  /                                                                          |
| O-> wave_dot_all_wv_wv_ks                        14          14      2.74s  |
|    /                                                                        |
|   o-> coeffs_dot_all_many_many                   14          14      2.73s  |
|   o-> comms_reduce_gv_complex                    14          14      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> model_initialise                              1           1      2.58s  |
|    /                                                                        |
|   o-> model_reset                                 1           1      0.00s  |
|   o-> cell_allocate                               2           2      0.00s  |
|   o-> cell_copy                                   2           2      0.00s  |
|   o-> wave_allocate_wv                            1           1      0.00s  |
|   o-> wave_initialise_wv                          1           1      1.06s  |
|   o-> density_allocate                            1           1      0.00s  |
|   o-> density_initialise_cell                     1           1      1.52s  |
|   o-> model_init_occ_eigenvalues_ef               1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           16                         |
|  /                                                                          |
| O-> hartree_calculate_potential                  16          16      2.51s  |
|    /                                                                        |
|   o-> hartree_check_inv_gsqr                     16          16      0.01s  |
|   o-> density_to_recip                           16          16      1.25s  |
|   o-> comms_reduce_gv_real                       16          16      0.00s  |
|   o-> basis_recip_grid_to_real                   16          16      1.13s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_check_locps_core                     1                         |
|  /o-- locps_calculate_stress                      1                         |
| |/                                                                          |
| O-> locps_calculate_potential                     2           2      2.49s  |
|    /                                                                        |
|   o-> pot_zero                                    2           2      0.00s  |
|   o-> basis_radial_to_recip_grid                  8           8      0.12s  |
|   o-> basis_scale_recip_grid                      8           8      0.03s  |
|   o-> locps_structure_factor                      8           8      2.01s  |
|   o-> basis_multiply_recip_grid                   8           8      0.05s  |
|   o-> basis_recip_grid_to_real                    2           2      0.15s  |
|   o-> basis_scale_real_grid                       2           2      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthog_lower_slice_s                276                         |
|  /                                                                          |
| O-> wave_Sdot_lower_slice                       276         276      2.42s  |
|    /                                                                        |
|   o-> coeffs_dot_lower_many_many                276         276      1.88s  |
|   o-> wave_beta_phi_slice                       552         552      0.07s  |
|   o-> wave_q_dot_lower_many_many_c              276         276      0.40s  |
|   o-> comms_reduce_gv_complex                   276         276      0.05s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthog_lower_slice_s                276                         |
|  /                                                                          |
| O-> wave_orthog_lower_over_slice_s              276         276      2.38s  |
|    /                                                                        |
|   o-> algor_matmul_cmplx_cmplx_3D               552         552      2.36s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                276                         |
|  /                                                                          |
| O-> wave_Sorthonormalise_slice                  276         276      2.22s  |
|    /                                                                        |
|   o-> wave_calc_Soverlap_slice                  276         276      1.21s  |
|   o-> wave_orthonormalise_over_slice            276         276      0.99s  |
|   o-> comms_reduce_bnd_integer                  276         276      0.00s  |
|   o-> comms_reduce_gv_integer                   276         276      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- wave_diagonalise_H_ks                      14                         |
|  /o-- hamiltonian_diagonalise_ks                276                         |
|  /o-- calc_cfdm_energy_freq                      28                         |
| |/                                                                          |
| O-> algor_diagonalise_complex                   318         318      2.18s  |
|    /                                                                        |
|   o-> algor_diagonalise_hermitian_lapack        318         318      2.18s  |
+-----------------------------------------------------------------------------+
|   o-- algor_diagonalise_complex                 318                         |
|  /                                                                          |
| O-> algor_diagonalise_hermitian_lapack          318         318      2.18s  |
|    /                                                                        |
|   o-> comms_copy_gv_complex                     318         318      0.04s  |
|   o-> comms_copy_gv_real                        318         318      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> basis_initialise                              1           1      2.06s  |
|    /                                                                        |
|   o-> cell_copy                                   1           1      0.00s  |
|   o-> basis_utils_prime_factors                   6           6      0.00s  |
|   o-> basis_distribute_grids                      1           1      1.89s  |
|   o-> basis_map_standard_to_fine                  1           1      0.00s  |
|   o-> basis_map_fine_recip_half_full              1           1      0.00s  |
|   o-> basis_assign_grid_coordinates               1           1      0.03s  |
|   o-> basis_count_plane_waves                     1           1      0.14s  |
|   o-> basis_assign_plane_wave_indexes             1           1      0.00s  |
|   o-> basis_assign_pw_gvectors                    1           1      0.00s  |
|   o-> basis_calculate_cut_off                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- calc_cfdm_energy                           28                         |
|  /                                                                          |
| O-> calc_cfdm_energy_freq                        28          28      1.94s  |
|    /                                                                        |
|   o-> algor_diagonalise_complex                  28          28      1.28s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_distribute_grids                        1           1      1.89s  |
|    /                                                                        |
|   o-> basis_utils_sort_columns                    1           1      0.00s  |
|   o-> comms_gather_gv_integer                    12          12      0.08s  |
|   o-> comms_gather_kp_integer                    12          12      0.02s  |
|   o-> comms_map_transpose                         2           2      0.25s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sdot_lower_slice                     276                         |
|  /                                                                          |
| O-> coeffs_dot_lower_many_many                  276         276      1.88s  |
|    /                                                                        |
|   o-> local_dot_all_many_many                   276         276      1.88s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /o-- electronic_minimisation                    12                         |
| |/                                                                          |
| O-> dm_mix_density                               13          13      1.88s  |
|    /                                                                        |
|   o-> dm_mix_density_pulay                       13          13      1.88s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density                             13                         |
|  /                                                                          |
| O-> dm_mix_density_pulay                         13          13      1.88s  |
|    /                                                                        |
|   o-> dm_initialise                               1           1      0.02s  |
|   o-> dm_density_to_mix_density                  12          12      0.88s  |
|   o-> dm_mix_density_kerker                       1           1      0.15s  |
|   o-> dm_mix_density_copy                        55          55      0.00s  |
|   o-> dm_mix_density_add                        176         176      0.01s  |
|   o-> dm_mix_density_dot                        363         363      0.01s  |
|   o-> dm_apply_kerker                            11          11      0.00s  |
|   o-> dm_mix_density_to_density                  11          11      0.81s  |
+-----------------------------------------------------------------------------+
|   o-- wave_q_dot_all_many_many_c                277                         |
|  /o-- wave_q_dot_lower_many_many_c              276                         |
| |/                                                                          |
| O-> local_q_dot_all_many_many_c                 553         553      1.81s  |
|    /                                                                        |
|   o-> weight_beta_phi_many_cmplx                553         553      1.29s  |
|   o-> algor_matmul_cmplx_cmplx                  553         553      0.37s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthonormalise_wv_ks                  1                         |
|  /o-- popn_calculate_matrices                     1                         |
| |/                                                                          |
| O-> wave_calc_Soverlap_wv_ks                      2           2      1.78s  |
|    /                                                                        |
|   o-> coeffs_dot_all_self                         2           2      0.33s  |
|   o-> wave_beta_phi_wv_ks                         2           2      1.42s  |
|   o-> wave_q_dot_all_self_c                       2           2      0.02s  |
|   o-> comms_reduce_gv_complex                     2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              4                         |
|  /o-- ion_generate_orbitals                       4                         |
| |/                                                                          |
| O-> ion_atom_pseudo_scf                           8           8      1.72s  |
|    /                                                                        |
|   o-> ion_atom_init_pseudo_basis                  8           8      0.02s  |
|   o-> ion_atom_init_pseudo_atom                   8           8      0.18s  |
|   o-> ion_atom_init_pseudo_H                      8           8      0.00s  |
|   o-> ion_atom_ps_diag                          194         194      1.10s  |
|   o-> ion_atom_set_pseudo_H                     194         194      0.04s  |
|   o-> ion_atom_regin                          30015       30015      0.05s  |
|   o-> ion_atom_basis_pseudo_dealloc               8           8      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_int_Q_at_origin_recip               43848                         |
|  /o-- ion_Q_apply_and_add_wsf                  1332                         |
| |/                                                                          |
| O-> ion_set_Q_at_origin_recip                 45180       45180      1.71s  |
|    /                                                                        |
|   o-> ion_generate_QLnm                         138         138      0.02s  |
|   o-> ion_atom_radial_transform                 142         142      0.35s  |
|   o-> ion_Q_recip_interpolation                 138         138      0.59s  |
|   o-> ion_apply_and_add_ylm                     462         462      0.25s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_mulliken                     1                         |
|  /                                                                          |
| O-> popn_calculate_matrices                       1           1      1.68s  |
|    /                                                                        |
|   o-> wave_calc_Soverlap_wv_ks                    1           1      0.88s  |
|   o-> wave_Sdot_all_wv_wv_ks                      1           1      0.78s  |
|   o-> algor_invert_complex                        1           1      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    14                         |
|  /                                                                          |
| O-> electronic_find_occupancies                  14          14      1.62s  |
|    /                                                                        |
|   o-> electronic_order_eigenvalues               14          14      0.00s  |
|   o-> electronic_permute_eigenstates             28          28      0.07s  |
|   o-> electronic_occupancy_update                14          14      1.55s  |
+-----------------------------------------------------------------------------+
|   o-- comms_reduce_gv_complex                  3312                         |
|  /o-- comms_reduce_bnd_complex                  342                         |
|  /o-- comms_reduce_kp_complex                    13                         |
| |/                                                                          |
| O-> comms_reduce_array_complex                 3667        3667      1.59s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_find_occupancies                14                         |
|  /                                                                          |
| O-> electronic_occupancy_update                  14          14      1.55s  |
|    /                                                                        |
|   o-> electronic_find_fermi_energy               14          14      1.54s  |
|   o-> comms_reduce_bnd_real                      42          42      0.00s  |
|   o-> comms_reduce_kp_real                       42          42      0.00s  |
|   o-> electronic_entropy_correction              14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_occupancy_update                14                         |
|  /                                                                          |
| O-> electronic_find_fermi_energy                 14          14      1.54s  |
|    /                                                                        |
|   o-> electronic_find_fermi_free                 14          14      1.54s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_find_fermi_energy               14                         |
|  /                                                                          |
| O-> electronic_find_fermi_free                   14          14      1.54s  |
|    /                                                                        |
|   o-> parameters_band_degeneracy                665         665      0.00s  |
|   o-> algor_broadening                         3262        3262      0.00s  |
|   o-> comms_reduce_bnd_real                     651         651      0.00s  |
|   o-> comms_reduce_kp_real                      651         651      1.52s  |
|   o-> comms_gcopy_logical                       637         637      0.00s  |
|   o-> comms_gcopy_real                          623         623      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_initialise                            1                         |
|  /                                                                          |
| O-> density_initialise_cell                       1           1      1.52s  |
|    /                                                                        |
|   o-> density_zero                                1           1      0.00s  |
|   o-> basis_radial_to_recip_grid                  4           4      0.05s  |
|   o-> basis_recip_grid_to_real                    1           1      0.07s  |
|   o-> comms_reduce_gv_real                        1           1      0.00s  |
|   o-> density_complex_to_real                     1           1      0.00s  |
|   o-> comms_copy_bnd_real                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_dbetadcell_phi_wv_iks                108                         |
|  /                                                                          |
| O-> ion_dbetadcell_recip_ion                    108         108      1.49s  |
|    /                                                                        |
|   o-> ion_cc_structure_factor                   108         108      0.02s  |
|   o-> ion_beta_recip_interpolation              600         600      0.48s  |
|   o-> ion_apply_ylmp                           3600        3600      0.07s  |
|   o-> basis_multiply_recip_reduced             3600        3600      0.07s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                      708                         |
|  /o-- wave_rotate_wv_ks                          28                         |
|  /o-- hamiltonian_diagonalise_ks                552                         |
|  /o-- electronic_permute_eigenstates            676                         |
| |/                                                                          |
| O-> wave_copy_slice_wv                         1964        1964      1.43s  |
+-----------------------------------------------------------------------------+
|   o-- wave_coeffs_Sdot_all_wv_slice             276                         |
|  /o-- wave_coeffs_Sdot_all_wv_ks                  1                         |
| |/                                                                          |
| O-> wave_q_dot_all_many_many_c                  277         277      1.41s  |
|    /                                                                        |
|   o-> local_q_dot_all_many_many_c               277         277      1.41s  |
+-----------------------------------------------------------------------------+
|   o-- wave_calc_Soverlap_wv_ks                    2                         |
|  /o-- wave_calc_Soverlap_slice                  276                         |
| |/                                                                          |
| O-> coeffs_dot_all_self                         278         278      1.33s  |
|    /                                                                        |
|   o-> local_dot_all_self                        278         278      1.33s  |
+-----------------------------------------------------------------------------+
|   o-- coeffs_dot_all_self                       278                         |
|  /                                                                          |
| O-> local_dot_all_self                          278         278      1.33s  |
+-----------------------------------------------------------------------------+
|   o-- locps_calculate_forces                    324                         |
|  /o-- locps_calculate_stress                     24                         |
| |/                                                                          |
| O-> basis_sum_recip_grid                        348         348      1.31s  |
|    /                                                                        |
|   o-> comms_reduce_gv_complex                   348         348      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- local_q_dot_all_many_many_c               553                         |
|  /                                                                          |
| O-> weight_beta_phi_many_cmplx                  553         553      1.29s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_mulliken                     1                         |
|  /                                                                          |
| O-> popn_create_lcao_basis_allkpt                 1           1      1.27s  |
|    /                                                                        |
|   o-> popn_check_lcao_basis                       1           1      0.00s  |
|   o-> ion_generate_orbitals                       4           4      1.20s  |
|   o-> wave_allocate_wv                            1           1      0.00s  |
|   o-> wave_initialise_wv                          1           1      0.02s  |
|   o-> basis_radial_to_recip_reduced             300         300      0.04s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           16                         |
|  /                                                                          |
| O-> locpot_check_locps_core                      16          16      1.25s  |
|    /                                                                        |
|   o-> comms_reduce_gv_logical                    32          32      0.01s  |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> pot_zero                                    1           1      0.00s  |
|   o-> locps_calculate_potential                   1           1      1.24s  |
+-----------------------------------------------------------------------------+
|   o-- wave_initialise_wv                          3                         |
|  /o-- wave_beta_phi_wv_ks                       297                         |
|  /o-- ion_beta_recip_set                       1942                         |
|  /o-- wave_calc_storage_wv                        3                         |
|  /o-- wave_calc_storage_slice                     3                         |
|  /o-- wave_calc_storage_bnd                       1                         |
|  /o-- wave_setup_slice                         1576                         |
|  /o-- wave_beta_phi_slice                      2916                         |
|  /o-- nlpot_prepare_precon_ks                    14                         |
| |/                                                                          |
| O-> ion_set_projectors                         6755        6755      1.25s  |
|    /                                                                        |
|   o-> comms_reduce_gv_logical                  6755        6755      1.21s  |
|   o-> comms_reduce_gv_integer                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_all_beta_multi_phi_recip              569                         |
|  /o-- wave_calc_Soverlap_wv_ks                    2                         |
|  /o-- nlpot_apply_add_slice_r                   984                         |
|  /o-- wave_dot_wv_wv_ks                          12                         |
|  /o-- wave_dot_all_wv_wv_ks                      14                         |
|  /o-- ion_beta_beta_recip_cmplx                   5                         |
|  /o-- nlpot_prepare_precon_ks                     5                         |
|  /o-- wave_dot_all_slice_slice                  388                         |
|  /o-- nlpot_apply_precon_ES_slice               276                         |
|  /o-- wave_Sdot_all_wv_slice                    276                         |
|  /o-- wave_Sdot_lower_slice                     276                         |
|  /o-- wave_calc_Soverlap_slice                  276                         |
|  /o-- ion_augment_charge_nospin_kp               12                         |
|  /o-- dm_mix_density_dot                        364                         |
|  /o-- basis_sum_recip_grid                      348                         |
|  /o-- wave_dbetadR_phi_wv                       108                         |
|  /o-- wave_dbetadcell_phi_wv_iks                108                         |
|  /o-- wave_Sdot_all_wv_wv_ks                      1                         |
| |/                                                                          |
| O-> comms_reduce_gv_complex                    4024        4024      1.24s  |
|    /                                                                        |
|   o-> comms_reduce_array_complex               3312        3312      1.19s  |
+-----------------------------------------------------------------------------+
|   o-- density_calc_soft_wvfn_real                12                         |
|  /                                                                          |
| O-> basis_real_std_to_fine_gamma                 12          12      1.24s  |
|    /                                                                        |
|   o-> basis_real_std_to_fine_grid                12          12      1.15s  |
+-----------------------------------------------------------------------------+
|   o-- pot_apply_slice                            12                         |
|  /                                                                          |
| O-> pot_interpolate                              12          12      1.23s  |
|    /                                                                        |
|   o-> pot_debug_check                            12          12      0.00s  |
|   o-> basis_real_fine_to_std_grid                12          12      1.16s  |
+-----------------------------------------------------------------------------+
|   o-- ion_set_projectors                       6755                         |
|  /o-- ewald_calculate_energy                      2                         |
|  /o-- locpot_check_locps_core                    32                         |
|  /o-- electronic_store_energy                    30                         |
|  /o-- hamiltonian_diagonalise_ks                 14                         |
|  /o-- electronic_check_occupancies               14                         |
|  /o-- sedc_calculate_properties                   3                         |
|  /o-- ewald_calculate_forces                      2                         |
|  /o-- ewald_calculate_stress                      2                         |
|  /o-- bib_output                                  1                         |
| |/                                                                          |
| O-> comms_reduce_gv_logical                    6855        6855      1.22s  |
|    /                                                                        |
|   o-> comms_lcopy                                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthonormalise_slice                276                         |
|  /                                                                          |
| O-> wave_calc_Soverlap_slice                    276         276      1.21s  |
|    /                                                                        |
|   o-> coeffs_dot_all_self                       276         276      1.00s  |
|   o-> wave_beta_phi_slice                       276         276      0.04s  |
|   o-> wave_q_dot_all_self_c                     276         276      0.14s  |
|   o-> comms_reduce_gv_complex                   276         276      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- popn_create_lcao_basis_allkpt               4                         |
|  /                                                                          |
| O-> ion_generate_orbitals                         4           4      1.20s  |
|    /                                                                        |
|   o-> ion_atom_allocate_pspot                     4           4      0.00s  |
|   o-> ion_set_psp                                 4           4      0.00s  |
|   o-> ion_atom_pseudo_scf                         4           4      1.19s  |
|   o-> ion_atom_deallocate_pspot                   4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_interpolate                            12                         |
|  /                                                                          |
| O-> basis_real_fine_to_std_grid                  12          12      1.16s  |
|    /                                                                        |
|   o-> comms_transpose                            48          48      0.82s  |
+-----------------------------------------------------------------------------+
|   o-- basis_real_std_to_fine_gamma               12                         |
|  /                                                                          |
| O-> basis_real_std_to_fine_grid                  12          12      1.15s  |
|    /                                                                        |
|   o-> comms_transpose                            48          48      0.80s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_pseudo_scf                       194                         |
|  /                                                                          |
| O-> ion_atom_ps_diag                            194         194      1.10s  |
|    /                                                                        |
|   o-> ion_atom_regin                         179100      179100      0.49s  |
|   o-> ion_atom_rectoreal                        890         890      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- model_initialise                            1                         |
|  /o-- electronic_initialise                       1                         |
|  /o-- popn_create_lcao_basis_allkpt               1                         |
| |/                                                                          |
| O-> wave_initialise_wv                            3           3      1.10s  |
|     +- section "initialisation", branch on value of method:-                |
|       R                                           1           1      1.06s  |
|         \                                                                   |
|          o-> wave_spin_type_wv                    1           1      0.00s  |
|          o-> wave_prepare_init_wvfn               1           1      0.00s  |
|          o-> algor_set_random_seed                1           1      0.00s  |
|          o-> algor_uniform_random_array         233         233      0.01s  |
|          o-> wave_Sorthonormalise_wv              1           1      1.03s  |
|       Z                                           2           2      0.04s  |
|    /                                                                        |
|   o-> ion_set_projectors                          3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_initialise_wv                          1                         |
|  /                                                                          |
| O-> wave_Sorthonormalise_wv                       1           1      1.03s  |
|    /                                                                        |
|   o-> wave_Sorthonormalise_wv_ks                  1           1      1.03s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthonormalise_wv                     1                         |
|  /                                                                          |
| O-> wave_Sorthonormalise_wv_ks                    1           1      1.03s  |
|    /                                                                        |
|   o-> wave_calc_Soverlap_wv_ks                    1           1      0.89s  |
|   o-> wave_orthonormalise_over_wv_ks              1           1      0.14s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthonormalise_slice                276                         |
|  /                                                                          |
| O-> wave_orthonormalise_over_slice              276         276      0.99s  |
|    /                                                                        |
|   o-> wave_beta_phi_slice                       276         276      0.00s  |
|   o-> algor_invert_complex                      276         276      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> write_local_pot_and_den                       1           1      0.98s  |
|    /                                                                        |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> locpot_calculate                            1           1      0.69s  |
|   o-> pot_write                                   1           1      0.28s  |
|   o-> pot_deallocate                              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_pulay                       12                         |
|  /o-- dm_mix_density_kerker                       1                         |
| |/                                                                          |
| O-> dm_density_to_mix_density                    13          13      0.96s  |
|    /                                                                        |
|   o-> basis_real_to_recip_grid                   13          13      0.89s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_kerker                       1                         |
|  /o-- dm_mix_density_pulay                       11                         |
| |/                                                                          |
| O-> dm_mix_density_to_density                    12          12      0.88s  |
|    /                                                                        |
|   o-> basis_recip_grid_to_real                   12          12      0.85s  |
+-----------------------------------------------------------------------------+
|   o-- model_write_all                             4                         |
|  /                                                                          |
| O-> density_write                                 4           4      0.85s  |
|    /                                                                        |
|   o-> density_allocate                            4           4      0.00s  |
|   o-> density_copy                                4           4      0.00s  |
|   o-> density_write_parallel                      4           4      0.84s  |
|   o-> density_deallocate                          4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_write                               4                         |
|  /                                                                          |
| O-> density_write_parallel                        4           4      0.84s  |
|    /                                                                        |
|   o-> comms_gather_gv_integer                    12          12      0.00s  |
|   o-> comms_gather_gv_real                        4           4      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate_stress                     1                         |
|  /                                                                          |
| O-> xc_calculate_stress_density                   1           1      0.82s  |
|    /                                                                        |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> xc_calculate_potential                      1           1      0.82s  |
|   o-> pot_deallocate                              1           1      0.00s  |
|   o-> xc_calculate_stress_correction              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_matrices                     1                         |
|  /                                                                          |
| O-> wave_Sdot_all_wv_wv_ks                        1           1      0.78s  |
|    /                                                                        |
|   o-> wave_coeffs_Sdot_all_wv_ks                  1           1      0.78s  |
|   o-> comms_reduce_gv_complex                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sdot_all_wv_wv_ks                      1                         |
|  /                                                                          |
| O-> wave_coeffs_Sdot_all_wv_ks                    1           1      0.78s  |
|    /                                                                        |
|   o-> coeffs_dot_all_many_many                    1           1      0.25s  |
|   o-> wave_beta_phi_wv_ks                         2           2      0.52s  |
|   o-> wave_q_dot_all_many_many_c                  1           1      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- comms_reduce_gv_real                      364                         |
|  /o-- comms_reduce_kp_real                      136                         |
|  /o-- comms_reduce_bnd_real                      48                         |
| |/                                                                          |
| O-> comms_reduce_array_real                     548         548      0.69s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_slice                   276                         |
|  /                                                                          |
| O-> wave_add_kinetic_energy_slice               276         276      0.65s  |
|    /                                                                        |
|   o-> local_kinetic_energy_many                 276         276      0.45s  |
|   o-> comms_reduce_gv_real                      276         276      0.19s  |
+-----------------------------------------------------------------------------+
|   o-- wave_initialise_slice                    1576                         |
|  /                                                                          |
| O-> wave_setup_slice                           1576        1576      0.64s  |
|    /                                                                        |
|   o-> ion_set_projectors                       1576        1576      0.63s  |
|   o-> wave_band_basis_initialise                514         514      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                      708                         |
|  /o-- hamiltonian_diagonalise_ks                224                         |
|  /o-- electronic_permute_eigenstates            726                         |
| |/                                                                          |
| O-> wave_copy_wv_slice                         1658        1658      0.62s  |
+-----------------------------------------------------------------------------+
|   o-- wave_add_kinetic_energy_wv_ks              12                         |
|  /o-- wave_add_kinetic_energy_slice             276                         |
| |/                                                                          |
| O-> local_kinetic_energy_many                   288         288      0.60s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> ion_initialise                                1           1      0.59s  |
|    /                                                                        |
|   o-> ion_allocate                                1           1      0.00s  |
|   o-> ion_atom_allocate_pspot                     8           8      0.01s  |
|   o-> ion_atom_read_usp                           4           4      0.04s  |
|   o-> ion_set_data                                4           4      0.01s  |
|   o-> ion_atom_deallocate_pspot                   9           9      0.00s  |
|   o-> ion_set_psp                                 4           4      0.00s  |
|   o-> ion_atom_pseudo_scf                         4           4      0.53s  |
|   o-> ion_clebsch_gordan                          1           1      0.00s  |
|   o-> ion_atom_radial_transform                   1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_mulliken                     1                         |
|  /                                                                          |
| O-> popn_calculate_populations                    1           1      0.58s  |
|    /                                                                        |
|   o-> comms_reduce_kp_complex                     1           1      0.00s  |
|   o-> cell_format_atom_label                   1710        1710      0.01s  |
|   o-> cell_verify_mixture_components          11556       11556      0.04s  |
|   o-> comms_reduce_kp_real                        1           1      0.00s  |
|   o-> popn_sort                                   1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_basis                784                         |
|  /o-- ion_atom_init_pseudo_H                    112                         |
|  /o-- ion_atom_pseudo_hartree                   202                         |
|  /o-- ion_atom_pseudo_xc                        202                         |
|  /o-- ion_atom_ps_diag                       179100                         |
|  /o-- ion_atom_calc_pseudo_rho                  194                         |
|  /o-- ion_atom_set_pseudo_H                    2914                         |
|  /o-- ion_atom_pseudo_scf                     30015                         |
| |/                                                                          |
| O-> ion_atom_regin                           213523      213523      0.56s  |
+-----------------------------------------------------------------------------+
|   o-- ion_beta_recip_set                        600                         |
|  /o-- ion_dbetadcell_recip_ion                  600                         |
| |/                                                                          |
| O-> ion_beta_recip_interpolation               1200        1200      0.56s  |
+-----------------------------------------------------------------------------+
|   o-- density_symmetrise_gv_parallel             25                         |
|  /                                                                          |
| O-> density_symmetrise_grid                      25          25      0.55s  |
+-----------------------------------------------------------------------------+
|   o-- basis_calculate_cut_off                     1                         |
|  /o-- density_initialise_cell                     1                         |
|  /o-- castep_report_storage                       2                         |
|  /o-- ewald_calculate_energy                      3                         |
|  /o-- hartree_calculate_potential                16                         |
|  /o-- pot_calc_energy_really_real                43                         |
|  /o-- xc_gga                                     51                         |
|  /o-- nlpot_calculate_d_real                     15                         |
|  /o-- wave_kinetic_eigenvalues_wv_ks             17                         |
|  /o-- nlpot_calc_eigenvals_nkns_r                15                         |
|  /o-- wave_add_kinetic_energy_wv_ks              12                         |
|  /o-- wave_add_kinetic_energy_slice             276                         |
|  /o-- hirshfeld_calculate                       866                         |
|  /o-- scs_calculate                               2                         |
|  /o-- ewald_calculate_forces                      3                         |
|  /o-- ion_int_Q_at_origin_recip                7308                         |
|  /o-- nlpot_calculate_forces_r                    1                         |
|  /o-- ewald_calculate_stress                      3                         |
|  /o-- wave_kinetic_stress_wv                      1                         |
|  /o-- hartree_calculate_stress                    1                         |
|  /o-- nlpot_calculate_dddcell_real                1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
| |/                                                                          |
| O-> comms_reduce_gv_real                       8639        8639      0.51s  |
|    /                                                                        |
|   o-> comms_reduce_array_real                   364         364      0.26s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sdot_lower_slice                     276                         |
|  /                                                                          |
| O-> wave_q_dot_lower_many_many_c                276         276      0.40s  |
|    /                                                                        |
|   o-> local_q_dot_all_many_many_c               276         276      0.40s  |
+-----------------------------------------------------------------------------+
|   o-- wave_orthonormalise_over_wv_ks              1                         |
|  /o-- ion_beta_beta_recip_cmplx                   5                         |
|  /o-- hamiltonian_diagonalise_ks                276                         |
|  /o-- ion_augment_charge_nospin_kp               12                         |
|  /o-- ion_Q_apply_and_add_wsf                    48                         |
| |/                                                                          |
| O-> comms_reduce_bnd_complex                    342         342      0.39s  |
|    /                                                                        |
|   o-> comms_reduce_array_complex                342         342      0.39s  |
+-----------------------------------------------------------------------------+
|   o-- model_write                                 1                         |
|  /                                                                          |
| O-> io_delete_file                                1           1      0.39s  |
+-----------------------------------------------------------------------------+
|   o-- density_initialise_cell                     4                         |
|  /o-- locps_calculate_potential                   8                         |
|  /o-- hirshfeld_calculate                         8                         |
|  /o-- locps_calculate_forces                      4                         |
|  /o-- locps_calculate_stress                      4                         |
| |/                                                                          |
| O-> basis_radial_to_recip_grid                   28          28      0.38s  |
|    /                                                                        |
|   o-> basis_utils_interpolation                  28          28      0.28s  |
|   o-> basis_utils_apply_s_harmonic               28          28      0.10s  |
+-----------------------------------------------------------------------------+
|   o-- ion_all_beta_multi_phi_recip              569                         |
|  /o-- ion_beta_add_multi_recip_all             1260                         |
|  /o-- ion_beta_beta_recip_cmplx                   5                         |
|  /o-- ion_dbetadR_recip_ion                     108                         |
| |/                                                                          |
| O-> ion_beta_recip_set                         1942        1942      0.37s  |
|    /                                                                        |
|   o-> ion_set_projectors                       1942        1942      0.12s  |
|   o-> ion_cc_structure_factor                   108         108      0.01s  |
|   o-> ion_beta_recip_interpolation              600         600      0.08s  |
|   o-> basis_multiply_recip_reduced              600         600      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                       24                         |
|  /o-- wave_rotate_wv_ks                          28                         |
|  /o-- wave_rotate_slice                         552                         |
|  /o-- hamiltonian_diagonalise_ks                 70                         |
|  /o-- electronic_permute_eigenstates             56                         |
| |/                                                                          |
| O-> wave_deallocate_slice                       730         730      0.36s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              1                         |
|  /o-- ion_set_Q_at_origin_recip                 142                         |
|  /o-- ion_set_dQdcell_at_origin_recip          2436                         |
| |/                                                                          |
| O-> ion_atom_radial_transform                  2579        2579      0.36s  |
+-----------------------------------------------------------------------------+
|   o-- density_symmetrise_gv_parallel             25                         |
|  /                                                                          |
| O-> comms_copy_kp_complex                        25          25      0.34s  |
+-----------------------------------------------------------------------------+
|   o-- basis_radial_to_recip_grid                 28                         |
|  /o-- basis_radial_to_recip_reduced             300                         |
| |/                                                                          |
| O-> basis_utils_interpolation                   328         328      0.31s  |
+-----------------------------------------------------------------------------+
|   o-- density_symmetrise_gv_parallel             25                         |
|  /                                                                          |
| O-> comms_gather_all_gv_complex                  25          25      0.30s  |
+-----------------------------------------------------------------------------+
|   o-- write_local_pot_and_den                     1                         |
|  /                                                                          |
| O-> pot_write                                     1           1      0.28s  |
|    /                                                                        |
|   o-> pot_write_parallel                          1           1      0.28s  |
+-----------------------------------------------------------------------------+
|   o-- pot_write                                   1                         |
|  /                                                                          |
| O-> pot_write_parallel                            1           1      0.28s  |
|    /                                                                        |
|   o-> comms_gather_gv_integer                     3           3      0.00s  |
|   o-> comms_gather_gv_real                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_apply_H_energy_local             1                         |
|  /o-- electronic_apply_H_energy_eigen            14                         |
|  /o-- hamiltonian_diagonalise_ks                  2                         |
| |/                                                                          |
| O-> wave_kinetic_eigenvalues_wv_ks               17          17      0.27s  |
|    /                                                                        |
|   o-> comms_reduce_gv_real                       17          17      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_distribute_grids                      2                         |
|  /                                                                          |
| O-> comms_map_transpose                           2           2      0.25s  |
|    /                                                                        |
|   o-> comms_map_transpose_n                       2           2      0.25s  |
+-----------------------------------------------------------------------------+
|   o-- comms_map_transpose                         2                         |
|  /                                                                          |
| O-> comms_map_transpose_n                         2           2      0.25s  |
|    /                                                                        |
|   o-> comms_local_map                             2           2      0.25s  |
+-----------------------------------------------------------------------------+
|   o-- comms_map_transpose_n                       2                         |
|  /                                                                          |
| O-> comms_local_map                               2           2      0.25s  |
+-----------------------------------------------------------------------------+
|   o-- density_symmetrise                         25                         |
|  /                                                                          |
| O-> density_real_to_complex                      25          25      0.23s  |
|    /                                                                        |
|   o-> density_allocate                           25          25      0.09s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_pseudo_scf                         8                         |
|  /                                                                          |
| O-> ion_atom_init_pseudo_atom                     8           8      0.18s  |
|    /                                                                        |
|   o-> ion_atom_resolve_pseudo_cfg                 8           8      0.00s  |
|   o-> ion_atom_locate                          4152        4152      0.01s  |
|   o-> ion_atom_interpolate                     4152        4152      0.01s  |
|   o-> ion_atom_set_pseudo_occ                     8           8      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_augment_charge_nospin_kp               12                         |
|  /                                                                          |
| O-> basis_recip_half_to_fine_grid                12          12      0.18s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> electronic_initialise                         1           1      0.17s  |
|    /                                                                        |
|   o-> density_allocate                            1           1      0.00s  |
|   o-> density_zero                                1           1      0.00s  |
|   o-> wave_allocate_wv                            1           1      0.00s  |
|   o-> wave_initialise_wv                          1           1      0.02s  |
|   o-> ion_real_initialise                         1           1      0.00s  |
|   o-> pot_allocate                                1           1      0.00s  |
|   o-> hubbard_initialise                          1           1      0.00s  |
|   o-> nlpot_allocate_nl_d                         1           1      0.00s  |
|   o-> electronic_restore                          1           1      0.00s  |
|   o-> dm_flush_history                            1           1      0.00s  |
|   o-> dm_mix_density                              1           1      0.10s  |
|   o-> ewald_calculate_energy                      1           1      0.06s  |
|   o-> locps_calculate_non_coulomb                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_initialise                            1                         |
|  /o-- nlxc_initialise                             1                         |
|  /o-- electronic_initialise                       1                         |
|  /o-- dm_initialise                               1                         |
|  /o-- density_real_to_complex                    25                         |
|  /o-- density_complex_to_real                    25                         |
|  /o-- density_calculate_soft_wvfn                12                         |
|  /o-- density_write                               4                         |
| |/                                                                          |
| O-> density_allocate                             70          70      0.17s  |
+-----------------------------------------------------------------------------+
|   o-- density_initialise_cell                     1                         |
|  /o-- density_symmetrise                         25                         |
| |/                                                                          |
| O-> density_complex_to_real                      26          26      0.16s  |
|    /                                                                        |
|   o-> density_allocate                           25          25      0.04s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> ewald_calculate_stress                        1           1      0.16s  |
|    /                                                                        |
|   o-> comms_reduce_kp_logical                     2           2      0.00s  |
|   o-> comms_reduce_bnd_logical                    2           2      0.00s  |
|   o-> comms_reduce_gv_logical                     2           2      0.00s  |
|   o-> ewald_calculate_num_cells                   1           1      0.01s  |
|   o-> comms_reduce_gv_real                        3           3      0.03s  |
|   o-> comms_reduce_bnd_real                       3           3      0.02s  |
|   o-> comms_reduce_kp_real                        3           3      0.06s  |
+-----------------------------------------------------------------------------+
|   o-- ewald_calc_storage                          1                         |
|  /o-- ewald_calculate_energy                      1                         |
|  /o-- ewald_dipole_corr                          18                         |
|  /o-- ewald_calculate_forces                      1                         |
|  /o-- ewald_calculate_stress                      1                         |
| |/                                                                          |
| O-> ewald_calculate_num_cells                    22          22      0.15s  |
|    /                                                                        |
|   o-> comms_gcopy_real                          110         110      0.01s  |
|   o-> comms_gcopy_logical                        66          66      0.00s  |
|   o-> comms_gcopy_character                      22          22      0.00s  |
|   o-> comms_gcopy_integer                       132         132      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_calc_Soverlap_wv_ks                    2                         |
|  /o-- wave_calc_Soverlap_slice                  276                         |
| |/                                                                          |
| O-> wave_q_dot_all_self_c                       278         278      0.15s  |
|    /                                                                        |
|   o-> local_q_dot_all_self_c                    278         278      0.15s  |
+-----------------------------------------------------------------------------+
|   o-- wave_write_all_par                       3728                         |
|  /                                                                          |
| O-> comms_recv_complex                         3728        3728      0.15s  |
+-----------------------------------------------------------------------------+
|   o-- wave_q_dot_all_self_c                     278                         |
|  /                                                                          |
| O-> local_q_dot_all_self_c                      278         278      0.15s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_pulay                        1                         |
|  /                                                                          |
| O-> dm_mix_density_kerker                         1           1      0.15s  |
|    /                                                                        |
|   o-> dm_density_to_mix_density                   1           1      0.07s  |
|   o-> dm_mix_density_copy                         1           1      0.00s  |
|   o-> dm_mix_density_add                          2           2      0.00s  |
|   o-> dm_mix_density_dot                          1           1      0.00s  |
|   o-> dm_apply_kerker                             1           1      0.00s  |
|   o-> dm_mix_density_to_density                   1           1      0.07s  |
+-----------------------------------------------------------------------------+
|   o-- basis_count_plane_waves                     3                         |
|  /o-- ion_set_projectors                          1                         |
|  /o-- wave_prepare_init_wvfn                      2                         |
|  /o-- dm_count_plane_waves                        6                         |
|  /o-- wave_Sorthonormalise_slice                276                         |
| |/                                                                          |
| O-> comms_reduce_gv_integer                     288         288      0.15s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                       12                         |
|  /                                                                          |
| O-> wave_add_kinetic_energy_wv_ks                12          12      0.15s  |
|    /                                                                        |
|   o-> local_kinetic_energy_many                  12          12      0.15s  |
|   o-> comms_reduce_gv_real                       12          12      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_count_plane_waves                       1           1      0.14s  |
|    /                                                                        |
|   o-> comms_reduce_gv_integer                     3           3      0.14s  |
|   o-> comms_reduce_kp_integer                     3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_Sorthonormalise_wv_ks                  1                         |
|  /                                                                          |
| O-> wave_orthonormalise_over_wv_ks                1           1      0.14s  |
|    /                                                                        |
|   o-> comms_reduce_bnd_complex                    1           1      0.00s  |
|   o-> algor_invert_complex                        1           1      0.00s  |
|   o-> comms_copy_gv_complex                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_forces                     1                         |
|  /                                                                          |
| O-> ewald_calculate_forces                        1           1      0.14s  |
|    /                                                                        |
|   o-> comms_reduce_kp_logical                     2           2      0.00s  |
|   o-> comms_reduce_bnd_logical                    2           2      0.00s  |
|   o-> comms_reduce_gv_logical                     2           2      0.00s  |
|   o-> ewald_calculate_num_cells                   1           1      0.01s  |
|   o-> comms_reduce_gv_real                        3           3      0.02s  |
|   o-> comms_reduce_bnd_real                       3           3      0.01s  |
|   o-> comms_reduce_kp_real                        3           3      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           16                         |
|  /o-- electronic_minimisation                     1                         |
|  /o-- firstd_calculate_forces                     1                         |
| |/                                                                          |
| O-> ewald_dipole_corr                            18          18      0.13s  |
|    /                                                                        |
|   o-> ewald_calculate_num_cells                  18          18      0.13s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           32                         |
|  /o-- locpot_check_locps_core                     1                         |
|  /o-- locps_calculate_potential                   2                         |
|  /o-- xc_calculate_potential                     17                         |
|  /o-- locpot_calculate_forces                     1                         |
|  /o-- locpot_calculate_stress                     1                         |
| |/                                                                          |
| O-> pot_zero                                     54          54      0.12s  |
|    /                                                                        |
|   o-> pot_debug_check                            54          54      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> write_eigenvalues                             1           1      0.11s  |
|    /                                                                        |
|   o-> comms_reduce_bnd_real                       1           1      0.00s  |
|   o-> comms_gather_kp_integer                     1           1      0.00s  |
|   o-> global_kpoint_index                         9           9      0.00s  |
|   o-> comms_recv_real                            24          24      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_dbetadR_phi_wv                       108                         |
|  /                                                                          |
| O-> ion_dbetadR_recip_ion                       108         108      0.10s  |
|    /                                                                        |
|   o-> ion_beta_recip_set                        108         108      0.01s  |
|   o-> ion_proj_index                            600         600      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_radial_to_recip_grid                 28                         |
|  /o-- basis_radial_to_recip_reduced             108                         |
| |/                                                                          |
| O-> basis_utils_apply_s_harmonic                136         136      0.10s  |
+-----------------------------------------------------------------------------+
|   o-- locps_calculate_potential                   8                         |
|  /o-- locps_calculate_forces                      4                         |
|  /o-- locps_calculate_stress                      4                         |
| |/                                                                          |
| O-> basis_multiply_recip_grid                    16          16      0.09s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           42                         |
|  /o-- locps_calculate_stress                      1                         |
| |/                                                                          |
| O-> pot_calc_energy_real                         43          43      0.09s  |
|    /                                                                        |
|   o-> pot_debug_check                            43          43      0.00s  |
|   o-> pot_calc_energy_really_real                43          43      0.09s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate_stress                     1                         |
|  /                                                                          |
| O-> hartree_calculate_stress                      1           1      0.09s  |
|    /                                                                        |
|   o-> hartree_check_inv_gsqr                      1           1      0.00s  |
|   o-> density_to_recip                            1           1      0.08s  |
|   o-> comms_reduce_gv_real                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_calc_energy_real                       43                         |
|  /                                                                          |
| O-> pot_calc_energy_really_real                  43          43      0.09s  |
|    /                                                                        |
|   o-> pot_debug_check                            43          43      0.00s  |
|   o-> comms_reduce_gv_real                       43          43      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_distribute_grids                     12                         |
|  /o-- density_write_parallel                     12                         |
|  /o-- wave_write_all_par                          8                         |
|  /o-- pot_write_parallel                          3                         |
| |/                                                                          |
| O-> comms_gather_gv_integer                      35          35      0.08s  |
+-----------------------------------------------------------------------------+
|   o-- ion_beta_recip_set                        600                         |
|  /o-- ion_dbetadcell_recip_ion                 3600                         |
| |/                                                                          |
| O-> basis_multiply_recip_reduced               4200        4200      0.08s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_apply_ks                       12                         |
|  /                                                                          |
| O-> wave_dot_wv_wv_ks                            12          12      0.08s  |
|    /                                                                        |
|   o-> local_dot_many                             12          12      0.08s  |
|   o-> comms_reduce_gv_complex                    12          12      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_dot_wv_wv_ks                          12                         |
|  /                                                                          |
| O-> local_dot_many                               12          12      0.08s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> ion_read                                      1           1      0.07s  |
|    /                                                                        |
|   o-> ion_atom_inquire_usp                        4           4      0.07s  |
+-----------------------------------------------------------------------------+
|   o-- ion_read                                    4                         |
|  /                                                                          |
| O-> ion_atom_inquire_usp                          4           4      0.07s  |
|    /                                                                        |
|   o-> comms_gcopy_integer                        20          20      0.00s  |
|   o-> comms_gcopy_real                            4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_find_occupancies                28                         |
|  /                                                                          |
| O-> electronic_permute_eigenstates               28          28      0.07s  |
|    /                                                                        |
|   o-> wave_allocate_slice                        56          56      0.00s  |
|   o-> wave_initialise_slice                      56          56      0.03s  |
|   o-> wave_copy_wv_slice                        726         726      0.02s  |
|   o-> wave_copy_slice_wv                        676         676      0.01s  |
|   o-> wave_copy_slice_slice                     676         676      0.01s  |
|   o-> wave_deallocate_slice                      56          56      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_dbetadcell_recip_ion                 3600                         |
|  /                                                                          |
| O-> ion_apply_ylmp                             3600        3600      0.07s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> cell_read_wrapped                             1           1      0.06s  |
|    /                                                                        |
|   o-> cell_setup_keywords                         1           1      0.00s  |
|   o-> comms_gcopy_integer                        63          63      0.00s  |
|   o-> cell_check_keywords                         1           1      0.00s  |
|   o-> cell_read_line_real                       324         324      0.01s  |
|   o-> cell_read_line_char                       108         108      0.00s  |
|   o-> cell_allocate                               1           1      0.00s  |
|   o-> cell_calculate_volume                       2           2      0.00s  |
|   o-> cell_recip_lattice                          1           1      0.00s  |
|   o-> cell_analyse_symmetry_wrapped               1           1      0.00s  |
|   o-> cell_detect_MP                              1           1      0.00s  |
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
|   o-> cell_find_related_atoms                     1           1      0.00s  |
|   o-> cell_detect_primitive                       1           1      0.00s  |
|   o-> cell_check_group                            1           1      0.00s  |
|   o-> cell_count_symmetry_translations_wrap       1           1      0.00s  |
|   o-> cell_check_group_equiv                      1           1      0.00s  |
|   o-> cell_symmetry_test                          2           2      0.00s  |
|   o-> cell_check_symmetry_ops_spin                1           1      0.00s  |
|   o-> cell_generate_ionic_constraints             1           1      0.00s  |
|   o-> comms_gcopy_logical                        11          11      0.00s  |
|   o-> comms_gcopy_real                           67          67      0.00s  |
|   o-> comms_gcopy_character                     116         116      0.00s  |
|   o-> cell_generate_qpoints_local                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           32                         |
|  /o-- xc_calculate_potential                     17                         |
| |/                                                                          |
| O-> pot_add                                      49          49      0.06s  |
|    /                                                                        |
|   o-> pot_debug_check                            98          98      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locps_calculate_potential                   8                         |
|  /o-- locps_calculate_forces                      4                         |
|  /o-- locps_calculate_stress                      4                         |
| |/                                                                          |
| O-> basis_scale_recip_grid                       16          16      0.06s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_d_real                     15                         |
|  /o-- nlpot_calculate_dddR                        1                         |
|  /o-- nlpot_calculate_dddcell_real                1                         |
| |/                                                                          |
| O-> basis_recip_fine_to_half_grid                17          17      0.06s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /                                                                          |
| O-> ewald_calculate_energy                        1           1      0.06s  |
|    /                                                                        |
|   o-> comms_reduce_kp_logical                     2           2      0.00s  |
|   o-> comms_reduce_bnd_logical                    2           2      0.00s  |
|   o-> comms_reduce_gv_logical                     2           2      0.00s  |
|   o-> ewald_calculate_num_cells                   1           1      0.01s  |
|   o-> comms_reduce_gv_real                        3           3      0.00s  |
|   o-> comms_reduce_bnd_real                       3           3      0.00s  |
|   o-> comms_reduce_kp_real                        3           3      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ewald_calculate_energy                      3                         |
|  /o-- nlpot_calculate_d_real                     15                         |
|  /o-- electronic_apply_H_energy_local             2                         |
|  /o-- hamiltonian_diagonalise_ks                 14                         |
|  /o-- electronic_find_fermi_free                651                         |
|  /o-- electronic_occupancy_update                42                         |
|  /o-- electronic_entropy_correction              14                         |
|  /o-- electronic_apply_H_energy_eigen            42                         |
|  /o-- density_calc_soft_wvfn_real                12                         |
|  /o-- scs_calculate                               2                         |
|  /o-- model_write_occ_eigenvalues                 8                         |
|  /o-- ewald_calculate_forces                      3                         |
|  /o-- nlpot_calculate_forces_r                    1                         |
|  /o-- ewald_calculate_stress                      3                         |
|  /o-- wave_kinetic_stress_wv                      1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
|  /o-- write_eigenvalues                           1                         |
| |/                                                                          |
| O-> comms_reduce_bnd_real                       815         815      0.05s  |
|    /                                                                        |
|   o-> comms_reduce_array_real                    48          48      0.05s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_mulliken                     1                         |
|  /                                                                          |
| O-> popn_calculate_charge_spilling                1           1      0.04s  |
|    /                                                                        |
|   o-> popn_nbands_occupied                        1           1      0.00s  |
|   o-> comms_reduce_kp_complex                     1           1      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_pseudo_scf                       194                         |
|  /                                                                          |
| O-> ion_atom_set_pseudo_H                       194         194      0.04s  |
|    /                                                                        |
|   o-> ion_atom_calc_pseudo_rho                  194         194      0.00s  |
|   o-> ion_atom_pseudo_hartree                   194         194      0.00s  |
|   o-> ion_atom_pseudo_xc                        194         194      0.02s  |
|   o-> ion_atom_regin                           2914        2914      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                 70                         |
|  /o-- hamiltonian_apply_ks                       24                         |
|  /o-- wave_rotate_wv_ks                          28                         |
|  /o-- wave_rotate_slice                         552                         |
|  /o-- electronic_permute_eigenstates             56                         |
| |/                                                                          |
| O-> wave_allocate_slice                         730         730      0.04s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              4                         |
|  /                                                                          |
| O-> ion_atom_read_usp                             4           4      0.04s  |
|    /                                                                        |
|   o-> ion_atom_derivative                         4           4      0.00s  |
|   o-> ion_atom_radin                             52          52      0.00s  |
|   o-> comms_gcopy_real                           72          72      0.00s  |
|   o-> comms_gcopy_integer                        28          28      0.00s  |
|   o-> comms_gcopy_logical                         8           8      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_create_lcao_basis_allkpt             300                         |
|  /                                                                          |
| O-> basis_radial_to_recip_reduced               300         300      0.04s  |
|    /                                                                        |
|   o-> basis_utils_interpolation                 300         300      0.03s  |
|   o-> basis_utils_apply_s_harmonic              108         108      0.00s  |
|   o-> basis_utils_apply_p_harmonic              192         192      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_populations              11556                         |
|  /                                                                          |
| O-> cell_verify_mixture_components            11556       11556      0.04s  |
+-----------------------------------------------------------------------------+
|   o-- wave_orthonormalise_over_wv_ks              1                         |
|  /o-- algor_diagonalise_hermitian_lapack        318                         |
| |/                                                                          |
| O-> comms_copy_gv_complex                       319         319      0.04s  |
+-----------------------------------------------------------------------------+
|   o-- cell_check_group                            4                         |
|  /o-- cell_symmetry_symbol_wrapped_wrapped        2                         |
|  /o-- cell_supercell                              2                         |
|  /o-- cell_generate_supercell_origins             2                         |
|  /o-- cell_set_supercell_symmetry                 2                         |
|  /o-- cell_find_reduced_cell                     64                         |
|  /o-- cell_sort_kpoints_with_recip              128                         |
|  /o-- scs_calculate                               4                         |
| |/                                                                          |
| O-> algor_invert_real                           208         208      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- ion_augment_charge_nospin_kp               12                         |
|  /o-- popn_calculate_charge_spilling              1                         |
|  /o-- popn_calculate_populations                  1                         |
| |/                                                                          |
| O-> comms_reduce_kp_complex                      14          14      0.03s  |
|    /                                                                        |
|   o-> comms_reduce_array_complex                 13          13      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_assign_grid_coordinates                 1           1      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- model_write_all                             8                         |
|  /                                                                          |
| O-> cell_dump                                     8           8      0.03s  |
|    /                                                                        |
|   o-> cell_dump_cell                              8           8      0.01s  |
|   o-> cell_dump_global                            8           8      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- wave_orthonormalise_over_wv_ks              1                         |
|  /o-- wave_orthonormalise_over_slice            276                         |
|  /o-- popn_calculate_matrices                     1                         |
| |/                                                                          |
| O-> algor_invert_complex                        278         278      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- ion_beta_recip_set                        108                         |
|  /o-- ion_dbetadcell_recip_ion                  108                         |
| |/                                                                          |
| O-> ion_cc_structure_factor                     216         216      0.03s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> electronic_apply_H_energy_local               1           1      0.02s  |
|    /                                                                        |
|   o-> wave_kinetic_eigenvalues_wv_ks              1           1      0.02s  |
|   o-> nlpot_calc_eigenvalues_nkns                 1           1      0.00s  |
|   o-> comms_reduce_bnd_real                       2           2      0.00s  |
|   o-> comms_reduce_kp_real                        2           2      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_mulliken                     1                         |
|  /                                                                          |
| O-> popn_calculate_density_matrix                 1           1      0.02s  |
|    /                                                                        |
|   o-> popn_nbands_occupied                        1           1      0.00s  |
|   o-> popn_invert_complex                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_distribute_grids                     12                         |
|  /o-- model_write_occ_eigenvalues                 4                         |
|  /o-- wave_write_all_par                          2                         |
|  /o-- write_eigenvalues                           1                         |
| |/                                                                          |
| O-> comms_gather_kp_integer                      19          19      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_pseudo_scf                         8                         |
|  /                                                                          |
| O-> ion_atom_init_pseudo_basis                    8           8      0.02s  |
|    /                                                                        |
|   o-> ion_atom_find_root                       1584        1584      0.01s  |
|   o-> ion_atom_regin                            784         784      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    15                         |
|  /                                                                          |
| O-> electronic_write_scf_energies                15          15      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- model_write_all                             4                         |
|  /                                                                          |
| O-> parameters_dump                               4           4      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /o-- cell_read_wrapped                          67                         |
|  /o-- cell_generate_qpoints_local                 4                         |
|  /o-- ion_atom_inquire_usp                        4                         |
|  /o-- parameters_bcast                           94                         |
|  /o-- ion_atom_read_usp                          72                         |
|  /o-- ewald_calculate_num_cells                 110                         |
|  /o-- electronic_find_fermi_free                623                         |
|  /o-- hirshfeld_calculate                         4                         |
|  /o-- mbd2_calculate                             16                         |
| |/                                                                          |
| O-> comms_gcopy_real                            995         995      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> parameters_read                               1           1      0.02s  |
|    /                                                                        |
|   o-> parameters_keywords_setup                   1           1      0.00s  |
|   o-> comms_gcopy_integer                         3           3      0.00s  |
|   o-> parameters_read_xc_block                    1           1      0.00s  |
|   o-> parameters_validate                         1           1      0.00s  |
|   o-> parameters_bcast                            1           1      0.00s  |
|   o-> algor_set_random_seed                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_write_all_par                        466                         |
|  /                                                                          |
| O-> comms_gather_gv_complex                     466         466      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_H                      8                         |
|  /o-- ion_atom_set_pseudo_H                     194                         |
| |/                                                                          |
| O-> ion_atom_pseudo_xc                          202         202      0.02s  |
|    /                                                                        |
|   o-> ion_atom_derivative                       404         404      0.00s  |
|   o-> ion_atom_regin                            202         202      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_dump                                   8                         |
|  /                                                                          |
| O-> cell_dump_global                              8           8      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- ion_set_Q_at_origin_recip                 138                         |
|  /                                                                          |
| O-> ion_generate_QLnm                           138         138      0.02s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_pulay                        1                         |
|  /                                                                          |
| O-> dm_initialise                                 1           1      0.02s  |
|    /                                                                        |
|   o-> dm_mix_density_initialise                   1           1      0.00s  |
|   o-> dm_mix_density_allocate                    46          46      0.00s  |
|   o-> dm_mix_density_zero                        46          46      0.00s  |
|   o-> density_allocate                            1           1      0.00s  |
|   o-> density_zero                                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              8                         |
|  /o-- ion_generate_orbitals                       4                         |
| |/                                                                          |
| O-> ion_atom_allocate_pspot                      12          12      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> wave_kinetic_stress_wv                        1           1      0.01s  |
|    /                                                                        |
|   o-> comms_reduce_gv_real                        1           1      0.00s  |
|   o-> comms_reduce_bnd_real                       1           1      0.00s  |
|   o-> comms_reduce_kp_real                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_kerker                       1                         |
|  /o-- dm_mix_density_pulay                      363                         |
| |/                                                                          |
| O-> dm_mix_density_dot                          364         364      0.01s  |
|    /                                                                        |
|   o-> comms_reduce_gv_complex                   364         364      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_atom                4152                         |
|  /                                                                          |
| O-> ion_atom_interpolate                       4152        4152      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_atom                4152                         |
|  /                                                                          |
| O-> ion_atom_locate                            4152        4152      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      5                         |
|  /o-- dftd_sedc_initialize                        1                         |
|  /o-- cell_kpoints_spacing_wrapped               16                         |
|  /o-- cell_generate_MP_set_wrapped               16                         |
| |/                                                                          |
| O-> bib_add                                      38          38      0.01s  |
|    /                                                                        |
|   o-> bib_setup                                   1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> parameters_output                             1           1      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- density_write_parallel                      4                         |
|  /o-- pot_write_parallel                          1                         |
| |/                                                                          |
| O-> comms_gather_gv_real                          5           5      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- wave_beta_phi_wv_ks                        17                         |
|  /o-- wave_beta_phi_slice                      2916                         |
| |/                                                                          |
| O-> wave_calc_ps_q_nonzero                     2933        2933      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- hirshfeld_calculate                         2                         |
|  /                                                                          |
| O-> cell_supercell                                2           2      0.01s  |
|    /                                                                        |
|   o-> cell_deallocate                             2           2      0.00s  |
|   o-> cell_num_supercells                         2           2      0.00s  |
|   o-> algor_invert_real                           2           2      0.00s  |
|   o-> cell_allocate                               2           2      0.00s  |
|   o-> cell_recip_lattice                          2           2      0.00s  |
|   o-> cell_calculate_volume                       2           2      0.00s  |
|   o-> cell_generate_supercell_origins             2           2      0.00s  |
|   o-> cell_set_supercell_symmetry                 2           2      0.00s  |
|   o-> cell_copy_kpoints                           2           2      0.00s  |
|   o-> cell_kpoints_spacing_wrapped                2           2      0.00s  |
|   o-> cell_generate_ionic_constraints             2           2      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              4                         |
|  /                                                                          |
| O-> ion_set_data                                  4           4      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /o-- cell_supercell                              2                         |
| |/                                                                          |
| O-> cell_generate_ionic_constraints               3           3      0.01s  |
|    /                                                                        |
|   o-> algor_uniform_random                     1872        1872      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_dump                                   8                         |
|  /                                                                          |
| O-> cell_dump_cell                                8           8      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           16                         |
|  /                                                                          |
| O-> pot_complex_to_real                          16          16      0.01s  |
|    /                                                                        |
|   o-> pot_deallocate_cmplx                        3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hartree_calculate_potential                16                         |
|  /o-- hartree_calculate_stress                    1                         |
| |/                                                                          |
| O-> hartree_check_inv_gsqr                       17          17      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- xc_calculate_potential                     17                         |
|  /                                                                          |
| O-> pot_copy                                     17          17      0.01s  |
|    /                                                                        |
|   o-> pot_debug_check                            34          34      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_kerker                       2                         |
|  /o-- dm_mix_density_pulay                      176                         |
| |/                                                                          |
| O-> dm_mix_density_add                          178         178      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks               1052                         |
|  /o-- wave_Sorthonormalise_slice                276                         |
| |/                                                                          |
| O-> comms_reduce_bnd_integer                   1328        1328      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- cell_supercell                              2                         |
|  /o-- mbd2_calculate                             14                         |
| |/                                                                          |
| O-> cell_kpoints_spacing_wrapped                 16          16      0.01s  |
|    /                                                                        |
|   o-> bib_add                                    16          16      0.00s  |
|   o-> cell_generate_MP_set_wrapped               16          16      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_initialise_cell                     1                         |
|  /o-- electronic_initialise                       1                         |
|  /o-- dm_initialise                               1                         |
| |/                                                                          |
| O-> density_zero                                  3           3      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- check_elec_ground_state                     1                         |
|  /                                                                          |
| O-> castep_calc_storage                           1           1      0.01s  |
|    /                                                                        |
|   o-> cell_calc_storage_global_wrapped            1           1      0.00s  |
|   o-> basis_calc_storage                          1           1      0.00s  |
|   o-> ion_calc_storage                            1           1      0.00s  |
|   o-> model_calc_storage                          1           1      0.00s  |
|   o-> ewald_calc_storage                          1           1      0.01s  |
|   o-> density_calc_storage                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> comms_parallel_strategy                       1           1      0.01s  |
|    /                                                                        |
|   o-> find_strategy                               1           1      0.00s  |
|   o-> assign_nodes                                1           1      0.00s  |
|   o-> reassign_nodes                              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locps_calculate_potential                   2                         |
|  /                                                                          |
| O-> basis_scale_real_grid                         2           2      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- check_forces_stresses                       1                         |
|  /                                                                          |
| O-> firstd_output_forces                          1           1      0.01s  |
|    /                                                                        |
|   o-> cell_format_atom_label                    108         108      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep_calc_storage                         1                         |
|  /                                                                          |
| O-> ewald_calc_storage                            1           1      0.01s  |
|    /                                                                        |
|   o-> ewald_calculate_num_cells                   1           1      0.01s  |
|   o-> algor_sizeof_real                          12          12      0.00s  |
|   o-> algor_sizeof_cmplx                          9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_symmetrise_gv_parallel            100                         |
|  /                                                                          |
| O-> comms_gather_all_gv_integer                 100         100      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> nlxc_initialise                               1           1      0.01s  |
|    /                                                                        |
|   o-> density_allocate                            1           1      0.00s  |
|   o-> density_copy                                1           1      0.00s  |
|   o-> cell_generate_qpoints_local                 1           1      0.00s  |
|   o-> density_deallocate                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                         324                         |
|  /                                                                          |
| O-> cell_read_line_real                         324         324      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- cell_generate_cell_constraints              1                         |
|  /o-- wave_initialise_wv                        233                         |
| |/                                                                          |
| O-> algor_uniform_random_array                  234         234      0.01s  |
|    /                                                                        |
|   o-> algor_set_random_seed                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_prepare_precon_ks                    14                         |
|  /o-- hamiltonian_diagonalise_ks                828                         |
| |/                                                                          |
| O-> comms_copy_gv_logical                       842         842      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_output_forces                      108                         |
|  /o-- popn_calculate_populations               1710                         |
| |/                                                                          |
| O-> cell_format_atom_label                     1818        1818      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              4                         |
|  /o-- ion_generate_orbitals                       4                         |
| |/                                                                          |
| O-> ion_set_psp                                   8           8      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- <parent(s) not traced>                      1                         |
|  /                                                                          |
| O-> basis_deallocate                              1           1      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_basis               1584                         |
|  /                                                                          |
| O-> ion_atom_find_root                         1584        1584      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- basis_radial_to_recip_reduced             192                         |
|  /                                                                          |
| O-> basis_utils_apply_p_harmonic                192         192      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_ps_diag                          890                         |
|  /                                                                          |
| O-> ion_atom_rectoreal                          890         890      0.01s  |
+-----------------------------------------------------------------------------+
|   o-- cell_generate_ionic_constraints          1872                         |
|  /                                                                          |
| O-> algor_uniform_random                       1872        1872      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_apply_slice                           984                         |
|  /o-- pot_nongamma_apply_slice                  984                         |
|  /o-- nlpot_apply_add_slice_r                   984                         |
| |/                                                                          |
| O-> wave_spin_type_slice                       2952        2952      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> hirshfeld_output                              1           1      0.00s  |
|    /                                                                        |
|   o-> parameters_nspins                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_kpoints_spacing_wrapped               16                         |
|  /                                                                          |
| O-> cell_generate_MP_set_wrapped                 16          16      0.00s  |
|    /                                                                        |
|   o-> bib_add                                    16          16      0.00s  |
|   o-> cell_reduce_kpoints_internal               16          16      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_initialise                               1                         |
|  /                                                                          |
| O-> dm_mix_density_initialise                     1           1      0.00s  |
|    /                                                                        |
|   o-> dm_count_plane_waves                        1           1      0.00s  |
|   o-> dm_assign_plane_wave_indices                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_zero                                   54                         |
|  /o-- pot_add                                    98                         |
|  /o-- pot_calc_energy_real                       43                         |
|  /o-- pot_calc_energy_really_real                43                         |
|  /o-- pot_copy                                   34                         |
|  /o-- pot_apply_slice                           984                         |
|  /o-- pot_interpolate                            12                         |
| |/                                                                          |
| O-> pot_debug_check                            1268        1268      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> cell_output_wrapped                           1           1      0.00s  |
|    /                                                                        |
|   o-> cell_cart_lattice_to_abc                    1           1      0.00s  |
|   o-> cell_check_group                            1           1      0.00s  |
|   o-> cell_factor_group_symmetry_wrapped          1           1      0.00s  |
|   o-> cell_symmetry_symbol_wrapped_wrapped        2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_kerker                       1                         |
|  /o-- dm_mix_density_pulay                       55                         |
| |/                                                                          |
| O-> dm_mix_density_copy                          56          56      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_set_pseudo_H                     194                         |
|  /                                                                          |
| O-> ion_atom_calc_pseudo_rho                    194         194      0.00s  |
|    /                                                                        |
|   o-> ion_atom_regin                            194         194      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- parameters_read                             1                         |
|  /                                                                          |
| O-> parameters_bcast                              1           1      0.00s  |
|    /                                                                        |
|   o-> comms_gcopy_logical                        64          64      0.00s  |
|   o-> comms_gcopy_character                      93          93      0.00s  |
|   o-> comms_gcopy_integer                        77          77      0.00s  |
|   o-> parameters_reallocate_xc                    7           7      0.00s  |
|   o-> comms_gcopy_real                           94          94      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> dftd_sedc_initialize                          1           1      0.00s  |
|    /                                                                        |
|   o-> dftd_new_calc                               1           1      0.00s  |
|   o-> sedc_get_default_params                     1           1      0.00s  |
|   o-> dftd_read_custom_params                     1           1      0.00s  |
|   o-> dftd_check_params                           1           1      0.00s  |
|   o-> bib_add                                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                          63                         |
|  /o-- cell_generate_qpoints_local                 4                         |
|  /o-- ion_atom_inquire_usp                       20                         |
|  /o-- parameters_read                             3                         |
|  /o-- parameters_bcast                           77                         |
|  /o-- parameters_reallocate_xc                    7                         |
|  /o-- ion_atom_read_usp                          28                         |
|  /o-- ewald_calculate_num_cells                 132                         |
|  /o-- hirshfeld_calculate                         2                         |
|  /o-- wave_write_all_par                          2                         |
| |/                                                                          |
| O-> comms_gcopy_integer                         338         338      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /o-- locpot_calculate                           16                         |
|  /o-- locpot_check_locps_core                     1                         |
|  /o-- xc_calculate_potential                     17                         |
|  /o-- locpot_calculate_forces                     1                         |
|  /o-- firstd_calculate_forces                     1                         |
|  /o-- locpot_calculate_stress                     1                         |
|  /o-- locps_calculate_stress                      1                         |
|  /o-- xc_calculate_stress_density                 1                         |
|  /o-- firstd_calculate_stress                     1                         |
|  /o-- write_local_pot_and_den                     1                         |
| |/                                                                          |
| O-> pot_allocate                                 42          42      0.00s  |
|    /                                                                        |
|   o-> pot_deallocate_real                         3           3      0.00s  |
|   o-> pot_deallocate_nc                          42          42      0.00s  |
|   o-> pot_deallocate_cmplx                       39          39      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlxc_initialise                             1                         |
|  /o-- density_write                               4                         |
| |/                                                                          |
| O-> density_copy                                  5           5      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_find_fermi_free               3262                         |
|  /                                                                          |
| O-> algor_broadening                           3262        3262      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                          11                         |
|  /o-- cell_generate_qpoints_local                 2                         |
|  /o-- parameters_bcast                           64                         |
|  /o-- ion_atom_read_usp                           8                         |
|  /o-- ewald_calculate_num_cells                  66                         |
|  /o-- electronic_minimisation                     1                         |
|  /o-- hubbard_initialise                          3                         |
|  /o-- electronic_find_fermi_free                637                         |
|  /o-- hirshfeld_calculate                         4                         |
| |/                                                                          |
| O-> comms_gcopy_logical                         796         796      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_write_all                             4                         |
|  /                                                                          |
| O-> model_write_occ_eigenvalues                   4           4      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_bnd_real                       8           8      0.00s  |
|   o-> comms_gather_kp_integer                     4           4      0.00s  |
|   o-> comms_recv_real                            96          96      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_calc_storage                             1                         |
|  /o-- dm_mix_density_initialise                   1                         |
| |/                                                                          |
| O-> dm_count_plane_waves                          2           2      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_gv_integer                     6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_initialise                              46                         |
|  /                                                                          |
| O-> dm_mix_density_zero                          46          46      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_density_matrix               1                         |
|  /                                                                          |
| O-> popn_invert_complex                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_initialise                   1                         |
|  /                                                                          |
| O-> dm_assign_plane_wave_indices                  1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- mbd2_calculate                             14                         |
|  /o-- comms_restore_strategy                     14                         |
| |/                                                                          |
| O-> comms_reassign_strategy                      28          28      0.00s  |
|    /                                                                        |
|   o-> find_strategy                              14          14      0.00s  |
|   o-> assign_nodes                               28          28      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_write_all_par                         64                         |
|  /                                                                          |
| O-> comms_recv_integer                           64          64      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hamiltonian_diagonalise_ks                 14                         |
|  /                                                                          |
| O-> wave_calc_precon                             14          14      0.00s  |
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
|   o-- electronic_find_fermi_free                665                         |
|  /                                                                          |
| O-> parameters_band_degeneracy                  665         665      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                         108                         |
|  /                                                                          |
| O-> cell_read_line_char                         108         108      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- check_forces_stresses                       1                         |
|  /                                                                          |
| O-> firstd_output_stress                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /o-- model_initialise                            2                         |
|  /o-- hirshfeld_calculate                         8                         |
|  /o-- mbd2_fd                                    14                         |
|  /o-- mbd2_calculate                             14                         |
| |/                                                                          |
| O-> cell_copy                                    39          39      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_analyse_symmetry_wrapped                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_map_fine_recip_half_full                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_read_usp                           4                         |
|  /o-- ion_atom_pseudo_xc                        404                         |
| |/                                                                          |
| O-> ion_atom_derivative                         408         408      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_H                      8                         |
|  /o-- ion_atom_set_pseudo_H                     194                         |
| |/                                                                          |
| O-> ion_atom_pseudo_hartree                     202         202      0.00s  |
|    /                                                                        |
|   o-> ion_atom_regin                            202         202      0.00s  |
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
|   o-- dftd_sedc_initialize                        1                         |
|  /o-- sedc_calculate_properties                   3                         |
|  /o-- mbd2_fd                                    12                         |
| |/                                                                          |
| O-> dftd_new_calc                                16          16      0.00s  |
|    /                                                                        |
|   o-> get_xc_functional_name                     16          16      0.00s  |
|   o-> cell_pack_vec_real                         32          32      0.00s  |
|   o-> cell_frac_to_cart_cell                     16          16      0.00s  |
|   o-> generate_periodic_table                    16          16      0.00s  |
|   o-> dftd_atomic_number                         64          64      0.00s  |
|   o-> dftd_scheme_to_int                         16          16      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_pseudo_scf                         8                         |
|  /                                                                          |
| O-> ion_atom_init_pseudo_H                        8           8      0.00s  |
|    /                                                                        |
|   o-> ion_atom_regin                            112         112      0.00s  |
|   o-> ion_atom_pseudo_hartree                     8           8      0.00s  |
|   o-> ion_atom_pseudo_xc                          8           8      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_generate_MP_set_wrapped               16                         |
|  /o-- cell_unfold_kpoints_arg_trace              14                         |
| |/                                                                          |
| O-> cell_reduce_kpoints_internal                 30          30      0.00s  |
|    /                                                                        |
|   o-> cell_sort_kpoints_with_recip               64          64      0.00s  |
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
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> electronic_scf_footer                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                         116                         |
|  /o-- parameters_bcast                           93                         |
|  /o-- ewald_calculate_num_cells                  22                         |
|  /o-- electronic_minimisation                     1                         |
| |/                                                                          |
| O-> comms_gcopy_character                       232         232      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- mbd2_calculate                             14                         |
|  /                                                                          |
| O-> comms_restore_strategy                       14          14      0.00s  |
|    /                                                                        |
|   o-> comms_reassign_strategy                    14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_reduce_kpoints_internal               64                         |
|  /                                                                          |
| O-> cell_sort_kpoints_with_recip                 64          64      0.00s  |
|    /                                                                        |
|   o-> cell_find_reduced_cell                     64          64      0.00s  |
|   o-> algor_invert_real                         128         128      0.00s  |
|   o-> algor_sort                                 64          64      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_reset                                 2                         |
|  /o-- nlxc_initialise                             1                         |
|  /o-- dm_finalise                                 1                         |
|  /o-- electronic_finalise                         1                         |
|  /o-- density_write                               4                         |
| |/                                                                          |
| O-> density_deallocate                            9           9      0.00s  |
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
|   o-- algor_diagonalise_hermitian_lapack        318                         |
|  /o-- electronic_order_eigenvalues               14                         |
| |/                                                                          |
| O-> comms_copy_gv_real                          332         332      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- <parent(s) not traced>                      1                         |
|  /                                                                          |
| O-> ion_finalise                                  1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_map_standard_to_fine                    1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_allocate_wv                            3                         |
|  /o-- wave_setup_slice                          514                         |
| |/                                                                          |
| O-> wave_band_basis_initialise                  517         517      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              9                         |
|  /o-- ion_generate_orbitals                       4                         |
| |/                                                                          |
| O-> ion_atom_deallocate_pspot                    13          13      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_phonon_fine_data                    1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_dbetadR_recip_ion                     600                         |
|  /                                                                          |
| O-> ion_proj_index                              600         600      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    15                         |
|  /                                                                          |
| O-> electronic_store_energy                      15          15      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_gv_logical                    30          30      0.00s  |
|   o-> comms_reduce_bnd_logical                   30          30      0.00s  |
|   o-> comms_reduce_kp_logical                    30          30      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_assign_pw_gvectors                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_write_occ_eigenvalues                96                         |
|  /o-- wave_write_all_par                         16                         |
|  /o-- write_eigenvalues                          24                         |
| |/                                                                          |
| O-> comms_recv_real                             136         136      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_mix_density_kerker                       1                         |
|  /o-- dm_mix_density_pulay                       11                         |
| |/                                                                          |
| O-> dm_apply_kerker                              12          12      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_assign_plane_wave_indexes               1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_find_occupancies                14                         |
|  /                                                                          |
| O-> electronic_order_eigenvalues                 14          14      0.00s  |
|    /                                                                        |
|   o-> algor_sort                                 14          14      0.00s  |
|   o-> comms_copy_gv_integer                      14          14      0.00s  |
|   o-> comms_copy_gv_real                         14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- mbd2_calculate                             14                         |
|  /                                                                          |
| O-> cell_unfold_kpoints_inplace_wrapped          14          14      0.00s  |
|    /                                                                        |
|   o-> cell_unfold_kpoints_arg_trace              14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_prepare_init_wvfn                      1                         |
|  /o-- ewald_calculate_energy                      2                         |
|  /o-- electronic_store_energy                    30                         |
|  /o-- electronic_check_occupancies               14                         |
|  /o-- sedc_calculate_properties                   3                         |
|  /o-- ewald_calculate_forces                      2                         |
|  /o-- ewald_calculate_stress                      2                         |
|  /o-- bib_output                                  1                         |
| |/                                                                          |
| O-> comms_reduce_kp_logical                      55          55      0.00s  |
|    /                                                                        |
|   o-> comms_lcopy                                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_unfold_kpoints_inplace_wrapped        14                         |
|  /                                                                          |
| O-> cell_unfold_kpoints_arg_trace                14          14      0.00s  |
|    /                                                                        |
|   o-> cell_reduce_kpoints_internal               14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calculate                           16                         |
|  /o-- electronic_finalise                         1                         |
|  /o-- locpot_calculate_forces                     1                         |
|  /o-- firstd_calculate_forces                     1                         |
|  /o-- locps_calculate_stress                      1                         |
|  /o-- xc_calculate_stress_density                 1                         |
|  /o-- locpot_calculate_stress                     1                         |
|  /o-- firstd_calculate_stress                     1                         |
|  /o-- write_local_pot_and_den                     1                         |
| |/                                                                          |
| O-> pot_deallocate                               24          24      0.00s  |
|    /                                                                        |
|   o-> pot_deallocate_real                        24          24      0.00s  |
|   o-> pot_deallocate_cmplx                       24          24      0.00s  |
|   o-> pot_deallocate_nc                          24          24      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              1                         |
|  /                                                                          |
| O-> ion_clebsch_gordan                            1           1      0.00s  |
|    /                                                                        |
|   o-> init_factorial                              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ewald_calculate_energy                      2                         |
|  /o-- electronic_store_energy                    30                         |
|  /o-- hamiltonian_diagonalise_ks                290                         |
|  /o-- electronic_check_occupancies               14                         |
|  /o-- sedc_calculate_properties                   3                         |
|  /o-- ewald_calculate_forces                      2                         |
|  /o-- ewald_calculate_stress                      2                         |
|  /o-- bib_output                                  1                         |
| |/                                                                          |
| O-> comms_reduce_bnd_logical                    344         344      0.00s  |
|    /                                                                        |
|   o-> comms_lcopy                                 1           1      0.00s  |
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
|   o-- parameters_read                             1                         |
|  /                                                                          |
| O-> parameters_keywords_setup                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_phonon_data                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    14                         |
|  /                                                                          |
| O-> electronic_check_occupancies                 14          14      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_bnd_logical                   14          14      0.00s  |
|   o-> comms_reduce_gv_logical                    14          14      0.00s  |
|   o-> comms_reduce_kp_logical                    14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_occupancy_update                14                         |
|  /                                                                          |
| O-> electronic_entropy_correction                14          14      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_bnd_real                      14          14      0.00s  |
|   o-> comms_reduce_kp_real                       14          14      0.00s  |
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
|  /o-- density_symmetrise_calc_storage             3                         |
|  /o-- locps_calc_storage                          3                         |
| |/                                                                          |
| O-> algor_sizeof_cmplx                           59          59      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_real                          59          59      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- check_elec_ground_state                     1                         |
|  /                                                                          |
| O-> castep_report_storage                         1           1      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_gv_real                        2           2      0.00s  |
|   o-> comms_reduce_kp_real                        2           2      0.00s  |
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
|   o-- pot_allocate                                3                         |
|  /o-- pot_deallocate                             24                         |
| |/                                                                          |
| O-> pot_deallocate_real                          27          27      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_populations                  1                         |
|  /                                                                          |
| O-> popn_sort                                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_order_eigenvalues               14                         |
|  /o-- cell_sort_kpoints_with_recip               64                         |
| |/                                                                          |
| O-> algor_sort                                   78          78      0.00s  |
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
|   o-- cell_sort_kpoints_with_recip               64                         |
|  /                                                                          |
| O-> cell_find_reduced_cell                       64          64      0.00s  |
|    /                                                                        |
|   o-> algor_invert_real                          64          64      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_calc_storage_cell                     20                         |
|  /o-- cell_calc_storage_global_wrapped           34                         |
|  /o-- basis_calc_storage                          9                         |
|  /o-- algor_sizeof_cmplx                         59                         |
|  /o-- ion_calc_storage                           24                         |
|  /o-- model_calc_storage                          2                         |
|  /o-- ewald_calc_storage                         12                         |
|  /o-- phonon_store_mem_estimate                   2                         |
|  /o-- dm_calc_storage                             3                         |
|  /o-- electronic_calc_storage                     2                         |
|  /o-- hamiltonian_calc_storage                    6                         |
|  /o-- hartree_calc_storage                        3                         |
| |/                                                                          |
| O-> algor_sizeof_real                           176         176      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- parameters_read                             1                         |
|  /                                                                          |
| O-> parameters_read_xc_block                      1           1      0.00s  |
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
|   o-- check_elec_ground_state                     2                         |
|  /                                                                          |
| O-> firstd_calc_storage                           2           2      0.00s  |
|    /                                                                        |
|   o-> pot_calc_storage                            2           2      0.00s  |
|   o-> locpot_calc_storage                         2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_new_calc                              16                         |
|  /                                                                          |
| O-> generate_periodic_table                      16          16      0.00s  |
|    /                                                                        |
|   o-> dftd_atomic_number                         64          64      0.00s  |
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
|  /o-- density_symmetrise_calc_storage            11                         |
| |/                                                                          |
| O-> algor_sizeof_int                            139         139      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_setup_keywords                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_initialise                              46                         |
|  /                                                                          |
| O-> dm_mix_density_allocate                      46          46      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> dftd_sedc_print_corr_energies                 1           1      0.00s  |
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
|   o-- firstd_calculate_forces                     1                         |
|  /                                                                          |
| O-> sedc_forces                                   1           1      0.00s  |
|    /                                                                        |
|   o-> sedc_calculate_properties                   1           1      0.00s  |
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
|   o-- generate_periodic_table                    64                         |
|  /o-- dftd_new_calc                              64                         |
| |/                                                                          |
| O-> dftd_atomic_number                          128         128      0.00s  |
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
|   o-- electronic_minimisation                    15                         |
|  /                                                                          |
| O-> electronic_write_energies                    15          15      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- bib_add                                     1                         |
|  /                                                                          |
| O-> bib_setup                                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_new_calc                              16                         |
|  /                                                                          |
| O-> dftd_scheme_to_int                           16          16      0.00s  |
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
|   o-- electronic_calc_storage                     1                         |
|  /                                                                          |
| O-> density_symmetrise_calc_storage               1           1      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_cmplx                          3           3      0.00s  |
|   o-> algor_sizeof_int                           11          11      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_allocate                              13                         |
|  /o-- model_reset                                 8                         |
|  /o-- cell_supercell                              2                         |
|  /o-- hirshfeld_calculate                         6                         |
| |/                                                                          |
| O-> cell_deallocate                              29          29      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_atom                   8                         |
|  /                                                                          |
| O-> ion_atom_set_pseudo_occ                       8           8      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- mbd2_calculate                             28                         |
|  /                                                                          |
| O-> cell_set_current_kpoints_wrapped             28          28      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> electronic_scf_banner                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep_calc_storage                         1                         |
|  /                                                                          |
| O-> cell_calc_storage_global_wrapped              1           1      0.00s  |
|    /                                                                        |
|   o-> cell_calc_storage_cell                      1           1      0.00s  |
|   o-> algor_sizeof_real                          34          34      0.00s  |
|   o-> algor_sizeof_int                            7           7      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            1                         |
|  /                                                                          |
| O-> basis_calculate_cut_off                       1           1      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_gv_real                        1           1      0.00s  |
|   o-> comms_reduce_kp_real                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_distribute_grids                      1                         |
|  /                                                                          |
| O-> basis_utils_sort_columns                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- parameters_read                             1                         |
|  /                                                                          |
| O-> parameters_validate                           1           1      0.00s  |
|    /                                                                        |
|   o-> parameters_xc_allowed                       7           7      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_count_plane_waves                     3                         |
|  /o-- electronic_analyse_occ_all                  2                         |
| |/                                                                          |
| O-> comms_reduce_kp_integer                       5           5      0.00s  |
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
|   o-- pot_complex_to_real                         3                         |
|  /o-- pot_allocate                               39                         |
|  /o-- pot_deallocate                             24                         |
| |/                                                                          |
| O-> pot_deallocate_cmplx                         66          66      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_symmetrise_gv_parallel             25                         |
|  /                                                                          |
| O-> density_set_grid_lookup                      25          25      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_initialise                            1                         |
|  /o-- electronic_initialise                       1                         |
|  /o-- popn_create_lcao_basis_allkpt               1                         |
| |/                                                                          |
| O-> wave_allocate_wv                              3           3      0.00s  |
|    /                                                                        |
|   o-> wave_band_basis_initialise                  3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_write_all_par                         32                         |
|  /                                                                          |
| O-> comms_send_integer                           32          32      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_read_usp                          52                         |
|  /                                                                          |
| O-> ion_atom_radin                               52          52      0.00s  |
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
|   o-- electronic_initialise                       1                         |
|  /o-- nlpot_calculate_forces_r                    1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
| |/                                                                          |
| O-> nlpot_allocate_nl_d                           3           3      0.00s  |
|    /                                                                        |
|   o-> nlpot_deallocate_nl_d                       3           3      0.00s  |
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
|   o-- castep_calc_storage                         1                         |
|  /                                                                          |
| O-> ion_calc_storage                              1           1      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_int                           14          14      0.00s  |
|   o-> algor_sizeof_cmplx                          3           3      0.00s  |
|   o-> algor_sizeof_real                          24          24      0.00s  |
|   o-> algor_sizeof_logical                        2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- wave_initialise_wv                          1                         |
|  /                                                                          |
| O-> wave_prepare_init_wvfn                        1           1      0.00s  |
|    /                                                                        |
|   o-> comms_reduce_gv_integer                     2           2      0.00s  |
|   o-> comms_reduce_kp_logical                     1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_symmetrise_gv_parallel             25                         |
|  /                                                                          |
| O-> comms_copy_bnd_complex                       25          25      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- parameters_bcast                            7                         |
|  /                                                                          |
| O-> parameters_reallocate_xc                      7           7      0.00s  |
|    /                                                                        |
|   o-> comms_gcopy_integer                         7           7      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- pot_allocate                               42                         |
|  /o-- pot_deallocate                             24                         |
| |/                                                                          |
| O-> pot_deallocate_nc                            66          66      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- check_elec_ground_state                     3                         |
|  /                                                                          |
| O-> phonon_store_mem_estimate                     3           3      0.00s  |
|    /                                                                        |
|   o-> secondd_store_mem_estimate                  3           3      0.00s  |
|   o-> algor_sizeof_cmplx                          1           1      0.00s  |
|   o-> algor_sizeof_real                           2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /o-- hirshfeld_calculate                         2                         |
|  /o-- mbd2_calculate                             28                         |
| |/                                                                          |
| O-> cell_distribute_kpoints_wrapped              31          31      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_calc_storage                     1                         |
|  /                                                                          |
| O-> wave_calc_storage_bnd                         1           1      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_cmplx                          2           2      0.00s  |
|   o-> ion_set_projectors                          1           1      0.00s  |
|   o-> algor_sizeof_logical                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep_calc_storage                         1                         |
|  /                                                                          |
| O-> basis_calc_storage                            1           1      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_int                           25          25      0.00s  |
|   o-> algor_sizeof_real                           9           9      0.00s  |
|   o-> algor_sizeof_logical                        3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- comms_parallel_strategy                     1                         |
|  /o-- comms_reassign_strategy                    14                         |
| |/                                                                          |
| O-> find_strategy                                15          15      0.00s  |
|    /                                                                        |
|   o-> best_mixed_strategy                         1           1      0.00s  |
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
|   o-- wave_initialise_wv                          1                         |
|  /o-- model_init_occ_eigenvalues_ef               1                         |
|  /o-- density_calculate_soft_wvfn                24                         |
|  /o-- density_augment_complex                    12                         |
|  /o-- nlpot_calculate_forces                      1                         |
|  /o-- nlpot_calculate_forces_r                    1                         |
|  /o-- nlpot_calculate_stress                      1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
| |/                                                                          |
| O-> wave_spin_type_wv                            42          42      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_find_related_atoms                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_set_supercell_symmetry                 2                         |
|  /                                                                          |
| O-> cell_find_related_atoms_supercell             2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_initialise                              1                         |
|  /                                                                          |
| O-> ion_allocate                                  1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- comms_parallel_strategy                     1                         |
|  /o-- comms_reassign_strategy                    28                         |
| |/                                                                          |
| O-> assign_nodes                                 29          29      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                     1                         |
|  /                                                                          |
| O-> electronic_write_occupancies                  1           1      0.00s  |
|    /                                                                        |
|   o-> electronic_analyse_occ_all                  1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calc_storage                         3                         |
|  /                                                                          |
| O-> locps_calc_storage                            3           3      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_cmplx                          3           3      0.00s  |
|   o-> pot_calc_storage                            1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_write                                 2                         |
|  /                                                                          |
| O-> comms_barrier                                 2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_sedc_initialize                        1                         |
|  /o-- mbd2_calculate                             14                         |
| |/                                                                          |
| O-> sedc_get_default_params                      15          15      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_order_eigenvalues               14                         |
|  /                                                                          |
| O-> comms_copy_gv_integer                        14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /o-- cell_generate_qpoints_local                 2                         |
| |/                                                                          |
| O-> cell_detect_MP                                3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_pseudo_scf                         8                         |
|  /                                                                          |
| O-> ion_atom_basis_pseudo_dealloc                 8           8      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_new_calc                              32                         |
|  /                                                                          |
| O-> cell_pack_vec_real                           32          32      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> model_deallocate                              1           1      0.00s  |
|    /                                                                        |
|   o-> model_reset                                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    14                         |
|  /                                                                          |
| O-> electronic_write_spin_density                14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_write_occupancies                1                         |
|  /                                                                          |
| O-> electronic_analyse_occ_all                    1           1      0.00s  |
|    /                                                                        |
|   o-> comms_copy_bnd_integer                      2           2      0.00s  |
|   o-> comms_reduce_kp_integer                     2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    12                         |
|  /                                                                          |
| O-> hubbard_calculate_occ_matrix                 12          12      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_count_symmetry_translations_wrap       1                         |
|  /o-- cell_output_wrapped                         1                         |
|  /o-- density_symmetrise_gv_parallel             25                         |
| |/                                                                          |
| O-> cell_factor_group_symmetry_wrapped           27          27      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_calculate_forces_r                    1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
| |/                                                                          |
| O-> wave_beta_phi_wv                              2           2      0.00s  |
|    /                                                                        |
|   o-> wave_beta_phi_wv_ks                         2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- mbd2_fd                                    12                         |
|  /                                                                          |
| O-> cell_reduce_symmetry_perturb                 12          12      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_new_calc                              16                         |
|  /                                                                          |
| O-> cell_frac_to_cart_cell                       16          16      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    12                         |
|  /                                                                          |
| O-> hubbard_mix_occ_matrix                       12          12      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_check_cell_constraints                 1                         |
|  /o-- cell_output_wrapped                         1                         |
|  /o-- locpot_calculate                           16                         |
| |/                                                                          |
| O-> cell_cart_lattice_to_abc                     18          18      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_new_calc                              16                         |
|  /                                                                          |
| O-> get_xc_functional_name                       16          16      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- algor_uniform_random_array                  1                         |
|  /o-- parameters_read                             1                         |
|  /o-- wave_initialise_wv                          1                         |
| |/                                                                          |
| O-> algor_set_random_seed                         3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_generate_cell_constraints                1           1      0.00s  |
|    /                                                                        |
|   o-> algor_uniform_random_array                  1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /                                                                          |
| O-> dm_flush_history                              1           1      0.00s  |
|    /                                                                        |
|   o-> dm_finalise                                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /o-- cell_output_wrapped                         1                         |
| |/                                                                          |
| O-> cell_check_group                              2           2      0.00s  |
|    /                                                                        |
|   o-> algor_invert_real                           4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> hubbard_calculate_stress                      1           1      0.00s  |
|    /                                                                        |
|   o-> hubbard_initialise                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- mbd2_fd                                     3                         |
|  /                                                                          |
| O-> dftd_store_results                            3           3      0.00s  |
|    /                                                                        |
|   o-> cell_unpack_vec_real                        2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /o-- cell_supercell                              2                         |
|  /o-- mbd2_fd                                    12                         |
| |/                                                                          |
| O-> cell_recip_lattice                           15          15      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hirshfeld_calculate                         2                         |
|  /                                                                          |
| O-> hirshfeld_init                                2           2      0.00s  |
|    /                                                                        |
|   o-> parameters_nspins                           4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           2                         |
|  /                                                                          |
| O-> cell_symmetry_test                            2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_allocate_nl_d                         3                         |
|  /o-- electronic_finalise                         1                         |
|  /o-- nlpot_calculate_forces_r                    1                         |
|  /o-- nlpot_calculate_stress_r                    1                         |
| |/                                                                          |
| O-> nlpot_deallocate_nl_d                         6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- write_eigenvalues                           9                         |
|  /                                                                          |
| O-> global_kpoint_index                           9           9      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_forces                     1                         |
|  /                                                                          |
| O-> hubbard_calculate_forces                      1           1      0.00s  |
|    /                                                                        |
|   o-> hubbard_initialise                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_reset                                 2                         |
|  /o-- electronic_finalise                         4                         |
|  /o-- popn_calculate_mulliken                     1                         |
| |/                                                                          |
| O-> wave_deallocate_wv                            7           7      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_flush_history                            1                         |
|  /                                                                          |
| O-> dm_finalise                                   1           1      0.00s  |
|    /                                                                        |
|   o-> density_deallocate                          1           1      0.00s  |
|   o-> dm_mix_density_deallocate                   6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_supercell                              2                         |
|  /                                                                          |
| O-> cell_copy_kpoints                             2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_magres_data                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_atom_init_pseudo_atom                   8                         |
|  /                                                                          |
| O-> ion_atom_resolve_pseudo_cfg                   8           8      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_minimisation                    14                         |
|  /                                                                          |
| O-> electronic_dump                              14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- basis_initialise                            6                         |
|  /                                                                          |
| O-> basis_utils_prime_factors                     6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           2                         |
|  /o-- cell_supercell                              2                         |
|  /o-- mbd2_fd                                    12                         |
| |/                                                                          |
| O-> cell_calculate_volume                        16          16      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_supercell                              2                         |
|  /                                                                          |
| O-> cell_generate_supercell_origins               2           2      0.00s  |
|    /                                                                        |
|   o-> cell_num_supercells                         2           2      0.00s  |
|   o-> algor_invert_real                           2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- nlpot_prepare_precon_ks                    14                         |
|  /                                                                          |
| O-> comms_copy_bnd_logical                       14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_sedc_initialize                        1                         |
|  /o-- mbd2_calculate                             14                         |
| |/                                                                          |
| O-> dftd_read_custom_params                      15          15      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_set_supercell_symmetry                 2                         |
|  /                                                                          |
| O-> cell_reduce_symmetry_supercell                2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_spectral_data                       1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- locpot_calc_storage                         3                         |
|  /                                                                          |
| O-> hartree_calc_storage                          3           3      0.00s  |
|    /                                                                        |
|   o-> algor_sizeof_real                           3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- mbd2_calculate                             14                         |
|  /                                                                          |
| O-> comms_save_strategy                          14          14      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- phonon_store_mem_estimate                   3                         |
|  /                                                                          |
| O-> secondd_store_mem_estimate                    3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_bs_data                             1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_sedc_initialize                        1                         |
|  /o-- mbd2_calculate                             14                         |
| |/                                                                          |
| O-> dftd_check_params                            15          15      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_elnes_data                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- bib_output                                  1                         |
|  /                                                                          |
| O-> comms_reduce_farm_logical                     1           1      0.00s  |
|    /                                                                        |
|   o-> comms_lcopy                                 1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_initialise                            1                         |
|  /                                                                          |
| O-> model_init_occ_eigenvalues_ef                 1           1      0.00s  |
|    /                                                                        |
|   o-> wave_spin_type_wv                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hubbard_initialise                          3                         |
|  /                                                                          |
| O-> hubbard_is_on                                 3           3      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_supercell                              2                         |
|  /o-- cell_generate_supercell_origins             2                         |
|  /o-- cell_set_supercell_symmetry                 2                         |
| |/                                                                          |
| O-> cell_num_supercells                           6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /o-- locps_calculate_stress                      1                         |
| |/                                                                          |
| O-> locps_calculate_non_coulomb                   2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- hirshfeld_init                              4                         |
|  /o-- hirshfeld_output                            1                         |
| |/                                                                          |
| O-> parameters_nspins                             5           5      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- check_forces_stresses                       1                         |
|  /                                                                          |
| O-> firstd_symmetrise_forces                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- parameters_validate                         7                         |
|  /                                                                          |
| O-> parameters_xc_allowed                         7           7      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_output_wrapped                         2                         |
|  /                                                                          |
| O-> cell_symmetry_symbol_wrapped_wrapped          2           2      0.00s  |
|    /                                                                        |
|   o-> algor_invert_real                           2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- firstd_calculate_stress                     1                         |
|  /                                                                          |
| O-> tddft_calculate_stress                        1           1      0.00s  |
|    /                                                                        |
|   o-> tddft_set_tddft_on                          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_check_cell_constraints                   1           1      0.00s  |
|    /                                                                        |
|   o-> cell_cart_lattice_to_abc                    1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- model_write                                 2                         |
|  /                                                                          |
| O-> comms_barrier_farm                            2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_count_symmetry_translations_wrappe       1           1      0.00s  |
|    /                                                                        |
|   o-> cell_factor_group_symmetry_wrapped          1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_check_group_equiv                        1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /o-- tddft_calculate_stress                      1                         |
| |/                                                                          |
| O-> tddft_set_tddft_on                            2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dftd_store_results                          2                         |
|  /                                                                          |
| O-> cell_unpack_vec_real                          2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_detect_primitive                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- check_forces_stresses                       1                         |
|  /                                                                          |
| O-> firstd_symmetrise_stress                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- comms_reduce_farm_logical                   1                         |
|  /o-- comms_reduce_bnd_logical                    1                         |
|  /o-- comms_reduce_gv_logical                     1                         |
|  /o-- comms_reduce_kp_logical                     1                         |
| |/                                                                          |
| O-> comms_lcopy                                   4           4      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- ion_clebsch_gordan                          1                         |
|  /                                                                          |
| O-> init_factorial                                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_calculate_charge_spilling              1                         |
|  /o-- popn_calculate_density_matrix               1                         |
| |/                                                                          |
| O-> popn_nbands_occupied                          2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- density_initialise_cell                     1                         |
|  /                                                                          |
| O-> comms_copy_bnd_real                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- dm_finalise                                 6                         |
|  /                                                                          |
| O-> dm_mix_density_deallocate                     6           6      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_supercell_data                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_analyse_occ_all                  2                         |
|  /                                                                          |
| O-> comms_copy_bnd_integer                        2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_read_optics_data                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /o-- electronic_initialise                       1                         |
| |/                                                                          |
| O-> ion_real_initialise                           2           2      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- popn_create_lcao_basis_allkpt               1                         |
|  /                                                                          |
| O-> popn_check_lcao_basis                         1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- find_strategy                               1                         |
|  /                                                                          |
| O-> best_mixed_strategy                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_check_keywords                           1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- comms_parallel_strategy                     1                         |
|  /                                                                          |
| O-> reassign_nodes                                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- castep                                      1                         |
|  /                                                                          |
| O-> memory_system_initialise                      1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- cell_read_wrapped                           1                         |
|  /                                                                          |
| O-> cell_check_symmetry_ops_spin                  1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_finalise                         1                         |
|  /                                                                          |
| O-> hubbard_finalise                              1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- xc_calculate_stress_density                 1                         |
|  /                                                                          |
| O-> xc_calculate_stress_correction                1           1      0.00s  |
+-----------------------------------------------------------------------------+
|   o-- electronic_initialise                       1                         |
|  /                                                                          |
| O-> electronic_restore                            1           1      0.00s  |
+-----------------------------------------------------------------------------+
Class of operation                  Time spent
COMMS                                 230.79s
COMMS_GV                              226.32s
COMMS_KP                                3.57s
COMMS_BND                               0.45s
COMMS_FARM                              0.00s
     508 different subroutines and functions were traced
Hash collisions:         92
  
