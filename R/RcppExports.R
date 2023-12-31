# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

my_abs <- function(x) {
    .Call(`_evinf_my_abs`, x)
}

ell_nb_i_fun <- function(beta_nb, alpha_nb, x_nb_ext_i, y_i) {
    .Call(`_evinf_ell_nb_i_fun`, beta_nb, alpha_nb, x_nb_ext_i, y_i)
}

delldtheta_nb_i_fun <- function(beta_nb, alpha_nb, x_nb_ext_i, y_i) {
    .Call(`_evinf_delldtheta_nb_i_fun`, beta_nb, alpha_nb, x_nb_ext_i, y_i)
}

d2elldtheta2_nb_i_fun <- function(beta_nb, alpha_nb, x_nb_ext_i, y_i) {
    .Call(`_evinf_d2elldtheta2_nb_i_fun`, beta_nb, alpha_nb, x_nb_ext_i, y_i)
}

ell_pl_i_fun <- function(beta_pl, c_pl, x_pl_ext_i, y_i) {
    .Call(`_evinf_ell_pl_i_fun`, beta_pl, c_pl, x_pl_ext_i, y_i)
}

delldbeta_pl_i_fun_approx <- function(beta_pl, c_pl, x_pl_ext_i, y_i) {
    .Call(`_evinf_delldbeta_pl_i_fun_approx`, beta_pl, c_pl, x_pl_ext_i, y_i)
}

d2elldbeta2_pl_i_fun_approx <- function(beta_pl, c_pl, x_pl_ext_i, y_i) {
    .Call(`_evinf_d2elldbeta2_pl_i_fun_approx`, beta_pl, c_pl, x_pl_ext_i, y_i)
}

log_lik_fun <- function(gamma_z, gamma_pl, beta_nb, alpha_nb, beta_pl, c_pl, x_mult_z_ext, x_mult_pl_ext, x_nb_ext, x_pl_ext, y) {
    .Call(`_evinf_log_lik_fun`, gamma_z, gamma_pl, beta_nb, alpha_nb, beta_pl, c_pl, x_mult_z_ext, x_mult_pl_ext, x_nb_ext, x_pl_ext, y)
}

update_bfgs_fun <- function(gamma_z_in, gamma_pl_in, beta_nb_in, alpha_nb_in, beta_pl_in, c_pl, x_mult_z_ext, x_mult_pl_ext, x_nb_ext, x_pl_ext, y, max_upd_par, no_m_bfgs_steps) {
    .Call(`_evinf_update_bfgs_fun`, gamma_z_in, gamma_pl_in, beta_nb_in, alpha_nb_in, beta_pl_in, c_pl, x_mult_z_ext, x_mult_pl_ext, x_nb_ext, x_pl_ext, y, max_upd_par, no_m_bfgs_steps)
}

