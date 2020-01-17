
# regression funs ---------------------------------------------------------


#' NB regression wrapper
#'
#' It performs the negative bionomial regression from the MASS package.
#' It also calculates the confindence intervals using multiple methods
#'
#' @param data A dataframe with counts and factors
#' @param formula A formula (as string or formula class)
#' @param control a glm.control object
#' @param ci_parameters parameters to compute confidence intervals
#' @param ci_method method to compute confidence intervals
#' @param alpha signifocance level for confidence intervals
#' @param name name of the experiment (logging)
#'
#' @return A data frame with coefficients and confidence intervals
#' @export
#'
#' @examples
#'
#' # not yet
#'
glm_nb_wrapper <- function(data,
                           formula,
                           control = glm.control(),
                           ci_parameters = NULL,
                           ci_method = c("profile","z","none"),
                           alpha = 0.05,
                           name = "Regression"){
  # needs to always return a dataFrame

  try_catch_error_mssg = "no valid set of coefficients has been found: please supply starting values"

  stopifnot(is.numeric(alpha) && alpha < 1 )

  # defaults to profile
  if(length(ci_method)>1){
    ci_method = ci_method[1]
  }

  if (!is(formula,"formula")){
    formula = as.formula(formula)
  }

  fit = tryCatch(
    {MASS::glm.nb(formula = formula,
                  data = data,
                  control = control)
    },error = function(err){
      warning(glue::glue("{name} didn't converge"))
      warning(err)
      NULL
    },finally = {})

  if (is.null(fit)){
    fit_df =  NULL

    ci_method = "none"
  } else {
    fit_df = fit %>% broom::tidy()
  }

  switch(ci_method,
         profile = {
           tryCatch({

             ci_obj = confint(fit,
                              level = (1 - alpha))

             ci_obj %>%
               as.data.frame() %>%
               tibble::rownames_to_column("term") %>%
               magrittr::set_colnames(c("term", "ci_low", "ci_high")) -> CI_df

           }, error = function(err) {
             # this is not the ideal way to handle code but is there any
             # better?
             if (err == try_catch_error_mssg){
               warning(glue::glue("{name} CI didn't work"))
               warning(err)
               CI_df = NULL
             } else {
               stop(err)
             }
           },
           finally = {
           })
         },
         z = {
           # here I compute the confidence intervals based on the std.error,
           # based on the assumption that estimates will be normally distributed
           zVal_low = qnorm(p = alpha / 2)
           zVal_high = qnorm(p = 1 - (alpha / 2))
           fit_df %>% dplyr::mutate(ci_low = estimate + zVal_low * std.error,
                                    ci_high = estimate + zVal_high * std.error) %>%
             dplyr::select(term, ci_low, ci_high) -> CI_df # rather stupid
         },
         none = {
           CI_df = NULL
         })

  if (! "CI_df" %in% ls(envir = environment())){
    # this is really weird.. no idea why is this happening...
    CI_df = NULL
  }
  if (is.null(CI_df)){
    CI_df = data.frame(term = fit_df$term,
                       ci_low = NA,
                       ci_high = NA)
  }

  result_df = dplyr::left_join(fit_df,CI_df,by = "term")
  result_df
}

