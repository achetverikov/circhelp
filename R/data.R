#' Data from an orientation estimation task
#'
#' A dataset with the orientation estimation results from Experiment 2 in Pascucci et al. (2019, PLOS Biology, https://dx.doi.org/10.1371/journal.pbio.3000144) available from https://doi.org/10.5281/zenodo.2544946.
#'
#' @format A data frame with 4400 rows and 5 variables:
#' \describe{
#'   \item{observer}{observer ID}
#'   \item{orientation}{true orientation}
#'   \item{reported}{reported orientation}
#'   \item{rt}{response time}
#'   \item{err}{estimation error}
#' }
#' @source \url{https://zenodo.org/record/2544946/files/Experiment2_rawdata.csv?download=1}
"Pascucci_et_al_2019_data"

#' Data from a motion estimation task
#'
#' A dataset with the motion estimation results from Bae & Luck (2018, Neuroimage, https://doi.org/10.1016/j.neuroimage.2018.09.029) available from https://osf.io/4m2kb/ (some variables are removed, see the link for the full dataset).
#'
#' @format A data frame with 20480 rows and 8 variables:
#' \describe{
#'   \item{subject_Num}{observer ID}
#'   \item{trial_Num}{trial number}
#'   \item{TargetDirection}{true motion direction}
#'   \item{RespAngle}{reported motion direction}
#'   \item{motionCoh}{motion coherence}
#'   \item{Block}{block number}
#'   \item{Session}{session number}
#'   \item{err}{estimation error}
#' }
#' @source \url{https://osf.io/4m2kb/download}
"Bae_Luck_2018_data"
