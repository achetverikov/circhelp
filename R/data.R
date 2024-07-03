#' Data from an orientation estimation task
#'
#' A dataset with the orientation estimation results from Experiment 2 in Pascucci et al. (2019) available from https://doi.org/10.5281/zenodo.2544946.
#'
#' @format A data frame with 4400 rows and 5 variables:
#' \describe{
#'   \item{observer}{observer ID}
#'   \item{orientation}{true orientation}
#'   \item{reported}{reported orientation}
#'   \item{rt}{response time}
#'   \item{err}{estimation error}
#' }
#' @references {
#' Pascucci, D., Mancuso, G., Santandrea, E., Libera, C. D., Plomp, G., & Chelazzi, L. (2019). Laws of concatenated perception: Vision goes for novelty, decisions for perseverance. PLoS Biology, 17(3). \doi{10.1371/journal.pbio.3000144}
#' }
#'
#' @source \url{https://zenodo.org/record/2544946/files/Experiment2_rawdata.csv?download=1}
"Pascucci_et_al_2019_data"

#' Data from a motion estimation task
#'
#' A dataset with the motion estimation results from Bae & Luck (2018) available from https://osf.io/4m2kb/ (some variables are removed, see the link for the full dataset).
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
#' @references {
#' Bae, G.-Y., & Luck, S. J. (2018). Decoding motion direction using the topography of sustained ERPs and alpha oscillations. NeuroImage, 184(August 2018), 242–255. \doi{10.1016/J.NEUROIMAGE.2018.09.029}
#' }
#' @source \url{https://osf.io/4m2kb/download}
"Bae_Luck_2018_data"
