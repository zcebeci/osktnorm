\name{plot_functions}
\alias{plot_trajectories}
\alias{plot_tracking_errors}
\alias{plot_consensus_stats}
\title{Visualization of Simulation Results}
\description{
  These functions provide visualization tools to reproduce the graphical results presented in Chapter 3 of the thesis. They use \code{ggplot2} to plot agent trajectories, tracking errors, and statistical consensus rates.
}
\usage{
plot_trajectories(history_pos, leader_history, is_informed, T_step)

plot_tracking_errors(history_pos, leader_history, is_informed, T_step)

plot_consensus_stats(results_df)
}
\arguments{
  \item{history_pos}{A numeric array of dimensions \eqn{steps \times N \times 2} containing the position history of all agents over the simulation duration.}
  
  \item{leader_history}{A numeric matrix of dimensions \eqn{steps \times 2} containing the position history of the virtual leader.}
  
  \item{is_informed}{A logical vector of length \eqn{N} indicating which agents are informed (TRUE) and which are uninformed (FALSE).}
  
  \item{T_step}{A numeric scalar representing the sampling time step used in the simulation (e.g., 0.1).}
  
  \item{results_df}{A data frame summarizing Monte Carlo simulation results. Must contain columns: \code{Radius}, \code{AgentCount}, \code{ConsensusRate}, and \code{InformedPercentage}.}
}
\details{
  \itemize{
    \item \strong{plot_trajectories:} Visualizes the movement paths of the agents and the virtual leader in X and Y dimensions. This function reproduces the style of \strong{Figures 3.5 - 3.8} in the thesis. Informed agents are plotted in blue, uninformed in black/gray, and the virtual leader in dashed red.
    
    \item \strong{plot_tracking_errors:} Plots the position error (\eqn{x_i - x_\gamma}) over time for each agent. This function reproduces the style of \strong{Figures 3.1 - 3.4} in the thesis. It allows visual verification of whether the swarm successfully tracks the leader (error converges to zero).
    
    \item \strong{plot_consensus_stats:} Visualizes the relationship between consensus rate, sensing radius, and the percentage of informed agents. This reproduces the statistical analysis plots found in \strong{Figures 3.33 - 3.34}.
  }
}
\value{
  A \code{ggplot} object that can be printed or modified further.
}
\author{
  Cagatay Cebeci
}
\seealso{
  \code{\link{prepare_plot_data}}
}
\examples{
# --- Example for Trajectory and Error Plots ---
# Assuming 'history_pos' and 'leader_history' are generated from a simulation:
N <- 20
is_informed <- rep(FALSE, N); is_informed[1:2] <- TRUE # 10\% informed
T_step <- 0.1

\dontrun{
# Plot Trajectories
p1 <- plot_trajectories(history_pos, leader_history, is_informed, T_step)
print(p1)

# Plot Errors
p2 <- plot_tracking_errors(history_pos, leader_history, is_informed, T_step)
print(p2)
}

# --- Example for Statistical Plot ---
# Mock data representing thesis Table 3.1 results
results <- data.frame(
  Radius = rep(c(5, 7, 10, 20), each = 2),
  AgentCount = rep(c(20, 50), 4),
  ConsensusRate = runif(8, 0.5, 1.0),
  InformedPercentage = 10
)

p3 <- plot_consensus_stats(results)
print(p3)
}
