# Statistical-analysis-of-Galaxies-Clusters

This is an academic project as part of my course ST522. Computational Statistics at University of Southern Denmark, where we are required to apply statistical methods we’ve learnt on a self-sourced dataset, formulate a problem statement and find a solution to that problem.

Through conducted secondary research on scholarly articles and peer-reviewed journals, our group have decided to base our dataset on 820 galaxies in the Corona Borealis region. Our problem statement, thus, is whether these galaxies belongs to cluster(s) and if so, how many.

From preliminary data visualization, we concluded there are 3 major clusters in this dataset. We then decided to use Bayesian inference (A major topic in the course) to help determine which galaxies belong to which cluster in the following steps:
1. Apply Bayes’s Theorem to model the relationship between the observed data (galaxies) and the latent variable (clusters)
2. Use the posterior distribution derived from Bayes’s Theorem along with the observed data in Gibb’s sampling to simulate a Markov chain, allowing for inference of latent variables and mixing proportions.

While we managed to estimate and validate our results, we later realized this can also be solved using EM Algorithm, which works as followed:
1. In the Expectation (E) step, the Q function, which is the log likelihood function of observing the data given the parameters, is calculated
2. In the Maximization (M) step, the parameters that maximizes the above Q function is calculated.
   
The 2 steps are then repeated until convergence condition is met. With this method, we were able to achieve the same result but with less computational time.
