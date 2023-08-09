library(ggplot2)

#Helper function to search the file path of the RData file (So we won't have to manually input the path)
find_RData_path <- function(file_name, search_dirs = getwd()) {
  file_path = NULL
  
  for (dir in search_dirs) {
    files = list.files(dir, full.names = TRUE)
    matching_files = files[grep(file_name, files)]
    
    if (length(matching_files) > 0) {
      file_path = matching_files[1]
      break
    }
  }
  return(file_path)
}

# Part 1
set.seed(2023)
# Set parameters
n = 10^5;
μ1 = 4;
μ2 = 8;
Σ1 = 1;
Σ2 = 1;
a=1;
b=1;

# (b)

# Loading the data
file_name_1 = "x.RData";
path_1 = find_RData_path(file_name_1)
load(path_1)

# Initialize the parameters p and Z with respectively Beta function (1,1) 
# and 200 = numbers of x_i Binomial functions (1,p)
p_0 = rbeta(1,shape1 = 1, shape2 = 1);
Z_0 = rbinom(length(x),1,p_0);

# Define the markov chain matrix with the initialized values in first row
markov_chain = matrix(nrow = n,ncol = length(x) + 1);
markov_chain[1, ] = c(p_0, Z_0);

# Define a function which generates values of p using the derivation of p|X_1,...,X_n,Z_1,...,Z_n
# found in (a)
P_Generator <- function(x,z) {
  p = rbeta(1, shape1 = a + sum(z), shape2 = b + length(z) - sum(z));
  return (p);
}

# Define a function which generates values of Z using the derivation of Z_i|X_1,...,X_n,p
# found in (a)
Z_Generator <- function(n,x,p) {
  Z <- vector("numeric", length = n);
  for (i in 1:n) {
    # Calculate the probabilities of Z_i being 1 with the given parameter
    p_i <- (p * dnorm(x[i], μ1, Σ1)) / ((p * dnorm(x[i], μ1, Σ1)) + ((1 - p) * dnorm(x[i], μ2, Σ2)))
    # Simulate a value of Z_i based on the calculated probabilities (inverse method)
    Z[i] <- ifelse(runif(1) < p_i, 1, 0)
  }
  return (Z);
}

# Gibbs Sampling algorithm: we initialize values for p and Z and then we use them to update their next values
# Note that we always use the most recent value, that is why we generate the values of Z with the new value of p
for (i in 2:n) {
  # Generate the value of p in first column of each line with all the precedent values of Z 
  markov_chain[i,1] = P_Generator(x, markov_chain[i-1, 2:201]);
  # Calculate all the new values of Z given the p value obtained just above
  markov_chain[i, 2:201] = Z_Generator(length(x), x, markov_chain[i,1]);
}

# Visualize p with the burn-in period by using an histogram
p_values = markov_chain[ ,1];
hist(p_values, breaks = 20, col = "red", xlab = "p", main = "Histogram of p values")

# Discard the burn-in period
burn_in = 10^3;
markov_chain_updated = markov_chain[(burn_in+1): n, ];

# Visualize p without the burn-in period
p_values = markov_chain_updated[ ,1];
hist(p_values, breaks = 20, col = "blue", xlab = "p", main = "Histogram of p values")




# (c)

Z_values = markov_chain_updated[ , 2:201];

# Calculate the posterior mean as an estimator under square loss
posterior_mean = (a + sum(Z_values)) / (a + b + length(Z_values));

# Calculate the credible interval by using the quantile function
lower_bound = qbeta(0.025, 1 + sum(Z_values), 1+length(Z_values) - sum(Z_values))
upper_bound = qbeta(0.975, 1 + sum(Z_values), 1+length(Z_values) - sum(Z_values))

credible_interval <- c(lower_bound, upper_bound)

# Print the results in the console
print(paste("Posterior mean estimate of p:", posterior_mean))
print(paste("95% Credible Interval of p:", credible_interval[1], " and ", credible_interval[2]))




# (d)
# Stock the values of Z and calculate the means of the values of Z with the function colMeans
posterior_probs = colMeans(Z_values);

# Define the classes as vectors
class_1 = c();
class_2 = c();

# Classify the probabilities
for (i in 1:200) {
  if (posterior_probs[i] >= 0.5) {
    class_1 = c(class_1, x[i]);
  } else {
    class_2 = c(class_2, x[i]);
  }
}

# Print the means of the classes
print(paste(length(class_1)," elements are in class 1 while ",length(class_2)," elements are in class 2"));
print(paste(mean(class_1)," is the mean of elements in class 1 and  ",mean(class_2)," is the one for class 2"));


# Part 2

# Parameter setup
file_name_2 = "galaxies.RData";
path_2 = find_RData_path(file_name_2)
load(path_2)
n = length(galaxies);
k = 3;

# Function to randomize initial parameters, given the observations and number of components
set_initial_parameter<- function(X, k) {
  random_p = runif(k);
  
  # We need sum(p_j) = 1
  normalized_p = random_p/sum(random_p);
  
  # Generate µ and σ
  indices <- sample(length(X), k, replace = FALSE)
  random_µ <- X[indices]
  
  random_selection = sample(X, replace = TRUE);
  σ = sd(random_selection);

  return  (list("p" = normalized_p, "µ" = random_µ, "σ" = σ))
}

# Function to calculate Gamma function given the parameters
gamma_function <- function(p, µ, σ) {
  gamma_matrix = matrix(nrow = n, ncol = k);
  for (i in 1:n) {
    for (j in 1:k) {
      gamma_matrix[i,j] = (p[j] * dnorm(galaxies[i], mean = µ[j], sd = σ))/ sum(p* dnorm(galaxies[i], µ, sd = σ));
    }
  }
  return (gamma_matrix)
}

# Implementation of EM Algorithm
EM <-function(p, µ, σ){
  # Initialize θ and total_iter
  prev_θ = rep(0, length(c(p, µ, σ)));
  total_iter = 0
  
  while (TRUE){
    # Total_iter to keep track of the iterations
    total_iter = total_iter + 1;
    
    # E Step to calculate gamma function and Q function using the parameters
    gamma_matrix = gamma_function(p, µ, σ);
    current_Q = 0;
    for (i in 1:n) {
      for (j in 1:k) {
        current_Q = current_Q + gamma_matrix[i,j]*(log(p[j]) + log(dnorm(galaxies[i],µ[j], σ)));
      }
    }
    
    # M Step to calculate the parameters that maximize Q
    σ_square = 0;
    for (j in 1:k) {
      µ[j] =  sum(gamma_matrix[ ,j] * galaxies) / sum(gamma_matrix[ ,j]);
      p[j] = sum(gamma_matrix[ ,j])/ n;
      σ_square = σ_square + sum(gamma_matrix[ ,j]*(galaxies - µ[j])^2)/sum(gamma_matrix);
    }
    σ = sqrt(σ_square);
    
    
    # Assess convergence using relative criterion
    current_θ = c(p,µ, σ);
    if  (sqrt(sum((current_θ - prev_θ)^2)) / (sqrt(sum(prev_θ^2))) < 1e-6) {
      break;
    } else {
      prev_θ = current_θ;
    }
  }
  
  #Sort µ from smallest to largest and its corresponding p
  sorted_indexes = order(µ)
  p = p[sorted_indexes];
  µ = µ[sorted_indexes];
  gamma_matrix = gamma_matrix[ ,sorted_indexes];
  
  
  #Check which component observations of X belong to based on the parameters
  component_1 = c();
  component_2 = c();
  component_3 = c();
  for (i in 1:n) {
    value = max(gamma_matrix[i, ]);
    indices <- which(gamma_matrix == value, arr.ind = TRUE)
    if (indices[,2] == 1) {
      component_1 = c(component_1, galaxies[i])
    } else if (indices[,2] == 2) {
      component_2 = c(component_2, galaxies[i])
    } else {
      component_3 = c(component_3, galaxies[i])
    }
  }
  
  return (list("p" =p, "µ" =µ, "σ" =σ, "iter" =total_iter,
               "Component_1" = component_1,
               "Component_2" = component_2,
               "Component_3" = component_3));
}


#Get the initial parameter and run the EM algorithm
p_0 = set_initial_parameter(galaxies, k)$p;
µ_0 = set_initial_parameter(galaxies, k)$µ;
σ_0 = set_initial_parameter(galaxies, k)$σ;

EM(p_0, µ_0, σ_0)


#Extract the components
obs_1 = EM(p_0, µ_0, σ_0)$Component_1
obs_2 = EM(p_0, µ_0, σ_0)$Component_2
obs_3 = EM(p_0, µ_0, σ_0)$Component_3

# Combine the observations into a single data frame
df <- data.frame(
  observations = c(obs_1, obs_2, obs_3),
  group = factor(
    c(rep(1, length(obs_1)), rep(2, length(obs_2)), rep(3, length(obs_3))),
    levels = c(1, 2, 3)
  )
)

# Plot the density of the observations based on their components
ggplot(df, aes(x = observations, fill = group)) +
  geom_density(alpha = 0.5) +
  labs(x = "Observations", title = "Observations Plot based on EM Model") +
  scale_fill_manual(values = c("blue", "red", "green"), labels = c("Component 1", "Component 2", "Component 3"))




# Run the EM algorithm multiple times with randomly seclected starting points
# ans store the resulting parameters into a vector
iterations = 50;
table <- vector(mode="list", length = iterations);
j = 1;

for (i in 1:iterations) {
  
  #Randomnize the parameters
  p_0 = set_initial_parameter(galaxies, k)$p;
  µ_0 = set_initial_parameter(galaxies, k)$µ;
  σ_0 = set_initial_parameter(galaxies, k)$σ;
  
  #Store the resulting parameters and round them up
  result = EM(p_0, µ_0, σ_0)[1:3];
  rounded_result = lapply(result, function(x) round(x, digits = 2))
  
  #If the resulting parameters are not in the vector, we add them and set their count to 1
  # Else, if they already have an entry in the table, we increment their count
  if (!(list(rounded_result) %in% names(table)) ) {
    names(table)[j] = list(rounded_result);
    table[[j]] = 1;
    j = j + 1;
  } else {
    index = which(names(table) == list(rounded_result));
    table[[index]] = table[[index]] + 1;
  }
}








