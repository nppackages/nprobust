## Shared kernel helper + fixture dataset
w_fun <- function(u, kernel) {
  switch(kernel,
         "uni" = 0.5 * (abs(u) <= 1),
         "tri" = (1 - abs(u)) * (abs(u) <= 1),
         "epa" = 0.75 * (1 - u^2) * (abs(u) <= 1))
}

make_fixture <- function(n = 1000, seed = 7) {
  set.seed(seed)
  x <- sort(runif(n, -1, 1))
  y <- 1 + 2 * x + 3 * x^2 - 0.5 * x^3 + rnorm(n, sd = 0.3)
  list(x = x, y = y, n = n)
}
