# Base 2x2 table class
contingency <- function(a = integer(),
                        b = integer(),
                        c = integer(),
                        d = integer(),
                        ...,
                        class = character()) {
  # Primary data structure is that of a 2x2 array
  #object <- array(c(a, b, c, d), dim = c(2, 2))
  # Primary data structure is that of a 2x2 matrix
  object <- matrix(c(a, b, c, d), nrow=2)
  
  structure(
    object,
    class = c(class, "contingency")
  )
}

# Convert any contingency table into a dataframe and label the rows + columns
as.data.frame.contingency <- function(object,
                                      top.row = 1,
                                      bottom.row = 2,
                                      ...) {
  # calculation of column totals
  exposed <- object[,1]
  exposed[3] <- sum(exposed)
  unexposed <- object[,2]
  unexposed[3] <- sum(unexposed)
  
  # building the dataframe
  df <- data.frame(
    exposed = exposed,
    unexposed = unexposed,
    total = exposed + unexposed # calculation of row totals
  )
  
  # rownames are assigned based on what kind of contingency table we have
  if (inherits(object, "count")) {
    row.names <- c("cases", "noncases", "total")
  } else if (inherits(object, "persontime")) {
    row.names <- c("cases", "person.time", "total")
  } else if (inherits(object, "casecontrol")) {
    row.names <- c("cases", "controls", "total")
  } else {
    row.names <- c(top.row, bottom.row, "total")
  }
  rownames(df) <- row.names

  # if we are returning a person-time table then there are no columnwise totals
  if (inherits(object, "persontime")) {
    df[-3,]
  } else {
    df
  }
}

# 2x2 table subclass for count data
contingency.count <- function(a = integer(),
                              b = integer(),
                              c = integer(),
                              d = integer()) {
  object <- contingency(a, b, c, d, class="count")
  object
}

# 2x2 table subclass for person-time data
contingency.persontime <- function(a = integer(),
                              b = integer(),
                              c = integer(),
                              d = integer()) {
  object <- contingency(a, b, c, d, class="persontime")
  object
}

# 2x2 table subclass for case-control data
contingency.casecontrol <- function(a = integer(),
                                    b = integer(),
                                    c = integer(),
                                    d = integer()) {
  object <- contingency(a, b, c, d, class="casecontrol")
  object
}
summary.casecontrol <- function(object,
                                alpha = 0.05,
                                ...) {
  print(as.data.frame(object))
  
  # constants
  a <- test[1,1]
  b <- test[1,2]
  c <- test[2,1]
  d <- test[2,2]
  m.1 <- a + b
  m.0 <- c + d
  n.1 <- a + c
  n.0 <- b + d
  total <- m.1 + m.0
  
  z <- qnorm(1 - 0.05/2)
  
  # odds ratio and confidence interval
  odds.ratio <- (a*d) / (b*c)
  x <- log(odds.ratio)
  variance <- 1/a + 1/b + 1/c + 1/d
  error <- z * variance
  odds.ratio.ci <- exp(c(x - error, x + error))
  
  # test of association
  x <- a
  expected <- n.1 * m.1 / total
  variance <- (m.1 * m.0 * n.1 * n.0) / (total^2 * (total - 1))
  z.sq <- (x - expected)^2 / variance
  p.association <- 1 - pchisq(z.sq, 1)
  
  # print
  cat(
    sprintf("\nUsing a critical value of alpha = %f", alpha),
    "\n",
    sprintf("\nOdds ratio:\t%f (%f, %f)", odds.ratio, odds.ratio.ci[1], odds.ratio.ci[2]),
    sprintf("\nChi-sq test for association:\tp"),
    ifelse(
      p.association < 1e-6,
      "< 0.000001", 
      sprintf("= %f", p.assoc)
    )
  )
}

# testing
test <- contingency.casecontrol(623, 519, 1062, 1641)
summary(test)
#as.data.frame(test)