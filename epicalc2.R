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
  object <- matrix(c(a, c, b, d), nrow=2)
  
  structure(
    object,
    class = c(class, "contingency")
  )
}

# convert a dataframe into a contingency table
as.contingency <- function(object,
                           exposures = character(),
                           cases = character(),
                           ...,
                           class = character()) {
  if (inherits(object, "data.frame")) {
    if (length(exposures) == 0) {
      stop("name of exposure column must be provided if initializing from a data.frame")
    } else if (length(cases) == 0) {
      stop("name of cases column must be provided if initializing from a data.frame")
    }
    
    contingency(
      a = sum(object[exposures] == 1 & object[cases] == 1),
      b = sum(object[exposures] == 0 & object[cases] == 1),
      c = sum(object[exposures] == 1 & object[cases] == 0),
      d = sum(object[exposures] == 0 & object[cases] == 0),
      class = class
    )
  } else {
    stop("unsupported class")
  }
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
contingency.count <- function(...) {
  contingency(..., class="count")
}

# 2x2 table subclass for person-time data
contingency.persontime <- function(...) {
  contingency(..., class="persontime")
}

# 2x2 table subclass for case-control data
contingency.casecontrol <- function(...) {
  contingency(..., class = "casecontrol")
}
as.casecontrol <- function(...) as.contingency(..., class = "casecontrol")
summary.casecontrol <- function(object,
                                alpha = 0.05,
                                ...) {
  # constants
  a <- object[1,1]
  b <- object[1,2]
  c <- object[2,1]
  d <- object[2,2]
  m.1 <- a + b
  m.0 <- c + d
  n.1 <- a + c
  n.0 <- b + d
  total <- m.1 + m.0
  
  z <- qnorm(1 - alpha/2)
  
  # odds ratio and confidence interval
  odds.ratio <- (a*d) / (b*c)
  x <- log(odds.ratio)
  variance <- 1/a + 1/b + 1/c + 1/d
  error <- z * variance
  odds.ratio.ci <- exp(c(x - error, x + error))
  names(odds.ratio.ci) <- c("lower", "upper")
  
  # test of association
  x <- a
  expected <- n.1 * m.1 / total
  variance <- (m.1 * m.0 * n.1 * n.0) / (total^2 * (total - 1))
  z.sq <- (x - expected)^2 / variance
  p.value <- 1 - pchisq(z.sq, 1)
  
  structure(
    list(
      matrix = object,
      alpha = alpha,
      odds.ratio = odds.ratio,
      odds.ratio.ci = odds.ratio.ci,
      z.sq = z.sq,
      p.value = p.value
    ),
    class = "summary.casecontrol"
  )
}
print.summary.casecontrol <- function(x, ...) {
  print(as.data.frame(x$matrix))
  cat(
    sprintf("\nUsing a critical value of alpha = %f", x$alpha),
    "\n",
    sprintf(
      "\nOdds ratio:\t%f (%f, %f)",
      x$odds.ratio,
      x$odds.ratio.ci["lower"],
      x$odds.ratio.ci["upper"]
    ),
    sprintf("\nChi-sq test for association:\tp"),
    ifelse(
      x$p.value < 1e-6,
      "< 0.000001", 
      sprintf("= %f", x$p.value)
    )
  )
  invisible(x)
}

# testing
#test <- contingency.casecontrol(623, 519, 1062, 1641)
#summary(test)
#as.data.frame(test)
#as.data.frame(as.contingency(data, "htn", "dead"))