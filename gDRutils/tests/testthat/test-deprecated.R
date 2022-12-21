test_that("logistic_4parameters works as expected", {
    c <- 1
    Vinf <- 0.1
    V0 <- 1
    h <- 2
    EC50 <- 0.5

    # Non-numeric values cause an error.
    expect_error(predict_efficacy_from_conc(
        c = "non-numeric_entry",
        Vinf = Vinf,
        V0 = V0,
        EC50 = EC50,
        h = h
    ))

    # Normal fit.
    v <- predict_efficacy_from_conc(
        c = c,
        x_inf = Vinf,
        x_0 = V0,
        ec50 = EC50,
        h = h
    )
    expect_equal(v, 0.28)

    # Flat fit.
    EC50 <- 0
    v <- predict_efficacy_from_conc(
        c = c,
        x_inf = Vinf,
        x_0 = V0,
        ec50 = EC50,
        h = h
    )
    expect_equal(v, Vinf)
})

test_that("support deprecated get_identifiers", {
  obs <- get_env_identifiers("drug_name")
  expect_equal(obs, "DrugName")
})
