# Regression test for factor_scores() score-alignment helper.
#
# Background: prior to 0.5.x .fs_apply_to_group set colnames(Mq) only when
# .fs_detect_perm() returned a non-NULL permutation match. For a genuine
# non-permutation Procrustes Q (any column with max|Q[,j]| <= 0.95) the
# returned matrix kept NULL column names, the downstream
# `for (cn in colnames(Mq)) scores[[cn]] <- Mq[, cn]` loop iterated zero
# times, and the rotation was silently discarded. The fix sets colnames(Mq)
# to colnames(M) in that branch so the rotated values are actually written
# back at the original factor positions.

test_that(".fs_detect_perm recognises identity and column swaps", {
  Q_id  <- diag(3)
  expect_identical(bifactory:::.fs_detect_perm(Q_id, c("A", "B", "C")),
                   c("A", "B", "C"))

  Q_swap <- diag(3)[, c(2, 1, 3)]
  expect_identical(bifactory:::.fs_detect_perm(Q_swap, c("A", "B", "C")),
                   c("B", "A", "C"))
})

test_that(".fs_detect_perm returns NULL for non-permutation rotation", {
  theta <- pi / 6  # cos(theta) = 0.866 < 0.95 threshold
  Q <- diag(3)
  Q[1:2, 1:2] <- matrix(c(cos(theta), -sin(theta),
                          sin(theta),  cos(theta)), 2, 2)
  expect_null(bifactory:::.fs_detect_perm(Q, c("A", "B", "C")))
})

test_that(".fs_apply_to_group applies identity Q without changing values", {
  M <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2,
              dimnames = list(NULL, c("A", "B", "C")))
  Mq <- bifactory:::.fs_apply_to_group(M, diag(3))
  expect_identical(colnames(Mq), c("A", "B", "C"))
  expect_equal(unname(Mq), unname(M))
})

test_that(".fs_apply_to_group relabels columns under a permutation Q", {
  M <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2,
              dimnames = list(NULL, c("A", "B", "C")))
  Q_swap <- diag(3)[, c(2, 1, 3)]
  Mq <- bifactory:::.fs_apply_to_group(M, Q_swap)
  # Columns of Mq are relabeled so that the downstream `scores[[cn]] <- Mq[, cn]`
  # loop puts factor B's rotated values back into scores$B etc.
  expect_identical(colnames(Mq), c("B", "A", "C"))
  expect_equal(unname(Mq), unname(M %*% Q_swap))
})

test_that(".fs_apply_to_group keeps colnames and applies sign-flip under signed permutation", {
  M <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2,
              dimnames = list(NULL, c("A", "B", "C")))
  Q_flip <- diag(c(-1, 1, 1))
  Mq <- bifactory:::.fs_apply_to_group(M, Q_flip)
  expect_identical(colnames(Mq), c("A", "B", "C"))
  expect_equal(unname(Mq[, "A"]), -unname(M[, "A"]))
  expect_equal(unname(Mq[, "B"]),  unname(M[, "B"]))
})

test_that(".fs_apply_to_group preserves colnames for non-permutation Q so the rotation is applied", {
  # Regression: silent-drop bug where colnames(Mq) was left NULL and the
  # downstream `for (cn in colnames(Mq))` loop wrote nothing back.
  set.seed(1L)
  n <- 50L
  M <- matrix(rnorm(n * 3), nrow = n,
              dimnames = list(NULL, c("A", "B", "C")))

  theta <- pi / 6
  Q <- diag(3)
  Q[1:2, 1:2] <- matrix(c(cos(theta), -sin(theta),
                          sin(theta),  cos(theta)), 2, 2)
  dimnames(Q) <- NULL  # mimic svd-derived Procrustes Q

  Mq <- bifactory:::.fs_apply_to_group(M, Q)

  expect_identical(colnames(Mq), c("A", "B", "C"))
  # Confirm the rotation is genuinely non-trivial in the score space.
  expect_gt(max(abs(Mq - M)), 0.1)
  # And that the rotated matrix matches M %*% Q exactly.
  expect_equal(unname(Mq), unname(M %*% Q))
})

test_that(".fs_apply_to_group downstream assignment loop writes rotated values back", {
  # Full pipeline check: simulate the assignment loop used inside
  # .fs_apply_alignment to confirm that with the fix, rotated values
  # actually reach the scores data.frame.
  set.seed(2L)
  n <- 30L
  M <- matrix(rnorm(n * 3), nrow = n,
              dimnames = list(NULL, c("A", "B", "C")))
  theta <- pi / 4  # cos = 0.707, well below 0.95
  Q <- diag(3)
  Q[1:2, 1:2] <- matrix(c(cos(theta), -sin(theta),
                          sin(theta),  cos(theta)), 2, 2)
  dimnames(Q) <- NULL

  scores <- as.data.frame(M)
  Mq <- bifactory:::.fs_apply_to_group(M, Q)
  for (cn in colnames(Mq)) scores[[cn]] <- Mq[, cn]

  expect_equal(scores[["A"]], unname(Mq[, "A"]))
  expect_equal(scores[["B"]], unname(Mq[, "B"]))
  expect_equal(scores[["C"]], unname(Mq[, "C"]))
  # Rotation actually changed the values (vs. the silent-drop bug where
  # scores would equal the un-rotated M).
  expect_gt(max(abs(as.matrix(scores) - M)), 0.1)
})
