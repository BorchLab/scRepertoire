# lifecycle throttles a deprecation warning to once per session per message id.
# .deprecate_arg() reuses the same id for every call to a given argument, so
# whichever test reaches it first would be the only one to see the warning and
# any later assertion on the same argument would fail depending on test order.
# "warning" makes every deprecation signal every time.
withr::local_options(lifecycle_verbosity = "warning", .local_envir = teardown_env())
