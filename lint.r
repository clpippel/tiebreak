# lint.r

library(lintr)
setwd(paste0(Sys.getenv("R_USER"), "/Work"))

lint("gf.r",
  linters_with_defaults(
    object_name_linter = NULL,
    commented_code_linter = NULL,
    line_length_linter(132),
    semicolon_linter(allow_compound = TRUE),
    brace_linter(allow_single_line = TRUE)
  )
) -> ttt
ttt
