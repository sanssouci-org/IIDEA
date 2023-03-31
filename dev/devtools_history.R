devtools::document()
devtools::check()


usethis::use_package("shiny")
usethis::use_package("sanssouci")
usethis::use_package("sanssouci.data")

usethis::use_build_ignore("dev")
usethis::use_build_ignore("requirements.txt")

usethis::use_test("fct_import")
usethis::use_test("fct_output")
usethis::use_test("fct_post-hoc")

