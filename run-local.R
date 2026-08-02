options(shiny.maxRequestSize = 200 * 1024^2)
shiny::runApp(".", host = "127.0.0.1", port = 3838, launch.browser = FALSE)
