# This script renders multiple R Markdown files into HTML format and saves them in the specified output directory.
# The directory "./html" is triggers a GitHub workflow that deploys the HTML files to the project's GitHub Pages site.
# HTML generation is now being handed by calling rmarkdown::render_site() from directory md.

rmarkdown::render("./md/sankey.Rmd", output_file = "sankey.html", output_dir = "./html")
rmarkdown::render("./md/renda.Rmd", output_file = "renda.html", output_dir = "./html")
rmarkdown::render("./md/age_distribution.Rmd", output_file = "age_distribution.html", output_dir = "./html")
rmarkdown::render("./md/uncertainty.Rmd", output_file = "uncertainty.html", output_dir = "./html")
rmarkdown::render("./md/uncertainty.Rmd", output_file = "uncertainty.html", output_dir = "./html")
rmarkdown::render("./md/office_correlation.Rmd", output_file = "office_correlation.html", output_dir = "./html")