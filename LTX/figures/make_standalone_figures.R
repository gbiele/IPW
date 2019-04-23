

figurenames = c("SelectionBias","rg","estimates+raw", "decision_trees")
sfigurenames = c("covariation","IPW","logRRs","ROPE_plots")

library(tools)
library(magick)

for (f in c(figurenames)) {
  
  tmptex = readLines(paste0(f,".tex"))
  if (any(grepl("subfigure",tmptex))) {
    tmptex = gsub("4cm","\\\\textwidth",tmptex)
    tmptex = c("\\documentclass[a4]{article}",
               "\\usepackage[textwidth=5cm]{geometry}",
               "\\pagenumbering{gobble}",
               "\\usepackage{tikz}",
               "\\usepackage{subcaption}",
               "\\usepackage[labelformat=parens,labelsep=quad,skip=3pt]{caption}",
               "\\usetikzlibrary{arrows,automata,calc}",
               "\\begin{document}",
               "\\begin{figure}",
               tmptex,
               "\\end{figure}",
               "\\end{document}",
               "")
  } else {
    tmptex = c("\\documentclass[a4paper]{article}",
               "\\usepackage{tikz}",
               "\\pagenumbering{gobble}",
               "\\begin{document}",
               "\\begin{figure}",
               "\\centering",
               tmptex,
               "\\end{figure}",
               "\\end{document}",
               "")
  }
  
  cat(paste(tmptex,collpase = "\n"),
      file = "tmptex.tex")
  texi2dvi("tmptex.tex", pdf = T, clean = T)
   
  
  
  
  if (f %in% figurenames) {
    epsf = paste0("eps/Figure",which(f == figurenames),"_",f,".eps")
  } else {
    epsf = paste0("eps/SF",which(f == sfigurenames),"_",f,".eps")
  }
  img = image_read_pdf("tmptex.pdf", density = 600)
  image_write(image_trim(img),
              path = paste0("pngs/",f,".png"),
              format = "png",
              quality = 100)  
  
  system(paste("pdftops -eps tmptex.pdf",epsf) )
  
  file.remove(c("tmptex.pdf","tmptex.tex"))
}
