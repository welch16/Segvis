# build package documentation
doc:
	R -e 'roxygen2::roxygenize()'

clean:
	rm -f *~
	rm -f */*~
	rm -f .*~
	rm -f inst/rscripts/*~
	rm -f vignettes/.Rhistory
	rm -f vignettes/.RData
	rm -f vignettes/vignette.aux
	rm -f vignettes/vignette.log
	rm -f vignettes/vignette.out
	rm -f vignettes/vignette.toc
	rm -fr vignettes/auto

# knit the vignettes
inst/doc/%.pdf:vignettes/%.Rnw
	cd vignettes;R CMD Sweave --engine=knitr::knitr --pdf $(<F);cp -f $(<F) ../inst/doc;mv -f $(<F:.Rnw=.pdf) ../inst/doc ;cd ..

# create the poster
poster/%.pdf:poster/%.tex poster/figs/fig0.pdf poster/figs/fig3.pdf poster/figs/fig2.pdf
	cd poster; pdflatex $(<F);pdflatex $(<F);pdflatex $(<F);cd ..

# figs poster
poster/figs/%.pdf:poster/figs/codes/%.R
	cd poster/figs/codes; R CMD BATCH $(<F);cd ../../..

