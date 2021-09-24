all: Figure2.png Figure3/Figure3.png Figure4.png Figure1.png

Figure2A/Figure_2A.png: Figure2A/plot.R
	cd Figure2A/ && Rscript plot.R

Figure2B/Figure_2B.png: Figure2B/plot.R
	cd Figure2B/ && Rscript plot.R

Figure2.png:  Figure2A/Figure_2A.png  Figure2B/Figure_2B.png Figure2B/single-sample-happy-evaluation.csv
	convert -font helvetica -fill black -pointsize 150  -draw "gravity northwest text 0,0 'a'" Figure2A/Figure_2A.png tmp1.png
	convert -font helvetica -fill black -pointsize 150  -draw "gravity northwest text 0,0 'b'" Figure2B/Figure_2B.png tmp2.png
	convert +append tmp1.png tmp2.png $@
	rm tmp1.png tmp2.png

Figure3/Figure3.png: Figure3/plot.R
	cd Figure3/ && Rscript plot.R

Figure4.png: Figure4/plot.R
	cd Figure4/ && Rscript plot.R
	convert -font helvetica -fill black -pointsize 150  -draw "gravity northwest text 0,0 'a'" Figure4/minimac4-five-panel-bakeoff-hard.png 4_1.png
	convert -font helvetica -fill black -pointsize 150  -draw "gravity northwest text 0,0 'b'" Figure4/minimac4-five-panel-bakeoff.png 4_2.png
	convert -font helvetica -fill black -pointsize 150  -draw "gravity northwest text 0,0 'c'" Figure4/gtex-sensivity.png 4_3.png
	convert +append 4_1.png 4_2.png 4_3.png Figure4.png
	rm 4_1.png 4_2.png 4_3.png

Figure1.png: Figure_1A.png Figure1B/Figure_1B.png
	convert -font helvetica -fill black -pointsize 150  -draw "gravity northwest text 0,0 'a'" afam-african-proportions.png 1_a.png
	convert -font helvetica -fill black -pointsize 150  -draw "gravity northwest text 0,0 'b'" Figure1B/Figure_1B_v2.png 1_b.png
	convert -append 1_a.png 1_b.png $@
	rm 1_a.png 1_b.png
