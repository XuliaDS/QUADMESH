set term postscript enhanced eps solid color 
load '~/gnuplot-palettes/acdl2.pal' 
# stats 'finalMeshProperties_1.txt' u 5
# set cbtics add (STATS_min , STATS_max)
# set cbrange [STATS_min : STATS_max]
set cbrange [90:180]
set cbtics 30
set colorbox horizontal user origin .35, 0.22 size .3,.02
set colorbox  bdefault

unset border; unset xtic; unset ytics; unset ztics; unset key; set view 0,0; set colorbox noborder; 
set view equal; 
set view 120,90;

 set style line 10 lt rgb '#4d4d4d' lw 0.5 pt 6
set style line 11 lt rgb '#2c0710' lw 1 pt 6
 set format cb '%.f'; set cblabel font '18' 
 splot 'finalMeshProperties_2.txt' u 1:2:3:5 w l lw 3.5 palette; 
 replot 'gnuFin_2.txt' u 1:2:3 w l ls 10; 
 replot 'Edges_gnuFin_2.txt' u 1:2:3 w l ls 11;
unset colorbox
 set output 'colorplot.eps' ; replot; set output;
