Starting point for TTbarZ analysis


--COMPILATION

To compile the code use MAKE
cd TTbarZ_analysis
make -f makefiles/makeTTZ clean; make -f makefiles/makeTTZ

--Plots/Tables/Datacards production for particular analysis



////////alternativeWay////////

cd alternativeWay
run CloneTree_bkg.C -o mc1.exe
./mc1.exe


for Ploting 

cd output
python plot1D.py -f test.list -v "lep_1stPt" -c " lep_1stPt > 100 " -t tree_4lep --xmin 0 --xmax  500 --nbins 50 --xlabel "lep_1stPt" -o tets.root --log


 
