clear
close all
run AYfigprops.m

nbeads=3;
par_id=0;
MH_id=2;
test_id=2;

mhdoc = MH_doctor_data(nbeads,par_id,MH_id,test_id);
mhdoc.test_block = mhdoc.read_test_block; 
