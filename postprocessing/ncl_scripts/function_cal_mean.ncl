procedure cal_mean (ncases, cases, nsite,lats, lons,area, filepath)
begin

 do im = 0, ncases -1
system("rm -f "+ filepath(im)+cases(im)+"/"+cases(im)+"_*_climo.nc")

; AMWG outputs
system("/usr/common/software/nco/4.7.4/intel/bin/ncea "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1980-01.nc" +"   "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1980-02.nc" + "   "\
 +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-12.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+"_DJF_climo.nc")
system("/usr/common/software/nco/4.7.4/intel/bin/ncea "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-03.nc" +"   "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-04.nc" + "   "\
 +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-05.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+"_MAM_climo.nc")
system("/usr/common/software/nco/4.7.4/intel/bin/ncea "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-06.nc" +"   "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-07.nc" + "   "\
 +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-08.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+"_JJA_climo.nc")
system("/usr/common/software/nco/4.7.4/intel/bin/ncea "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-09.nc" +"   "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-10.nc" + "   "\
 +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-11.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+"_SON_climo.nc")

system("/usr/common/software/nco/4.7.4/intel/bin/ncea "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-03.nc" +"   "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-04.nc" + "   "\
 +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-05.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-06.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-06.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-07.nc" \
+ "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-08.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-09.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-10.nc" + "   " \
+ "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-11.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-12.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1980-01.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1980-02.nc" + "   "  \
+filepath(im)+cases(im)+"/"+cases(im)+"_ANN_climo.nc")

; CLUBB budgets
system("/usr/common/software/nco/4.7.4/intel/bin/ncea "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1980-01.nc" +"   "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1980-02.nc" + "   "\
 +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-12.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+"_DJF_budget_climo.nc")
system("/usr/common/software/nco/4.7.4/intel/bin/ncea "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-03.nc" +"   "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-04.nc" + "   "\
 +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-05.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+"_MAM_budget_climo.nc")
system("/usr/common/software/nco/4.7.4/intel/bin/ncea "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-06.nc" +"   "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-07.nc" + "   "\
 +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-08.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+"_JJA_budget_climo.nc")
system("/usr/common/software/nco/4.7.4/intel/bin/ncea "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-09.nc" +"   "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-10.nc" + "   "\
 +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-11.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+"_SON_budget_climo.nc")

system("/usr/common/software/nco/4.7.4/intel/bin/ncea "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-03.nc" +"   "+filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-04.nc" + "   "\
 +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-05.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-06.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-06.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-07.nc" \
+ "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-08.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-09.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-10.nc" + "   " \
+ "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-11.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1979-12.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1980-01.nc" + "   " +filepath(im)+cases(im)+"/"+cases(im)+".cam.h2.1980-02.nc" + "   "  \
+filepath(im)+cases(im)+"/"+cases(im)+"_ANN_budget_climo.nc")
 end do


end


