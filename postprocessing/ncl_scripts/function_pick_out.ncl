procedure pick_out (ncases, cases, nsite,lats, lons,area, filepath)
; ncases, the number of models
; cases, the name of models
; casename, the name of cases
; filepath, model output filepath
; filepathobs, filepath for observational data

begin

 inptrs = new((/ncases/),  "file") ; pointers to files

 do im=0, ncases-1
    infile = filepath(im)+cases(im)+"/"+cases(im)+".cam.h0.1979-01.nc"
    print((/infile/))
    inptrs(im) = addfile(infile,"r")       ; pointer to file1

   inptr = inptrs(im)
   lat = inptr->lat
   nlat = dimsizes(lat)
   lon = inptr->lon
   nlon = dimsizes(lon)
   delete(inptr)
  print(nlat)
 idx_cols = new((/nsite,5/),  "integer")
 n =  new((/nsite/),  "integer")
 n(:)=0

   ; ========================================================================== 
   ; find out the cols and their numbers
   ; the axis of site is stored in idx_cols(site,n)
      do i = 0, nlat-1
       do j = 0, nsite-1 
        if (lon (i) .ge. lons(j)-area .and. lon (i) .lt. lons(j)+area \
        .and. lat (i) .ge. lats(j)-area .and. lat (i) .lt. lats(j)+area ) then 
         idx_cols(j,n(j))= i 
         n(j) = n(j)+1
        end if
       end do
      end do
   ; ==================================================================
     print (idx_cols (:,:))
     system ("mkdir ./data")
     system ("rm -f ./data/"+cases(im)+"_site_location.nc")
     ncdf     = addfile("./data/"+cases(im)+"_site_location.nc" ,"c")  ; open output netCDF file
     ncdf->n=n
     ncdf->ind_cols=idx_cols
 delete (lat)
 delete (lon)
 end do
end


