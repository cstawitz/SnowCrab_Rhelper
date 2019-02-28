#'
#' @title Get ocean times from a set of ROMS netCDF files as a dataframe
#'
#' @description Function to get ocean times from a set of ROMS netCDF files as a dataframe
#'
#' @param topdir - top-level directory at which to start cataloging files/times
#' @param recursive - TRUE to work recursively through folder structure, FALSE to process topdir only
#' @param verbose - print debugging info
#'
#' @details ROMS files are assumed to end in ".nc". Uses the ncdf4 package library.
#'
#' @export
#'
getOceanTimes<-function(topdir=".",
                        recursive=TRUE,
                        verbose=TRUE){

  currdir<-getwd();
  on.exit(setwd(currdir));

  mdfr<-NULL;
  fldrs<-normalizePath(list.dirs(topdir,recursive=recursive));
  if (verbose) cat("--Folders to check:\n",paste("\t",fldrs,collapse="\n"),"\n\n")
  for (fldr in fldrs){
    cat("--Processing folder '",fldr,"'\n",sep="");
    setwd(fldr);
    fns<-list.files(pattern="^.*\\.nc$")
    for (fn in fns){
      ncf<-ncdf4::nc_open(fn)
      if ("ocean_time" %in% names(ncf$dim)){
        times<-ncf$dim$ocean_time$vals;
        t0<-gsub("seconds since ","",ncf$dim$ocean_time$units,fixed=TRUE);
        dfr<-data.frame(path=fldr,
                        name=fn,
                        ocean_time=times,
                        date=as.character(as.POSIXct(times,origin=t0,tz="GMT")),
                        stringsAsFactors=FALSE);
        mdfr<-rbind(mdfr,dfr);
      }
      nc_close(ncf);
    }
  }

  return(mdfr);
}
